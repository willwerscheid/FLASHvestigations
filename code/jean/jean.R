library(flashier)
library(ggplot2)

mats <- readRDS("../../data/metabo3_gwas_mats.RDS")
mats$beta_hat[is.na(mats$se_hat)] <- NA

si.res <- softImpute::softImpute(mats$beta_hat, rank.max = 30, lambda = 0, type = "als")

base.args <- list(list(var.type = 0),
                  list(var.type = 1),
                  list(var.type = 2),
                  list(var.type = c(1, 2)),
                  list(S = mats$se_hat, var.type = NULL),
                  list(S = mats$se_hat, var.type = 0),
                  list(S = mats$se_hat, var.type = 1),
                  list(S = mats$se_hat, var.type = 2))
svd.args <- lapply(base.args, `c`, list(fit = "full"))
si.args <- lapply(base.args, `c`, list(fit = "full", init.fn = init.fn.softImpute))
affd.args <- lapply(base.args, `c`, list(fit = "backfit.only", EF.init = si.res))

all.args <- c(svd.args, si.args, affd.args)
all.args <- lapply(all.args, `c`, list(data = mats$beta_hat, output.lvl = 1))

flashier.res <- lapply(all.args, function(args) {
  cat(str(args))
  fl <- do.call(flashier, args)
  fl$call <- args
  return(fl)
})

saveRDS(flashier.res, "../../output/jean/flashier_res.rds")

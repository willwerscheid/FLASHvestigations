library(flashier)
library(ggplot2)

mats <- readRDS("../../data/metabo3_gwas_mats.RDS")
mats$beta_hat[is.na(mats$se_hat)] <- NA

si.res <- softImpute::softImpute(mats$beta_hat, rank.max = 30, lambda = 0, type = "als")

base.args <- list(list(var.type = 0),
                  list(var.type = 1),
                  list(var.type = 2),
                  list(var.type = c(1, 2)),
                  list(S = mats$se_hat, var.type = NULL))
svd.args <- lapply(base.args, `c`, list(fit = "full"))
si.args <- lapply(base.args, `c`, list(fit = "full", init.fn = init.fn.softImpute))
affd.args <- lapply(base.args, `c`, list(fit = "backfit.only", EF.init = si.res))

noisy.args <- list(list(S = mats$se_hat, var.type = 0),
                   list(S = mats$se_hat, var.type = 1),
                   list(S = mats$se_hat, var.type = 2))
noisy.svd.args <- lapply(noisy.args, `c`, list(fit = "full"))
noisy.si.args <- lapply(noisy.args, `c`, list(fit = "full", init.fn = init.fn.softImpute))
noisy.affd.args <- lapply(noisy.args, `c`, list(fit = "backfit.only", EF.init = si.res))

all.args <- c(svd.args, si.args, affd.args, noisy.svd.args, noisy.si.args, noisy.affd.args)
all.args <- lapply(all.args, `c`, list(data = mats$beta_hat, output.lvl = 1))

flashier.res <- lapply(all.args, function(args) {
  cat(str(args))
  fl <- do.call(flashier, args)
  fl$call <- args
  return(fl)
})

saveRDS(flashier.res, "../../output/jean/flashier_res.rds")

#
#
# # Some baselines.
# const.res <- flashier(jean$beta_hat, var.type = 0, fit = "full") # 16 factors; 3.016e+05
# byrow.res <- flashier(jean$beta_hat, var.type = 1, fit = "full") # 10 factors; 3.067e+05
# bycol.res <- flashier(jean$beta_hat, var.type = 2, fit = "full") # 8 factors; 3.135e+05
# kron.res  <- flashier(jean$beta_hat, var.type = c(1, 2), fit = "full") # 4 factors; 3.225e+05
# zero.res  <- flashier(jean$beta_hat, S = jean$se_hat, var.type = NULL, fit = "full") # 10 factors; 3.107e+05
#
# # Try softImpute.
# const.si <- flashier(jean$beta_hat, var.type = 0, fit = "full",
#                      init.fn = init.fn.softImpute) # 19 factors; 3.058e+05
# byrow.si <- flashier(jean$beta_hat, var.type = 1, fit = "full",
#                      init.fn = init.fn.softImpute) # 10 factors; 3.065e+05
# bycol.si <- flashier(jean$beta_hat, var.type = 2, fit = "full",
#                      init.fn = init.fn.softImpute) # 13 factors; 3.177e+05
# kron.si  <- flashier(jean$beta_hat, var.type = c(1, 2), fit = "full",
#                      init.fn = init.fn.softImpute) # 7 factors; 3.250e+05
# zero.si  <- flashier(jean$beta_hat, S = jean$se_hat, var.type = NULL, fit = "full",
#                      init.fn = init.fn.softImpute) # 14 factors; 3.172e+05
#
# # Check against flashr.
# flashr.const <- flashr::flash(jean$beta_hat, var_type = "constant",
#                               backfit = TRUE) # 21 factors, 3.084e+05; takes forever!
# flashr.byrow <- flashr::flash(jean$beta_hat, var_type = "by_row",
#                               backfit = TRUE) # 10 factors, 3.065e+05
# flashr.bycol <- flashr::flash(jean$beta_hat, var_type = "by_col",
#                               backfit = TRUE) # 10 factors; 3.153e+05
# flash.data   <- flashr::flash_set_data(jean$beta_hat, jean$se_hat)
# flashr.zero  <- flashr::flash(flash.data, var_type = "zero",
#                               backfit = TRUE) # 14 factors; 3.173e+05
#
# # Try the add_factors_from_data method.
# flashr.affd <- flashr::flash_add_factors_from_data(flash.data, K = 30, var_type = "zero",
#                                                    backfit = TRUE) # 26 factors; 3.196e+05
#
# si.init <- softImpute::softImpute(jean$beta_hat, rank.max = 30, lambda = 0, type = "als")
#
# const.affd <- flashier(jean$beta_hat, var.type = 0, EF.init = si.init,
#                        fit = "backfit.only") # 30 factors; 3.143e+05
# bycol.affd <- flashier(jean$beta_hat, var.type = 2, EF.init = si.init,
#                        fit = "backfit.only") # 30 factors; 3.203e+05
# kron.affd  <- flashier(jean$beta_hat, var.type = c(1, 2), EF.init = si.init,
#                        fit = "backfit.only") # 29 factors; 3.273e+05; slow to converge
# zero.affd  <- flashier(jean$beta_hat, S = jean$se_hat, var.type = NULL, EF.init = si.init,
#                        fit = "backfit.only", maxiter = 10) # 30 factors; 3.191e+05
#
# # Useful to compare estimated variances against supplied ones.
# var.df <- data.frame(est.se = rep(bycol.affd$residuals.sd,
#                                   each = nrow(jean$beta_hat))[!is.na(jean$beta_hat)],
#                      given.se = jean$se_hat[!is.na(jean$beta_hat)])
# ggplot(var.df, aes(x = given.se, y = est.se)) + geom_point(alpha = 0.2) +
#   geom_abline(slope = 1, linetype = "dashed")
#
# # Try "noisy" variance structures.
# const.noisy <- flashier(jean$beta_hat, S = jean$se_hat, var.type = 0,
#                         fit = "full") # 9 factors; 3.115e+05
# bycol.noisy <- flashier(jean$beta_hat, S = jean$se_hat, var.type = 1,
#                         fit = "full") # 6 factors; 3.136e+05
# const.noisy.si <- flashier(jean$beta_hat, S = jean$se_hat, var.type = 0,
#                            fit = "full", init.fn = init.fn.softImpute) # 9 factors; 3.136e+05
# bycol.noisy.si <- flashier(jean$beta_hat, S = jean$se_hat, var.type = 1,
#                            fit = "full", init.fn = init.fn.softImpute) # 6 factors; 3.145e+05
# const.noisy.affd <- flashier(jean$beta_hat, S = jean$se_hat, var.type = 0,
#                              EF.init = si.init, fit = "backfit.only") # 29 factors; 3.194e+05
# bycol.noisy.affd <- flashier(jean$beta_hat, S = jean$se_hat, var.type = 1,
#                              EF.init = si.init, fit = "backfit.only") # 29 factors; 3.212e+05
# byrow.noisy.affd <- flashier(jean$beta_hat, S = jean$se_hat, var.type = 2,
#                              EF.init = si.init, fit = "backfit.only") # 29 factors; 3.192e+05
#
# # "noisy kron" doesn't work with missing data.
# # kron.noisy <- flashier(jean$beta_hat, S = jean$se_hat, var.type = c(1, 2),
# #                        EF.init = EF.init, fit = "backfit.only")

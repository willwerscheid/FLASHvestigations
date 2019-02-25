devtools::load_all("~/Github/ashr")
devtools::load_all("~/Github/flashier")
library(Matrix)

# Anscombe transform --------------------------------------------------------

trachea <- read.table("~/Downloads/GSE103354_Trachea_droplet_UMIcounts.txt")
trachea <- as.matrix(trachea)
# trachea <- Matrix(trachea)
trachea <- 2 * sqrt(trachea + 0.375)

# ncells <- rowSums(trachea > 0)
# gene.idx <- which(ncells > 4)
# trachea <- trachea[gene.idx, ]

# saveRDS(trachea, "./data/tmpdata.rds")
trachea <- readRDS("./data/tmpdata.rds")

# Let's flash.
flash.fit <- flashier(trachea,
                      S = 1, var.type = NULL,
                      prior.type = c("normal.mix", "nonnegative"),
                      greedy.Kmax = 10,
                      backfit = "alternating", backfit.order = "dropout",
                      backfit.maxiter = 5, backfit.reltol = 10,
                      verbose = 3)
# saveRDS(flash.fit, "./data/tmp.rds")

flash.fit2 <- flashier(trachea, flash.init = flash.fit,
                       S = 1, var.type = NULL,
                       prior.type = c("normal.mix", "nonnegative"),
                       greedy.Kmax = 25,
                       backfit = "alternating", backfit.order = "dropout",
                       backfit.maxiter = 3, backfit.reltol = 10,
                       verbose = 3)
# saveRDS(flash.fit2, "./data/tmp2.rds")
# rm(flash.fit)

flash.fit3 <- flashier(trachea, flash.init = flash.fit2,
                       S = 1, var.type = NULL,
                       prior.type = c("normal.mix", "nonnegative"),
                       greedy.Kmax = 50,
                       backfit.every = 3, backfit.order = "dropout",
                       backfit.maxiter = 5, backfit.reltol = 10,
                       verbose = 3)
# saveRDS(flash.fit3, "./data/tmp3.rds")
# rm(flash.fit2)

flash.fit4 <- flashier(trachea, flash.init = flash.fit3,
                       S = 1, var.type = NULL,
                       prior.type = c("normal.mix", "nonnegative"),
                       greedy.Kmax = 50,
                       final.backfit = TRUE, backfit.order = "dropout",
                       backfit.maxiter = 100, backfit.reltol = 10,
                       verbose = 3)
# saveRDS(flash.fit4, "./data/Anscombe50.rds")
# rm(flash.fit3)

n.rounds <- 25
samps.per.round <- 20

fl <- flash.fit4
gene.signs <- sign(fl$loadings$normalized.loadings[[1]])
lfsr <- array(0, dim = dim(gene.signs))
set.seed(666)
for (i in 1:n.rounds) {
  cat(paste("Sampling round", i, "\n"))
  samp <- fl$sampler(samps.per.round)
  samp <- lapply(samp, `[[`, 1)
  samp <- lapply(samp, sign)
  samp <- lapply(samp, function(x) x == gene.signs)
  lfsr <- lfsr + Reduce(`+`, samp)
}
lfsr <- lfsr / (n.rounds * samps.per.round)
# saveRDS(lfsr, "./data/AnscombeLFSR.rds")

# log1p transform -----------------------------------------------------------

trachea <- read.table("~/Downloads/GSE103354_Trachea_droplet_UMIcounts.txt")
trachea <- as.matrix(trachea)
trachea <- log1p(trachea)

# saveRDS(trachea, "./data/tmpdata.rds")
# trachea <- readRDS("./data/tmpdata.rds")

S <- sqrt(trachea) / (trachea + 1)
S[S == 0] <- 0.5

trachea <- set.flash.data(data = trachea, S = S, var.type = NULL)

rm(S)

flash.fit <- flashier(trachea, var.type = NULL,
                      prior.type = c("normal.mix", "nonnegative"),
                      greedy.Kmax = 10,
                      backfit = "alternating", backfit.order = "dropout",
                      backfit.maxiter = 5, backfit.reltol = 10,
                      verbose = 3)
# saveRDS(flash.fit, "./tmp.rds")

flash.fit2 <- flashier(trachea, flash.init = flash.fit,
                       var.type = NULL,
                       prior.type = c("normal.mix", "nonnegative"),
                       greedy.Kmax = 25,
                       backfit = "alternating", backfit.order = "dropout",
                       backfit.maxiter = 3, backfit.reltol = 10,
                       verbose = 3)
# saveRDS(flash.fit2, "./tmp2.rds")
# rm(flash.fit)

flash.fit3 <- flashier(trachea, flash.init = flash.fit2,
                       var.type = NULL,
                       prior.type = c("normal.mix", "nonnegative"),
                       greedy.Kmax = 50,
                       backfit.every = 3, backfit.order = "dropout",
                       backfit.maxiter = 5, backfit.reltol = 10,
                       verbose = 3)
saveRDS(flash.fit3, "./log1p39.rds")

fl <- flash.fit3
gene.signs <- sign(fl$loadings$normalized.loadings[[1]])
lfsr <- array(0, dim = dim(gene.signs))
set.seed(666)
for (i in 1:n.rounds) {
  cat(paste("Sampling round", i, "\n"))
  samp <- fl$sampler(samps.per.round)
  samp <- lapply(samp, `[[`, 1)
  samp <- lapply(samp, sign)
  samp <- lapply(samp, function(x) x == gene.signs)
  lfsr <- lfsr + Reduce(`+`, samp)
}
lfsr <- lfsr / (n.rounds * samps.per.round)
# saveRDS(lfsr, "./data/log1pLFSR.rds")

# "pseudodata" approach -----------------------------------------------------

Y <- read.table("~/Downloads/GSE103354_Trachea_droplet_UMIcounts.txt")
Y <- as.matrix(Y)
log.pm <- flashier:::lowrank.expand(flash.fit3$fit$EF)
pm <- exp(log.pm)
X <- log.pm + (Y - pm) / pm
S <- 1 / sqrt(pm)

rm(Y)
rm(log.pm)
rm(pm)

pseudodata <- set.flash.data(data = X, S = S, var.type = NULL)
rm(X)
rm(S)

flash.fit <- flashier(pseudodata, var.type = NULL,
                      prior.type = c("normal.mix", "nonnegative"),
                      greedy.Kmax = 10,
                      backfit = "alternating", backfit.order = "dropout",
                      backfit.maxiter = 5, backfit.reltol = 10,
                      verbose = 3)
# saveRDS(flash.fit, "./tmp.rds")

flash.fit2 <- flashier(pseudodata, flash.init = flash.fit,
                       prior.type = c("normal.mix", "nonnegative"),
                       greedy.Kmax = 25,
                       backfit = "alternating", backfit.order = "dropout",
                       backfit.maxiter = 3, backfit.reltol = 10,
                       verbose = 3)
# saveRDS(flash.fit2, "./tmp2.rds")
# rm(flash.fit)

flash.fit3 <- flashier(pseudodata, flash.init = flash.fit2,
                       prior.type = c("normal.mix", "nonnegative"),
                       greedy.Kmax = 50,
                       backfit.every = 3, backfit.order = "dropout",
                       backfit.maxiter = 5, backfit.reltol = 10,
                       verbose = 3)
# saveRDS(flash.fit3, "./data/tmp3.rds")
# rm(flash.fit2)

flash.fit4 <- flashier(pseudodata, flash.init = flash.fit3,
                       prior.type = c("normal.mix", "nonnegative"),
                       greedy.Kmax = 50,
                       final.backfit = TRUE, backfit.order = "dropout",
                       backfit.maxiter = 100, backfit.reltol = 10,
                       verbose = 3)
# saveRDS(flash.fit4, "./data/pseudo50.rds")
# rm(flash.fit3)

fl <- flash.fit4
gene.signs <- sign(fl$loadings$normalized.loadings[[1]])
lfsr <- array(0, dim = dim(gene.signs))
set.seed(666)
for (i in 1:n.rounds) {
  cat(paste("Sampling round", i, "\n"))
  samp <- fl$sampler(samps.per.round)
  samp <- lapply(samp, `[[`, 1)
  samp <- lapply(samp, sign)
  samp <- lapply(samp, function(x) x == gene.signs)
  lfsr <- lfsr + Reduce(`+`, samp)
}
lfsr <- lfsr / (n.rounds * samps.per.round)
# saveRDS(lfsr, "./data/pseudoLFSR.rds")

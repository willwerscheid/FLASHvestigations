library(flashier)

greedy.K <- 50
backfit.maxiter <- 500
n.cores <- parallel::detectCores()
data <- readRDS("../data/flashier_bench/sp_trachea.rds")

fl.g <- flashier(data,
                 prior.family = prior.normal.scale.mix(),
                 greedy.Kmax = greedy.K,
                 backfit = "none")

# Usual "dropout" backfit:
sink("../output/parallel_v2/dropout_res.txt")
# system.time() gives weird results for parallel fits, so use Sys.time() instead.
t0 <- Sys.time()
fl.dropout <- flashier(flash.init = fl.g,
                       backfit = "only",
                       backfit.maxiter = backfit.maxiter,
                       verbose.lvl = -1)
t.dropout <- Sys.time() - t0

# Parallel backfit:
sink("../output/parallel_v2/parallel_res.txt")
options(cl.type = "FORK")
options(cl.cores = n.cores)
t0 <- Sys.time()
fl.parallel <- flashier(flash.init = fl.g,
                        backfit = "only",
                        backfit.maxiter = backfit.maxiter,
                        backfit.order = "parallel",
                        verbose.lvl = -1)
t.parallel <- Sys.time() - t0

sink()
all.t <- list(dropout = t.dropout, parallel = t.parallel)
saveRDS(all.t, "../output/parallel_v2/timing.rds")

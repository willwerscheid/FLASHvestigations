library(flashier)

greedy.K <- 50
backfit.maxiter <- 500
data <- readRDS(paste0("../../data/flashier_bench/", data.name, ".rds"))

fl.g <- flashier(data,
                 prior.family = prior.normal.scale.mix(),
                 greedy.Kmax = greedy.K,
                 fit = "greedy")

# Usual backfit:
sink(paste0("../../output/parallel_v2/", data.name, "_serial.txt"))
# system.time() gives weird results for parallel fits, so use Sys.time() instead.
t0 <- Sys.time()
fl.serial <- flashier(init = fl.g,
                      fit = "backfit.only",
                      backfit.maxiter = backfit.maxiter,
                      verbose.lvl = -1)
t.serial <- Sys.time() - t0

# Parallel backfit:
sink(paste0("../../output/parallel_v2/", data.name, "_parallel.txt"))
options(cl.type = "FORK", cl.cores = n.cores)
t0 <- Sys.time()
fl.parallel <- flashier(init = fl.g,
                        fit = "backfit.only",
                        backfit.maxiter = backfit.maxiter,
                        backfit.reltol = parallel.backfit.reltol,
                        backfit.order = "parallel",
                        verbose.lvl = -1)
t.parallel <- Sys.time() - t0

sink()
all.t <- list(serial = t.serial, parallel = t.parallel)
saveRDS(all.t, paste0("../../output/parallel_v2/", data.name, "_timing.rds"))

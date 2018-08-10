niter <- 40
nfactors <- 20

# devtools::install_github("stephenslab/flashr")
devtools::load_all("/Users/willwerscheid/GitHub/flashr/")
# devtools::install_github("stephenslab/ebnm")
devtools::load_all("/Users/willwerscheid/GitHub/ebnm/")

source("./code/parallel.R")

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))

strong <- t(gtex$strong.z)
strong_data <- flash_set_data(strong, S = 1)

fl <- flash_add_greedy(strong_data, nfactors, var_type = "zero")
fl_greedy <- fl

fl <- fl_greedy # load
backfit_t <- rep(0, niter)
backfit_obj <- rep(0, niter)
for (i in 1:niter) {
   t <- system.time(
    for (k in 1:nfactors) {
      fl <- flashr:::flash_update_single_fl(strong_data,
                                            fl,
                                            k,
                                            "zero",
                                            "ebnm_pn", list(),
                                            "ebnm_pn", list())
    }
  )
  backfit_t[i] <- t[3] # elapsed time
  backfit_obj[i] <- flash_get_objective(strong_data, fl)
}

fl <- fl_greedy
parallel_t <- rep(0, niter)
parallel_obj <- rep(0, niter)
ebnm_param <- vector("list", nfactors)
for (k in 1:nfactors) {ebnm_param[[k]] <- list()}
for (i in 1:niter) {
  t <- system.time(
    fl <- flash_update_fl_parallel(strong_data,
                                 fl,
                                 1:nfactors,
                                 "zero",
                                 as.list(rep("ebnm_pn", nfactors)),
                                 ebnm_param,
                                 as.list(rep("ebnm_pn", nfactors)),
                                 ebnm_param)
  )
  parallel_t[i] <- t[3]
  parallel_obj[i] <- flash_get_objective(strong_data, fl)
}

res <- list(backfit_t = backfit_t,
            parallel_t = parallel_t,
            backfit_obj = backfit_obj,
            parallel_obj = parallel_obj)
saveRDS(res, "./data/parallel/nfactor20niter40.rds")

plot(backfit_obj, pch=1)
points(parallel_obj, pch=2)
plot(backfit_obj - parallel_obj, type = "l")

boxplot(cbind(backfit_t, parallel_t))
c(sum(backfit_t), sum(parallel_t))

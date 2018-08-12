# devtools::install_github("stephenslab/flashr")
devtools::load_all("/Users/willwerscheid/GitHub/flashr/")
# devtools::install_github("stephenslab/ebnm")
devtools::load_all("/Users/willwerscheid/GitHub/ebnm/")

library(parallel)
source("./code/parallel.R")

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))

strong <- t(gtex$strong.z)
strong_data <- flash_set_data(strong, S = 1)

run_test <- function(fl_init, niter) {
  nfactors <- flash_get_k(fl_init)

  message("Usual backfit...")
  fl <- fl_init
  backfit_t <- rep(0, niter)
  backfit_obj <- rep(0, niter)
  for (i in 1:niter) {
    message("  Iteration ", i)
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

  message("Parallel updates with lapply...")
  fl <- fl_init
  parallel_t <- rep(0, niter)
  parallel_obj <- rep(0, niter)
  ebnm_param <- vector("list", nfactors)
  for (k in 1:nfactors) {ebnm_param[[k]] <- list()}
  for (i in 1:niter) {
    message("  Iteration ", i)
    t <- system.time(
      fl <- flash_update_fl_parallel(strong_data,
                                     fl,
                                     1:nfactors,
                                     "zero",
                                     as.list(rep("ebnm_pn", nfactors)),
                                     ebnm_param,
                                     as.list(rep("ebnm_pn", nfactors)),
                                     ebnm_param,
                                     lapply)
    )
    parallel_t[i] <- t[3]
    parallel_obj[i] <- flash_get_objective(strong_data, fl)
  }

  message("Parallel updates with mclapply...")
  fl <- fl_init
  multicore_t <- rep(0, niter)
  for (i in 1:niter) {
    message("  Iteration ", i)
    t <- system.time(
      fl <- flash_update_fl_parallel(strong_data,
                                     fl,
                                     1:nfactors,
                                     "zero",
                                     as.list(rep("ebnm_pn", nfactors)),
                                     ebnm_param,
                                     as.list(rep("ebnm_pn", nfactors)),
                                     ebnm_param,
                                     mclapply)
    )
    multicore_t[i] <- t[3]
  }

  res <- list(backfit_t = backfit_t,
              parallel_t = parallel_t,
              multicore_t = multicore_t,
              backfit_obj = backfit_obj,
              parallel_obj = parallel_obj)
}

fl_greedy <- flash_add_greedy(strong_data, 20, var_type = "zero")
res_greedy <- run_test(fl_greedy, 100)
saveRDS(res_greedy, "./data/parallel/greedy20niter100.rds")

fl_svd <- flash_add_factors_from_data(strong_data, 20, init_fn = "udv_svd")
res_svd <- run_test(fl_svd, 100)
saveRDS(res_svd, "./data/parallel/svd20niter100.rds")

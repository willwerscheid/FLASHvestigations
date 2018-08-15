# devtools::install_github("stephenslab/flashr")
devtools::load_all("/Users/willwerscheid/GitHub/flashr/")
# devtools::install_github("stephenslab/ebnm")
devtools::load_all("/Users/willwerscheid/GitHub/ebnm/")

library(parallel)
source("./code/parallel.R")

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))

strong <- t(gtex$strong.z)
strong_data <- flash_set_data(strong, S = 1)

run_test <- function(data, fl_init, niter,
                     ebnm_fn_l = "ebnm_pn", ebnm_fn_f = "ebnm_pn",
                     ebnm_param_l = NULL, ebnm_param_f = NULL,
                     ksets = NULL) {
  nfactors <- flash_get_k(fl_init)

  if (is.null(ebnm_param_l)) {
    ebnm_param_l <- vector("list", nfactors)
    for (k in 1:nfactors) {ebnm_param_l[[k]] <- list(warmstart = TRUE)}
  }

  if (is.null(ebnm_param_f)) {
    ebnm_param_f <- vector("list", nfactors)
    for (k in 1:nfactors) {ebnm_param_f[[k]] <- list(warmstart = TRUE)}
  }

  if (is.null(ksets)) {
    ksets = list()
    ksets[[1]] = 1:nfactors
  }

  message("Usual backfit...")
  fl <- fl_init
  backfit_t <- rep(0, niter)
  backfit_obj <- rep(0, niter)
  for (i in 1:niter) {
    message("  Iteration ", i)
    t <- system.time(
      for (k in 1:nfactors) {
        fl <- flashr:::flash_update_single_fl(data,
                                              fl,
                                              k,
                                              "zero",
                                              ebnm_fn_l,
                                              ebnm_param_l[[k]],
                                              ebnm_fn_f,
                                              ebnm_param_f[[k]])
      }
    )
    backfit_t[i] <- t[3] # elapsed time
    backfit_obj[i] <- flash_get_objective(data, fl)
  }

  message("Parallel updates with lapply...")
  fl <- fl_init
  parallel_t <- rep(0, niter)
  parallel_obj <- rep(0, niter)
  for (i in 1:niter) {
    message("  Iteration ", i)
    t <- system.time({
      for (kset in ksets) {
        fl <- flash_update_fl_parallel(data,
                                       fl,
                                       kset,
                                       "zero",
                                       as.list(rep(ebnm_fn_l, nfactors)),
                                       ebnm_param_l,
                                       as.list(rep(ebnm_fn_f, nfactors)),
                                       ebnm_param_f,
                                       lapply)
      }
    })
    parallel_t[i] <- t[3]
    parallel_obj[i] <- flash_get_objective(data, fl)
  }

  message("Parallel updates with mclapply...")
  fl <- fl_init
  multicore_t <- rep(0, niter)
  for (i in 1:niter) {
    message("  Iteration ", i)
    t <- system.time({
      for (kset in ksets) {
        fl <- flash_update_fl_parallel(data,
                                       fl,
                                       kset,
                                       "zero",
                                       as.list(rep(ebnm_fn_l, nfactors)),
                                       ebnm_param_l,
                                       as.list(rep(ebnm_fn_f, nfactors)),
                                       ebnm_param_f,
                                       mclapply)
      }
    })
    multicore_t[i] <- t[3]
  }

  res <- list(backfit_t = backfit_t,
              parallel_t = parallel_t,
              multicore_t = multicore_t,
              backfit_obj = backfit_obj,
              parallel_obj = parallel_obj)
}

fl_greedy <- flash_add_greedy(strong_data, 20, var_type = "zero")
res_greedy <- run_test(strong_data, fl_greedy, 100)
saveRDS(res_greedy, "./data/parallel/greedy20niter100.rds")

fl_svd <- flash_add_factors_from_data(strong_data, 20, init_fn = "udv_svd")
res_svd <- run_test(strong_data, fl_svd, 100)
saveRDS(res_svd, "./data/parallel/svd20niter100.rds")


## Test MASH v FLASH backfits:

random <- t(gtex$random.z)
random_data <- flash_set_data(random, S = 1)

fpath <- "/Users/willwerscheid/GitHub/MASHvFLASH/output/"
nn <- readRDS(paste0(fpath, "MASHvFLASHnn/fl.rds"))
multi <- c(2, 5, 6, 8, 11:13, 17, 22:25, 31)
n <- nrow(strong)
dd <- nn$EL[, multi]
dd <- dd / rep(apply(dd, 2, max), each=n) # normalize
canonical <- cbind(rep(1, n), diag(rep(1, n)))
LL <- cbind(canonical, dd)

fl_random <- flash_add_fixed_loadings(random_data, LL)
res_random <- run_test(random_data, fl_random, 5)
saveRDS(res_random, "./data/parallel/MASHvFLASHrandom_bad.rds")

ksets=list(1, 2:45, c(46, 50), c(47:49, 51:flash_get_k(fl_final)))
res_random <- run_test(random_data, fl_random, 20,
                       ksets = ksets)
saveRDS(res_random, "./data/parallel/MASHvFLASHrandom.rds")

res_final <- run_test(strong_data, fl_final, 20,
                      ebnm_param_f = ebnm_param_f,
                      ksets = ksets)
saveRDS(res_final, "./data/parallel/MASHvFLASHfinal.rds")



ksets=list(1, 2:45, c(46, 50), c(47:49, 51:flash_get_k(fl_final)))
res_random <- run_test(random_data, fl_random, 20,
                       ksets=list(c(1, 46:flash_get_k(fl_random)), 2:45))
saveRDS(res_random, "./data/parallel/MASHvFLASHrandom.rds")

fl_final <- flash_add_fixed_loadings(strong_data, LL)
gf <- readRDS(paste0(fpath, "MASHvFLASHgtex3/flgf.rds"))
ebnm_param_f = lapply(gf, function(g) {list(g=g, fixg=TRUE)})
res_final <- run_test(strong_data, fl_final, 5,
                      ebnm_param_f = ebnm_param_f, ksets=ksets)
saveRDS(res_final, "./data/parallel/MASHvFLASHfinal_bad.rds")

ksets=list(1, 2:45, c(46, 50), c(47:49, 51:flash_get_k(fl_final)))

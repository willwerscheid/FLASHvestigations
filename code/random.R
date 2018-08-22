devtools::load_all("/Users/willwerscheid/GitHub/flashr/") # use "dev" branch
devtools::load_all("/Users/willwerscheid/GitHub/ebnm/")

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))
strong <- t(gtex$strong.z)
strong_data <- flash_set_data(strong, S = 1)

fl_init <- readRDS("/Users/willwerscheid/GitHub/MASHvFLASH/output/MASHvFLASHnn/fl_g.rds")


# Sequential backfit (from factor/loading 1 to factor/loading 34) -------
#
ebnm_param = list(f = list(), l = list(mixcompdist="+uniform"))
fl1 <- flash_backfit_workhorse(strong_data,
                               fl_init,
                               kset = 1:34,
                               var_type = "zero",
                               ebnm_fn = "ebnm_ash",
                               ebnm_param = ebnm_param,
                               verbose_output = "odLn",
                               nullcheck = FALSE)

# Use warmstarts to clamp it down:
ebnm_param = list(f = list(warmstart = TRUE),
                  l = list(mixcompdist = "+uniform", warmstart = TRUE))
fl1 <- flash_backfit_workhorse(strong_data,
                               fl1,
                               kset = 1:34,
                               var_type = "zero",
                               ebnm_fn = "ebnm_ash",
                               ebnm_param = ebnm_param,
                               verbose_output = "odLn",
                               nullcheck = FALSE)
## 77 + 375 iterations; -1257537
saveRDS(fl1, "./data/random/fl1.rds")


# Sequential backfit in reverse order (from 34 to 1) --------------------
#
ebnm_param = list(f = list(), l = list(mixcompdist="+uniform"))
fl2 <- flash_backfit_workhorse(strong_data,
                               fl_init,
                               kset = 34:1,
                               var_type = "zero",
                               ebnm_fn = "ebnm_ash",
                               ebnm_param = ebnm_param,
                               verbose_output = "odLn",
                               nullcheck = FALSE)

ebnm_param = list(f = list(warmstart = TRUE),
                  l = list(mixcompdist = "+uniform", warmstart = TRUE))
fl2 <- flash_backfit_workhorse(strong_data,
                               fl2,
                               kset = 34:1,
                               var_type = "zero",
                               ebnm_fn = "ebnm_ash",
                               ebnm_param = ebnm_param,
                               verbose_output = "odLn",
                               nullcheck = FALSE)
## 77 + 125 iterations; objective: -1257237
saveRDS(fl2, "./data/random/fl2.rds")


# Backfit the factor/loadings in a random order -------------------------
#
set.seed(666)
ebnm_param = list(f = list(), l = list(mixcompdist="+uniform"))
fl3 <- fl_init
old_obj <- flash_get_objective(strong_data, fl3)
diff <- Inf
iter <- 0
while (diff > .01) {
  iter <- iter + 1
  fl3 <- flash_backfit_workhorse(strong_data,
                                 fl3,
                                 kset=sample(1:34, 34),
                                 var_type = "zero",
                                 ebnm_fn = "ebnm_ash",
                                 ebnm_param = ebnm_param,
                                 verbose_output = "",
                                 maxiter = 1,
                                 nullcheck = FALSE)
  obj <- flash_get_objective(strong_data, fl3)
  message("Iteration ", iter, ": ", obj)
  diff <- obj - old_obj
  old_obj <- obj
}

ebnm_param = list(f = list(warmstart = TRUE),
                  l = list(mixcompdist = "+uniform", warmstart = TRUE))
old_obj <- flash_get_objective(strong_data, fl3)
diff <- Inf
iter <- 0
while (diff > .01) {
  iter <- iter + 1
  fl3 <- flash_backfit_workhorse(strong_data,
                                 fl3,
                                 kset=sample(1:34, 34),
                                 var_type = "zero",
                                 ebnm_fn = "ebnm_ash",
                                 ebnm_param = ebnm_param,
                                 verbose_output = "",
                                 maxiter = 1,
                                 nullcheck = FALSE)
  obj <- flash_get_objective(strong_data, fl3)
  message("Iteration ", iter, ": ", obj)
  diff <- obj - old_obj
  old_obj <- obj
}
# 52 + 231 iterations; final obj: -1257408
saveRDS(fl3, "./data/random/fl3.rds")

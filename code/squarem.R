# Load packages ---------------------------------------------------------

# devtools::install_github("stephenslab/flashr")
devtools::load_all("/Users/willwerscheid/GitHub/flashr/")
# devtools::install_github("stephenslab/ebnm")
devtools::load_all("/Users/willwerscheid/GitHub/ebnm/")
library(SQUAREM)
library(daarem)

# Load data -------------------------------------------------------------

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))
data <- flash_set_data(t(gtex$strong.z), S = 1)

# Create flash objects for backfitting ----------------------------------

fl_g5 <- flash_add_greedy(data, Kmax = 5,
                          var_type = "zero", init_fn = "udv_svd")

fl_g10 <- flash_add_greedy(data, Kmax = 5, f_init = fl_g5,
                           var_type = "zero", init_fn = "udv_svd")

fl_g20 <- flash_add_greedy(data, Kmax = 10, f_init = fl_g10,
                           var_type = "zero", init_fn = "udv_svd")

# Testing functions -----------------------------------------------------

logit <- function(x) {log(x / (1 - x))}
inv_logit <- function(x) {exp(x) / (1 + exp(x))}

fl.to.param <- function(f) {
  c(as.vector(f$EL), as.vector(f$EF),
    as.vector(f$EL2), as.vector(f$EF2),
    sapply(f$gl, function(k) {logit(k$pi0)}),
    sapply(f$gf, function(k) {logit(k$pi0)}),
    sapply(f$gl, function(k) {log(k$a)}),
    sapply(f$gf, function(k) {log(k$a)}))
}

param.to.fl <- function(param, m, n, k) {
  LL = matrix(param[1:(m*k)], ncol=k)
  curr_idx = m*k
  FF = matrix(param[(curr_idx + 1):(curr_idx + n*k)], ncol=k)
  curr_idx = curr_idx + n*k

  f = flashr:::flash_init_lf(LL, FF)

  f$EL2 = matrix(param[(curr_idx + 1):(curr_idx + m*k)], ncol=k)
  curr_idx = curr_idx + m*k
  f$EF2 = matrix(param[(curr_idx + 1):(curr_idx + n*k)], ncol=k)
  curr_idx = curr_idx + n*k

  f$gl = list()
  f$gf = list()
  for (i in 1:k) {
    f$gl[[i]] = list()
    f$gl[[i]]$pi0 = inv_logit(param[curr_idx + 1])
    curr_idx = curr_idx + 1
  }
  for (i in 1:k) {
    f$gf[[i]] = list()
    f$gf[[i]]$pi0 = inv_logit(param[curr_idx + 1])
    curr_idx = curr_idx + 1
  }
  for (i in 1:k) {
    f$gl[[i]]$a = exp(param[curr_idx + 1])
    curr_idx = curr_idx + 1
  }
  for (i in 1:k) {
    f$gf[[i]]$a = exp(param[curr_idx + 1])
    curr_idx = curr_idx + 1
  }

  return(f)
}

flash.iter <- function(p, data, m, n, k) {
  init_fl = param.to.fl(p, m, n, k)
  fl = flash_backfit(data, init_fl,
                     ebnm_fn = "ebnm_pn", var_type = "zero",
                     nullcheck = FALSE, verbose = FALSE,
                     maxiter = 1)
  obj = flash_get_objective(data, fl)
  message(obj)
  return(c(fl.to.param(fl), obj))
}

flash.obj <- function(p, data, m, n, k) {
  return(p[length(p)])
}

run_test <- function(f_init, data, niter) {
  m <- flash_get_n(f_init)
  n <- flash_get_p(f_init)
  k <- flash_get_k(f_init)

  backfit_obj <- rep(NA, niter)
  f <- f_init
  backfit_t <- system.time({
    for (i in 1:niter) {
      f <- flash_backfit(data, f,
                         ebnm_fn = "ebnm_pn", var_type = "zero",
                         nullcheck = FALSE, verbose = FALSE,
                         maxiter = 1)
      obj <- flash_get_objective(data, f)
      message(obj)
      backfit_obj[i] <- obj
    }
  })

  message("Sinking SQUAREM results to file...")
  zz <- file("tmp.txt", open="wt")
  sink(zz, type="message")
  squarem_t <- system.time(
    squarem_res <- squarem(c(fl.to.param(f_init),
                             flash_get_objective(data, f_init)),
                           flash.iter, flash.obj,
                           data = data, m = m, n = n, k = k,
                           control = (list(tol=1, maxiter=niter)))
  )
  sink(type="message")
  squarem_obj <- read.csv("tmp.txt", header = FALSE, sep = '\n')
  file.remove("tmp.txt")

  daarem_t <- system.time(
    daarem_res <- daarem(c(fl.to.param(f_init),
                           flash_get_objective(data, f_init)),
                         flash.iter, flash.obj,
                         data = data, m = m, n = n, k = k,
                         control = (list(tol=1, maxiter=niter)))
  )

  return(list(backfit_obj = backfit_obj, backfit_t = backfit_t,
              squarem_obj = squarem_obj, squarem_t = squarem_t,
              daarem_obj = daarem_res$objfn.track, daarem_t = daarem_t))
}

# Run tests -------------------------------------------------------------

fpath <- "./data/squarem/"

# Normal backfit of fl_g5 takes 304 iterations.
res_5 <- run_test(fl_g5, data, niter = 300)
saveRDS(res_5, paste0(fpath, "res5.rds"))

# Normal backfit of fl_g10 takes 228 iterations.
res_10 <- run_test(fl_g10, data, niter = 225)
saveRDS(res_10, paste0(fpath, "res10.rds"))

# Normal backfit of fl_g20 takes 497 iterations.
res_20 <- run_test(fl_g20, data, niter = 495)
saveRDS(res_20, paste0(fpath, "res20.rds"))


# Use this function if gf and gl parameters aren't converted to a
#   log/logit scale.
#
# flash.iter <- function(p) {
#   init_fl = param.to.fl(p[1:(length(p) - 1)])
#   fl = try(flash_backfit(data, init_fl, ebnm_fn="ebnm_pn",
#                          ebnm_param=list(warmstart=TRUE),
#                          var_type="zero", nullcheck=F, maxiter=1))
#   if (class(fl) == "try-error") {
#     fl = flash_backfit(data, init_fl, ebnm_fn="ebnm_pn",
#                        ebnm_param=list(warmstart=FALSE),
#                        var_type="zero", nullcheck=F, maxiter=1)
#   }
#   return(c(f.to.param(fl), flash_get_objective(data, fl)))
# }

plot_obj <- function(res) {
  data <- c(res$backfit_obj, res$squarem_obj$V1, res$daarem_obj)
  plot(1:length(res$backfit_obj), res$backfit_obj,
       type='l', col='red', ylim=c(min(data), max(data)),
       xlab="Iteration", ylab="Objective")
  lines(1:length(res$squarem_obj$V1), res$squarem_obj$V1,
        col='blue')
  lines(1:length(res$daarem_obj), res$daarem_obj,
        col='green')
}

plot_obj_zoom <- function(res, yrange) {
  data <- c(res$backfit_obj, res$squarem_obj$V1, res$daarem_obj)
  max_obj <- max(data)
  t <- max_obj - yrange
  begin_iter <- min(c(which(res$backfit_obj > t),
                      which(res$squarem_obj > t),
                      which(res$daarem_obj > t)))
  plot(1:length(res$backfit_obj), res$backfit_obj,
       type='l', col='red',
       xlim=c(begin_iter, length(res$backfit_obj)),
       ylim=c(t, max_obj),
       xlab="Iteration", ylab="Objective")
  lines(1:length(res$squarem_obj$V1), res$squarem_obj$V1,
        col='blue')
  lines(1:length(res$daarem_obj), res$daarem_obj,
        col='green')
}

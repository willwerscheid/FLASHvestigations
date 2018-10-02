devtools::load_all("~/GitHub/flashr/")
devtools::load_all("~/GitHub/ebnm/")

# Load data:
gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))
strong <- t(gtex$strong.z)

# Set global flash fit parameters:
ebnm_fn <- "ebnm_ash"
ebnm_param <- list(l = list(), f = list())
tol <- 0.1

# Fit the dataset in the usual way:
control_fl <- flash(strong,
                    var_type = "constant",
                    ebnm_fn = ebnm_fn,
                    ebnm_param = ebnm_param,
                    backfit = TRUE,
                    tol = tol)

# Set a seed and split the dataset into 10 minibatches:
set.seed(666)
p <- ncol(strong)
idx <- sample(1:p)

nbatch <- 10
batchsize <- ceiling(p / nbatch)
batch_idx <- list()
for (i in 1:(nbatch - 1)) {
  batch_idx[[i]] <- idx[((i - 1) * batchsize + 1):(i * batchsize)]
}
batch_idx[[nbatch]] <- idx[((nbatch - 1) * batchsize + 1):p]

# Fit an initial flash object on the first minibatch:
fl_data <- strong[, batch_idx[[1]]]
fl <- flash(fl_data,
            var_type = "constant",
            ebnm_fn = ebnm_fn,
            ebnm_param = ebnm_param,
            backfit = TRUE,
            nullcheck = TRUE,
            tol = tol)
# Save these results as an example of a fit using only subsampled data:
subsample_fl <- fl

# Extract normalized loadings and rescaled priors on factors:
LL <- fl$fit$EL
LL_norms <- sqrt(colSums(LL^2))
LL <- scale(LL, scale = LL_norms, center = FALSE)
gf_sds <- mapply(function(x, y) x$sd * y, fl$fit$gf, as.list(LL_norms))

fixgrid_param = list(l = list(),
                     f = lapply(gf_sds, function(x) {
                       list(mixsd = x, pi_thresh = -1)
                     }))

# Iterate over the minibatches:
for (i in 2:nbatch) {
  message("MINIBATCH ", i)

  # Fit the minibatch to the loadings we have so far:
  fl_data <- flash_set_data(strong[, batch_idx[[i]]])
  old_fl <- flash_add_fixed_loadings(fl_data,
                                     LL,
                                     var_type = "constant",
                                     ebnm_fn = ebnm_fn,
                                     ebnm_param = fixgrid_param,
                                     backfit = TRUE,
                                     nullcheck = FALSE,
                                     tol = tol)

  # Look for new loadings in the current minibatch:
  new_fl <- flash_add_greedy(fl_data,
                             f_init = old_fl,
                             var_type = "constant",
                             ebnm_fn = ebnm_fn,
                             ebnm_param = ebnm_param,
                             nullcheck = FALSE,
                             tol = tol)

  # Unfix all loadings:
  new_fl$fit$fixl[new_fl$fit$fixl] = FALSE
  # Fix the grid for old factors...
  fixgrid_param = list(l = list(),
                       f = lapply(old_fl$fit$gf, function(g) {
                         list(mixsd = g$sd, pi_thresh = -1)
                       }))
  # ...but allow the grid for new factors to change...
  if (ncol(new_fl$fit$EL) > ncol(LL)) {
    for (k in (ncol(LL) + 1):ncol(new_fl$fit$EL)) {
      fixgrid_param$f[[k]] <- list()
    }
  }
  # ...and backfit the flash object to update loadings values:
  new_fl <- flash_backfit(fl_data,
                          new_fl,
                          var_type = "constant",
                          ebnm_fn = ebnm_fn,
                          ebnm_param = fixgrid_param,
                          nullcheck = FALSE,
                          tol = tol)

  # Do a nullcheck here:
  nullcheck_res <- perform_nullcheck(fl_data,
                                     new_fl$fit,
                                     kset = 1:ncol(new_fl$fit$EL),
                                     var_type = "constant",
                                     verbose = TRUE)

  # Normalize loadings, fix, and refit (we can't normalize directly
  #   because the grid is fixed)...
  new_LL <- nullcheck_res$f$EL
  LL_norms <- sqrt(colSums(new_LL^2))
  new_LL <- scale(new_LL, scale = LL_norms, center = FALSE)
  # ...removing any newly added loadings that have been zeroed out...
  to_remove <- nullcheck_res$zeroed_out[nullcheck_res$zeroed_out > ncol(LL)]
  if (length(to_remove) > 0) {
    new_LL <- new_LL[, -to_remove]
    fixgrid_param$f <- fixgrid_param$f[-to_remove]
  }
  # ...and setting any old loadings that have been zeroed out to zero:
  new_LL[is.nan(new_LL)] <- 0
  new_fl <- flash_add_fixed_loadings(fl_data,
                                     new_LL,
                                     var_type = "constant",
                                     ebnm_fn = ebnm_fn,
                                     ebnm_param = fixgrid_param,
                                     backfit = TRUE,
                                     nullcheck = FALSE,
                                     tol = tol)

  # Take weighted average of old and new loadings:
  new_LL[, 1:ncol(LL)] <- ((i - 1) / i) * LL + (1 / i) * new_LL[, 1:ncol(LL)]
  LL <- scale(new_LL, scale = sqrt(colSums(new_LL^2)), center = FALSE)

  # Take weighted average of old and new priors:
  new_gf <- new_fl$fit$gf
  for (k in 1:length(old_fl$fit$gf)) {
    old_pi <- old_fl$fit$gf[[k]]$pi
    new_pi <- new_gf[[k]]$pi
    # Deal with zeroed-out factors separately:
    if (k %in% nullcheck_res$zeroed_out) {
      new_gf[[k]] <- old_fl$fit$gf[[k]] # copy grid, class, etc. over
      new_pi <- c(1, rep(0, length(old_pi) - 1))
    }
    new_gf[[k]]$pi <- ((i - 1) / i) * old_pi + (1 / i) * new_pi
  }
  # For newly added factors, take weighted average of new prior and zero:
  if (length(new_gf) > length(old_fl$fit$gf)) {
    for (k in (length(old_fl$fit$gf) + 1):length(new_gf)) {
      new_pi <- new_gf[[k]]$pi
      if (length(new_pi) > 1) {
        new_pi[1] <- (i - 1) / i + (1 / i) * new_pi[1]
        new_pi[2:length(new_pi)] <- (1 / i) * new_pi[2:length(new_pi)]
      }
      new_gf[[k]]$pi <- new_pi
    }
  }

  fixgrid_param <- list(l = list(),
                        f = lapply(new_gf, function(x) {
                          list(mixsd = x$sd, pi_thresh = -1)
                        }))
}
# Save results:
oneiter_LL <- LL
oneiter_gf <- new_gf

# A second iteration over the same minibatches:
for (i in 1:nbatch) {
  message("MINIBATCH ", i)

  old_gf <- new_gf

  # Fit the minibatch to the loadings without fixing:
  fl_data <- flash_set_data(strong[, batch_idx[[i]]])
  new_fl <- flash_add_fixed_loadings(fl_data,
                                     LL,
                                     fixl = FALSE,
                                     var_type = "constant",
                                     ebnm_fn = ebnm_fn,
                                     ebnm_param = fixgrid_param,
                                     backfit = TRUE,
                                     nullcheck = FALSE,
                                     tol = tol)

  # Do nullcheck here:
  nullcheck_res <- perform_nullcheck(fl_data,
                                     new_fl$fit,
                                     kset = 1:ncol(new_fl$fit$EL),
                                     var_type = "constant",
                                     verbose = TRUE)

  # Normalize loadings, fix, and refit:
  new_LL <- nullcheck_res$f$EL
  LL_norms <- sqrt(colSums(new_LL^2))
  new_LL <- scale(new_LL, scale = LL_norms, center = FALSE)
  new_LL[is.na(new_LL)] <- 0
  new_fl <- flash_add_fixed_loadings(fl_data,
                                     new_LL,
                                     fixl = TRUE,
                                     var_type = "constant",
                                     ebnm_fn = ebnm_fn,
                                     ebnm_param = fixgrid_param,
                                     backfit = TRUE,
                                     nullcheck = FALSE,
                                     tol = tol)

  # Take weighted average of old and new loadings:
  new_LL <- ((nbatch - 1) / nbatch) * LL + (1 / nbatch) * new_LL
  LL <- scale(new_LL, scale = sqrt(colSums(new_LL^2)), center = FALSE)

  # Take weighted average of old and new priors:
  new_gf <- new_fl$fit$gf
  for (k in 1:length(old_gf)) {
    # Deal with zeroed-out factors separately:
    if (k %in% nullcheck_res$zeroed_out) {
      new_gf[[k]] <- old_gf[[k]] # copy grid, etc.
      new_gf[[k]]$pi <- c(1, rep(0, length(old_gf[[k]]$pi) - 1))
    }
    new_gf[[k]]$pi <- ((nbatch - 1) / nbatch) * old_gf[[k]]$pi +
      (1 / nbatch) * new_gf[[k]]$pi
  }
}
# Save results
twoiter_LL <- LL
twoiter_gf <- new_gf

# Do some postprocessing and save the results to file.
#   "Control" is best case (fit all data at once):
control_LL <- control_fl$fit$EL
control_norms <- sqrt(colSums(control_LL^2))
control_LL <- scale(control_LL, scale = control_norms, center = FALSE)
control_gf <- control_fl$fit$gf
for (i in 1:length(control_gf)) {
  control_gf[[i]]$sd <- control_gf[[i]]$sd * control_norms[i]
}
control_ebnm_param <- list(l = list(),
                            f = lapply(control_gf, function(g) {
                              list(g = g, fixg = TRUE)
                            }))
control_fit <- flash_add_fixed_loadings(strong,
                                        LL = control_LL,
                                        var_type = "constant",
                                        ebnm_fn = ebnm_fn,
                                        ebnm_param = control_ebnm_param,
                                        backfit = TRUE,
                                        tol = tol)
saveRDS(control_fit, "./data/minibatch/control_fit.rds")

#   "Subsample" is worst case (only use first minibatch):
subsample_LL <- subsample_fl$fit$EL
subsample_norms <- sqrt(colSums(subsample_LL^2))
subsample_LL <- scale(subsample_LL, scale = subsample_norms, center = FALSE)
subsample_gf <- subsample_fl$fit$gf
for (i in 1:length(subsample_gf)) {
  subsample_gf[[i]]$sd <- subsample_gf[[i]]$sd * subsample_norms[i]
}
subsample_ebnm_param <- list(l = list(),
                             f = lapply(subsample_gf, function(g) {
                               list(g = g, fixg = TRUE)
                             }))
subsample_fit <- flash_add_fixed_loadings(strong,
                                          LL = subsample_LL,
                                          var_type = "constant",
                                          ebnm_fn = ebnm_fn,
                                          ebnm_param = subsample_ebnm_param,
                                          backfit = TRUE,
                                          tol = tol)
saveRDS(subsample_fit, "./data/minibatch/subsample_fit.rds")

oneiter_ebnm_param <- list(l = list(),
                           f = lapply(oneiter_gf, function(g) {
                             list(g = g, fixg = TRUE)
                           }))
oneiter_fit <- flash_add_fixed_loadings(strong,
                                        LL = oneiter_LL,
                                        var_type = "constant",
                                        ebnm_fn = ebnm_fn,
                                        ebnm_param = oneiter_ebnm_param,
                                        backfit = TRUE,
                                        tol = tol)
saveRDS(oneiter_fit, "./data/minibatch/oneiter_fit.rds")

twoiter_ebnm_param <- list(l = list(),
                           f = lapply(twoiter_gf, function(g) {
                             list(g = g, fixg = TRUE)
                           }))
twoiter_fit <- flash_add_fixed_loadings(strong,
                                        LL = twoiter_LL,
                                        var_type = "constant",
                                        ebnm_fn = ebnm_fn,
                                        ebnm_param = twoiter_ebnm_param,
                                        backfit = TRUE,
                                        tol = tol)
saveRDS(twoiter_fit, "./data/minibatch/twoiter_fit.rds")

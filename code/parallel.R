flash_update_fl_parallel = function(data,
                                    f,
                                    kset,
                                    var_type,
                                    ebnm_fn_l,
                                    ebnm_param_l,
                                    ebnm_fn_f,
                                    ebnm_param_f,
                                    parallel_fn) {
  f = flash_update_precision(data, f, var_type)
  f = flash_update_loadings_parallel(data,
                                     f,
                                     kset,
                                     ebnm_fn_l,
                                     ebnm_param_l,
                                     parallel_fn)
  f = flash_update_factors_parallel(data,
                                    f,
                                    kset,
                                    ebnm_fn_f,
                                    ebnm_param_f,
                                    parallel_fn)
}


flash_update_loadings_parallel = function(data,
                                          f,
                                          kset,
                                          ebnm_fn,
                                          ebnm_param,
                                          parallel_fn) {
  R = flash_get_R(data, f)
  subset = !f$fixl

  update_fn = function(k) {
    Rk = R + outer(f$EL[, k], f$EF[, k])
    calc_update_vals(data,
                     f,
                     k,
                     which(subset[, k]),
                     ebnm_fn[[k]],
                     ebnm_param[[k]],
                     loadings = TRUE,
                     Rk)
  }
  res = parallel_fn(as.list(kset), update_fn)

  # Deal with "failed" updates:
  null_idx = which(sapply(res, is.null))
  if (length(null_idx) > 0) {
    res = res[-null_idx]
    kset = kset[-null_idx]
  }
  subset[, -kset] = FALSE

  f$EL[subset] = unlist(lapply(res, function(k) {k$EX}))
  f$EL2[subset] = unlist(lapply(res, function(k) {k$EX2}))
  f$ebnm_fn_l[kset] = ebnm_fn[kset]
  f$ebnm_param_l[kset] = ebnm_param[kset]
  f$gl[kset] = lapply(res, function(k) {k$g})
  f$KL_l[kset] = lapply(res, function(k) {k$KL})

  return(f)
}


flash_update_factors_parallel = function(data,
                                         f,
                                         kset,
                                         ebnm_fn,
                                         ebnm_param,
                                         parallel_fn) {
  R = flash_get_R(data, f)
  subset = !f$fixf

  update_fn = function(k) {
    Rk = R + outer(f$EL[, k], f$EF[, k])
    calc_update_vals(data,
                     f,
                     k,
                     which(subset[, k]),
                     ebnm_fn[[k]],
                     ebnm_param[[k]],
                     loadings = FALSE,
                     Rk)
  }
  res = parallel_fn(as.list(kset), update_fn)

  # Deal with "failed" updates:
  null_idx = which(sapply(res, is.null))
  if (length(null_idx) > 0) {
    res = res[-null_idx]
    kset = kset[-null_idx]
  }

  subset[, -kset] = FALSE
  f$EF[subset] = unlist(lapply(res, function(k) {k$EX}))
  f$EF2[subset] = unlist(lapply(res, function(k) {k$EX2}))
  f$ebnm_fn_f[kset] = ebnm_fn[kset]
  f$ebnm_param_f[kset] = ebnm_param[kset]
  f$gf[kset] = lapply(res, function(k) {k$g})
  f$KL_f[kset] = lapply(res, function(k) {k$KL})

  return(f)
}

library(ashr)

nfactors <- 5

fname <- paste0("../../data/flashier_bench/", data.name, ".rds")
data <- readRDS(fname)

tol <- switch(data.name,
              gtex = 1e-2,
              trachea = 10,
              sp_trachea = 10,
              pulseseq = 10)

ebnm_fn <- switch(prior.family,
                  point.normal = "ebnm_pn",
                  normal.scale.mix = "ebnm_ash")

ebnm_param <- switch(prior.family,
                     point.normal = list(control = list(iterlim = 10)),
                     normal.scale.mix = list(mixcompdist = "normal",
                                             control = list(maxiter.sqp = 10)))

t.init <- system.time({
  fl <- flashr:::flash_greedy_workhorse(data,
                                        Kmax = 1,
                                        var_type = "constant",
                                        ebnm_fn = ebnm_fn,
                                        ebnm_param = ebnm_param,
                                        tol = tol,
                                        maxiter = 1,
                                        nullcheck = FALSE)
})

t.greedy <- system.time({
  fl <- flashr::flash_backfit(data,
                              f_init = fl,
                              var_type = "constant",
                              ebnm_fn = ebnm_fn,
                              ebnm_param = ebnm_param,
                              tol = tol,
                              nullcheck = FALSE)
})

for (k in 2:nfactors) {
  t.init <- t.init + system.time({
    fl <- flashr:::flash_greedy_workhorse(data,
                                          f_init = fl,
                                          Kmax = 1,
                                          var_type = "constant",
                                          ebnm_fn = ebnm_fn,
                                          ebnm_param = ebnm_param,
                                          tol = tol,
                                          maxiter = 1,
                                          nullcheck = FALSE)
  })

  t.greedy <- t.greedy + system.time({
    fl <- flashr::flash_backfit(data,
                                fl,
                                kset = k,
                                var_type = "constant",
                                ebnm_fn = ebnm_fn,
                                ebnm_param = ebnm_param,
                                tol = tol,
                                nullcheck = FALSE)
  })
}

t.backfit <- system.time({
  fl <- flashr::flash_backfit(data,
                              f_init = fl,
                              var_type = "constant",
                              ebnm_fn = ebnm_fn,
                              ebnm_param = ebnm_param,
                              tol = tol * nfactors,
                              nullcheck = FALSE)
})

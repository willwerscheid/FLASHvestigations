nfactors <- 5

fname <- paste0("../../data/flashier_bench/", data.name, ".rds")
data <- readRDS(fname)

tol <- switch(data.name,
              gtex = 1e-2,
              trachea = 10,
              sp_trachea = 10,
              pulseseq = 10)

prior <- switch(prior.family,
                point.normal = flashier::prior.point.normal(),
                normal.scale.mix = flashier::prior.normal.scale.mix())

t.init <- system.time({
  fl <- flashier::flashier(data,
                           greedy.Kmax = 1,
                           var.type = 0,
                           prior.family = prior,
                           greedy.tol = tol,
                           greedy.maxiter = 1,
                           backfit = "none",
                           final.nullchk = FALSE,
                           verbose.lvl = 3,
                           output.lvl = 0)
})

t.greedy <- system.time({
  fl <- flashier::flashier(flash.init = fl,
                           greedy.tol = tol,
                           backfit = "only",
                           warmstart.backfits = FALSE,
                           final.nullchk = FALSE,
                           verbose.lvl = 3,
                           output.lvl = 0)
})

for (k in 2:nfactors) {
  t.init <- t.init + system.time({
    fl <- flashier::flashier(flash.init = fl,
                             greedy.Kmax = 1,
                             greedy.tol = tol,
                             greedy.maxiter = 1,
                             backfit = "none",
                             final.nullchk = FALSE,
                             verbose.lvl = 3,
                             output.lvl = 0)
  })

  t.greedy <- t.greedy + system.time({
    fl <- flashier::flashier(flash.init = fl,
                             greedy.tol = tol,
                             backfit = "only",
                             backfit.kset = k,
                             warmstart.backfits = FALSE,
                             final.nullchk = FALSE,
                             verbose.lvl = 3,
                             output.lvl = 0)
  })
}

t.backfit <- system.time({
  fl <- flashier::flashier(flash.init = fl,
                           greedy.tol = tol,
                           backfit = "only",
                           backfit.maxiter = 500,
                           final.nullchk = FALSE,
                           verbose.lvl = 3,
                           output.lvl = 0)
})

fl.to.par <- function(fl) {
  par <- unlist(flashier:::get.EF(fl))
  par <- c(par, unlist(flashier:::get.EF2(fl)))
  par <- c(par, unlist(flashier:::get.tau(fl)))
  return(par)
}

add.par.to.fl <- function(fl, par, epsilon = 1e-10) {
  dims <- flashier:::get.dims(fl)
  n    <- dims[1]
  p    <- dims[2]
  k    <- flashier:::get.n.factors(fl)

  EFa <- matrix(par[1:(n * k)], nrow = n, ncol = k)
  idx <- n * k
  fl  <- flashier:::set.EF(fl, EFa, n = 1)

  EFb <- matrix(par[(idx + 1):(idx + p * k)], nrow = p, ncol = k)
  idx <- idx + p * k
  fl  <- flashier:::set.EF(fl, EFb, n = 2)

  # Do projections here (ensure that EF2 > EF^2 and tau > 0).
  EF2a <- matrix(par[(idx + 1):(idx + n * k)], nrow = n, ncol = k)
  EF2a <- pmax(EF2a, EFa^2 + epsilon)
  idx  <- idx + n * k
  fl   <- flashier:::set.EF2(fl, EF2a, n = 1)

  EF2b <- matrix(par[(idx + 1):(idx + p * k)], nrow = p, ncol = k)
  EF2b <- pmax(EF2b, EFb^2 + epsilon)
  idx  <- idx + p * k
  fl   <- flashier:::set.EF2(fl, EF2b, n = 2)

  tau <- par[(idx + 1):length(par)]
  tau <- pmax(tau, epsilon)
  fl  <- flashier:::set.tau(fl, tau)

  return(fl)
}

update.factors.serial <- function(flash, kset, cl) {
  for (k in kset) {
    factor <- flashier:::extract.factor(flash, k)
    factor <- flashier:::update.factor(factor, flash, update.tau = FALSE)

    flash <- flashier:::set.EFk(flash, k, flashier:::get.EF(factor))
    flash <- flashier:::set.EF2k(flash, k, flashier:::get.EF2(factor))
    flash <- flashier:::set.KLk(flash, k, flashier:::get.KL(factor))
    flash <- flashier:::set.gk(flash, k, flashier:::get.g(factor))
  }

  flash <- flashier:::init.tau(flash)
  flash <- flashier:::set.obj(flash, flashier:::calc.obj(flash))

  return(flash)
}

extrapolate.flash <- function(fl, maxiter = 100, tol = 0.01,
                              add.step = FALSE, parallel = FALSE) {
  extrapolate <- 2
  beta.init <- 0.5
  beta.increase <- 1.2
  beta.reduce <- 0.5
  betamax <- 2

  if (parallel) {
    update.fn <- flashier:::update.factors.parallel
    cl <- parallel::makeCluster(getOption("cl.cores", 8L),
                                type = getOption("cl.type", "FORK"),
                                useXDR = FALSE)
  } else {
    update.fn <- update.factors.serial
    cl <- NULL
  }

  fl <- flashier:::get.fit(fl)
  K  <- flashier:::get.n.factors(fl)

  all.obj <- rep(NA, maxiter)

  fl <- update.fn(fl, 1:K, cl)
  old.par <- fl.to.par(fl)
  all.obj[1] <- fl$obj

  fl <- update.fn(fl, 1:K, cl)
  new.par <- fl.to.par(fl)
  all.obj[2] <- fl$obj

  iter <- 2
  beta <- beta.init
  f0   <- fl$obj

  while (iter < maxiter) {
    iter <- iter + 1

    propose.par <- new.par + beta * (new.par - old.par)
    fl.propose  <- add.par.to.fl(fl, propose.par)
    fl.propose  <- update.fn(fl.propose, 1:K, cl)
    f <- fl.propose$obj

    if (f < f0) {
      # The objective did not improve.
      fl   <- update.fn(fl, 1:K, cl)
      f    <- fl$obj
      beta <- beta.reduce * beta
    } else {
      # The objective improved.
      if (add.step) {
        fl <- update.fn(fl.propose, 1:K, cl)
      } else {
        fl <- fl.propose
      }
      beta <- min(betamax, beta.increase * beta)
    }

    if (abs(f - f0) < tol)
      break

    all.obj[iter] <- f
    f0 <- f

    old.par <- new.par
    new.par <- fl.to.par(fl)
  }

  all.obj <- all.obj[!is.na(all.obj)]

  if (parallel) {
    parallel::stopCluster(cl)
  }

  return(list(fl = fl, obj = all.obj))
}

daarem.flash <- function(fl, update.fn, ..., control=list()) {
  # added --JW
  n.factors <- fl$n.factors
  fl <- flashier:::get.fit(fl)
  par <- fl.to.par(fl)

  control.default <- list(maxiter=2000, order=10, tol=1.e-08, mon.tol=0.01, cycl.mon.tol=0.0, kappa=25, alpha=1.2)
  namc <- names(control)
  if (!all(namc %in% names(control.default))) {
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  control <- modifyList(control.default, control)

  maxiter <- control$maxiter
  tol <- control$tol
  mon.tol <- control$mon.tol  ## monotonicity tolerance
  cycl.mon.tol <- control$cycl.mon.tol
  a1 <- control$alpha
  kappa <- control$kappa

  num.params <- length(par)
  nlag <- min(control$order, ceiling(num.params/2))

  Fdiff <- Xdiff <- matrix(0.0, nrow=num.params, ncol=nlag)
  obj_funvals <- rep(NA, maxiter + 2)

  # modified --JW
  fl.new <- update.fn(fl, kset = 1:n.factors, cl = NULL)
  xnew <- fl.to.par(fl.new)
  obj_funvals[1] <- fl$obj
  obj_funvals[2] <- fl.new$obj

  likchg <- obj_funvals[2] - obj_funvals[1]
  fn.evals <- 1

  # modified --JW
  fnew <- xnew - par

  k <- 1
  count <- 0
  shrink.count <- 0
  shrink.target <- 1/(1 + a1^kappa)
  lambda.ridge <- 100000
  r.penalty <- 0
  conv <- TRUE
  #num.em <- 0  ## number of EM fallbacks
  ell.star <- obj_funvals[2]

  # added --JW
  fl.propose <- fl.new
  x.propose <- xnew

  while(k < maxiter) {
    count <- count + 1

    # modified --JW
    fl.old <- fl.new
    fl.new <- update.fn(fl.propose, kset = 1:n.factors, cl = NULL)
    fn.evals <- fn.evals + 1
    if (fl.new$obj < fl.old$obj - mon.tol) {
      # Objective decreased too much. Backtrack.
      fl.new <- update.fn(fl.old, kset = 1:n.factors, cl = NULL)
      fn.evals <- fn.evals + 1
      count <- 1
    } else {
      # Objective increased, more or less. Good on you.
      shrink.count <- shrink.count + 1
    }

    # modified --JW
    fold <- fnew
    xold <- xnew
    xnew <- fl.to.par(fl.new)
    fnew <- xnew - xold
    obj_funvals[k+2] <- fl.new$obj
    if(obj_funvals[k+2] - obj_funvals[k+1] < tol & count==nlag) break

    Fdiff[,count] <- fnew - fold
    Xdiff[,count] <- xnew - xold

    np <- count
    Ftmp <- matrix(Fdiff[,1:np], nrow=length(fnew), ncol=np)
    Xtmp <- matrix(Xdiff[,1:np], nrow=length(fnew), ncol=np)  ## is this matrix function needed?

    tmp <- svd(Ftmp)
    dvec <- tmp$d
    uy <- crossprod(tmp$u, fnew)
    uy.sq <- uy*uy

    ### Still need to compute Ftf
    Ftf <- sqrt(sum(as.vector(crossprod(Ftmp, fnew))^2))
    tmp_lam <- daarem:::DampingFind(uy.sq, dvec, a1, kappa, shrink.count, Ftf, lambda.start=lambda.ridge, r.start=r.penalty)
    lambda.ridge <- tmp_lam$lambda
    r.penalty <- tmp_lam$rr

    dd <- (dvec*uy)/(dvec^2 + lambda.ridge)
    gamma_vec <- tmp$v%*%dd

    # modified --JW
    if(class(gamma_vec) == "try-error") {
      fl.propose <- fl.new
    } else {
      xbar <- xnew - drop(Xtmp%*%gamma_vec)
      fbar <- fnew - drop(Ftmp%*%gamma_vec)
      x.propose  <- fbar + xbar
      fl.propose <- add.par.to.fl(fl, x.propose)
    }

    if(count==nlag) {
      count <- 0
      ## restart count
      ## make comparison here l.star vs. obj_funvals[k+2]
      if(obj_funvals[k+2] < ell.star - cycl.mon.tol) {
        ## Decrease delta
        shrink.count <- max(shrink.count - nlag, -2*kappa)
      }
      ell.star <- obj_funvals[k+2]
    }

    shrink.target <-  1/(1 + a1^(kappa - shrink.count))
    k <- k+1
  }
  obj_funvals <- obj_funvals[!is.na(obj_funvals)]

  # modified --JW
  value.obj <- obj_funvals[length(obj_funvals)]

  if(k >= maxiter) {
    conv <- FALSE
    warning("Algorithm did not converge")
  }
  return(list(par=c(xnew), fpevals = k, value.objfn=value.obj, fn.evals=fn.evals, convergence=conv, objfn.track=obj_funvals))
}

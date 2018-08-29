devtools::load_all("/Users/willwerscheid/GitHub/flashr")

sim_mat <- function(n, p, type, missing) {
  if (type == 1) {
    out <- matrix(rnorm(n * p), nrow=n, ncol=p)
  } else if (type == 2) {
    out <- (outer(rep(1, n), rep(1, p))
            + matrix(rnorm(n * p), nrow=n, ncol=p))
  }
  if (missing) {
    out[rbinom(n * p, 1, 0.2) == 1] <- NA
  }
  return(out)
}

do_experiment <- function(nreps, ns, p, type, missing, verbose=TRUE) {
  res <- list()
  if (missing) {
    init_fns <- c("udv_si", "udv_si_svd")
  } else {
    init_fns <- c("udv_svd", "udv_si", "udv_si_svd")
  }

  for (init_fn in init_fns) {
    res[[init_fn]] <- rep(NA, nrow=length(ns))
  }

  for (i in 1:length(ns)) {
    if (verbose) {
      message("n = ", ns[i])
    }

    data <- list()
    for (rep in 1:nreps) {
      data[[rep]] <- sim_mat(ns[i], p, type, missing)
    }

    for (init_fn in init_fns) {
      t <- system.time({
        for (rep in 1:nreps) {
          fl <- flash_add_factors_from_data(data[[rep]],
                                            K=5,
                                            init_fn=init_fn,
                                            backfit=FALSE,
                                            verbose=FALSE)
        }
      })
      res[[init_fn]][i] <- t["elapsed"]
    }

  }

  return(res)
}

set.seed(666)
all_res <- list()

nreps <- 5
ns <- c(10, 25, 50, 100, 250, 500, 1000)
p <- 1000

all_res$null_noNA_p1000 <- do_experiment(nreps, ns, p, type=1,
                                         missing=FALSE)
all_res$r1_noNA_p1000 <- do_experiment(nreps, ns, p, type=2,
                                       missing=FALSE)
all_res$null_missing_p1000 <- do_experiment(nreps, ns, p, type=1,
                                            missing=TRUE)
all_res$r1_missing_p1000 <- do_experiment(nreps, ns, p, type=2,
                                          missing=TRUE)

nreps <- 1
ns <- c(10, 25, 50, 100, 250, 500, 1000)
p <- 10000

all_res$null_noNA_p10000 <- do_experiment(nreps, ns, p, type=1,
                                          missing=FALSE)
all_res$r1_noNA_p10000 <- do_experiment(nreps, ns, p, type=2,
                                        missing=FALSE)
all_res$null_missing_p10000 <- do_experiment(nreps, ns, p, type=1,
                                             missing=TRUE)
all_res$r1_missing_p10000 <- do_experiment(nreps, ns, p, type=2,
                                           missing=TRUE)

saveRDS(all_res, "./data/init_fn3/all_res.rds")

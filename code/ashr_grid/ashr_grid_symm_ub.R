## ashr grid approximations for symmetric unimodal priors
library(tidyverse)
set.seed(666)

log_add <- function(x, y) {
  C <- pmax(x, y)
  x <- x - C
  y <- y - C
  return(log(exp(x) + exp(y)) + C)
}

log_minus <- function(x, y) {
  C <- pmax(x, y)
  x <- x - C
  y <- y - C
  return(log(exp(x) - exp(y)) + C)
}

llik_UN <- function(a, samp) {
  llik1 <- pnorm(a, mean = samp, log.p = TRUE)
  llik2 <- pnorm(-a, mean = samp, log.p = TRUE)
  return(log_minus(llik1, llik2) - log(2 * a))
}

UN_KLdiv <- function(a, omega, aleft, aright, unif_samp, norm_samp) {
  samp <- a * unif_samp + norm_samp
  h_llik <- mean(llik_UN(a, samp))
  llik1 <- log(omega) + llik_UN(aleft, samp)
  llik2 <- log(1 - omega) + llik_UN(aright, samp)
  htilde_llik <- mean(log_add(llik1, llik2))
  return(h_llik - htilde_llik)
}

min_symmKLdiv <- function(a, aleft, aright, unif_samp, norm_samp) {
  optres <- optimize(function(omega) UN_KLdiv(a, omega, aleft, aright, unif_samp, norm_samp),
                     interval = c(0, 1), maximum = FALSE)
  return(optres$objective)
}

ub_symmKLdiv <- function(aleft, aright, unif_samp, norm_samp) {
  optres <- optimize(function(a) min_symmKLdiv(a, aleft, aright, unif_samp, norm_samp),
                     interval = c(aleft, aright), maximum = TRUE)
  return(optres$objective)
}

find_next_gridpt <- function(aleft, targetKL, unif_samp, norm_samp, max_space) {
  uniroot_fn <- function(aright) {ub_symmKLdiv(aleft, aright, unif_samp, norm_samp) - targetKL}
  optres <- uniroot(uniroot_fn, c(aleft + 1e-6, aleft + max_space))
  return(optres$root)
}

build_grid <- function(targetKL, unif_samp, norm_samp, len = 30L) {
  i <- 1
  startpt <- 1e-6
  cat("Gridpoint", i, ":", startpt, "\n")
  i <- 2
  nextpt <- find_next_gridpt(startpt, targetKL, unif_samp, norm_samp, max_space = 10)
  grid <- c(startpt, nextpt)
  cat("Gridpoint", i, ":", max(grid), "\n")
  last_ratio <- 2 * nextpt
  for (i in 3:len) {
    nextpt <- find_next_gridpt(max(grid), targetKL, unif_samp, norm_samp, max(grid) * last_ratio)
    last_ratio <- nextpt / max(grid)
    grid <- c(grid, nextpt)
    cat("Gridpoint", i, ":", max(grid), "\n")
  }
  return(grid)
}

# The sample size again needs to be relatively small.
sampsize <- 100000 # Each grid takes about a minute per grid point
unif_samp <- runif(sampsize, min = -1, max = 1)
norm_samp <- rnorm(sampsize)

df <- NULL
for (KL in 10^(seq(-2, -6, by = -1))) {
  cat("KL:", KL, "\n")
  grid <- build_grid(KL, unif_samp, norm_samp, len = 30L)
  tib <- tibble(KL = KL, grid = grid, idx = 1:length(grid))
  if (is.null(df)) {
    df <- tib
  } else {
    df <- df %>% bind_rows(tib)
  }
}

saveRDS(df, "../../output/ashr_grid/symm_ub_res.rds")

## ashr grid approximations for scale mixtures of normals
set.seed(666)

log_add <- function(x, y) {
  C <- pmax(x, y)
  x <- x - C
  y <- y - C
  return(log(exp(x) + exp(y)) + C)
}

# To "sample" from different a, sample once and then scale depending on a. This method ensures a
#   smooth optimization objective.
smn_KLdiv <- function(a, omega, m, samp) {
  samp <- samp * sqrt(a)
  h_llik <- mean(dnorm(samp, sd = sqrt(a), log = TRUE))
  llik1 <- log(omega) + dnorm(samp, sd = 1, log = TRUE)
  llik2 <- log(1 - omega) + dnorm(samp, sd = sqrt(m), log = TRUE)
  htilde_llik <- mean(log_add(llik1, llik2))
  return(h_llik - htilde_llik)
}

min_smnKLdiv <- function(a, m, samp) {
  optres <- optimize(function(omega) smn_KLdiv(a, omega, m, samp), interval = c(0, 1), maximum = FALSE)
  return(optres$objective)
}

ub_smnKLdiv <- function(m, samp) {
  optres <- optimize(function(a) min_smnKLdiv(a, m, samp), interval = c(1, m), maximum = TRUE)
  return(optres$objective)
}

m <- seq(1.1, 2, by = .05)

cat("UPPER BOUNDS", "\n")

sampsize <- 5000000 # Takes about a minute per value of m
samp <- rnorm(sampsize)

ub <- numeric(length(m))
for (i in 1:length(m)) {
  cat("m:", m[i], "\n")
  ub[i] <- ub_smnKLdiv(m[i], samp)
}

cat("EXPECTED VALUES", "\n")

a <- runif(1000)

sampsize <- 500000 # Takes about 10 minutes per value of m
samp <- rnorm(sampsize)

ev <- numeric(length(m))
for (i in 1:length(m)) {
  cat("m:", m[i], "\n")
  ascaled <- a * (m[i] - 1) + 1
  ev[i] <- mean(sapply(ascaled, function(a) min_smnKLdiv(a, m[i], samp)))
}

res <- list(m = m, ub = ub, ev = ev)

saveRDS(res, "../../output/ashr_grid/smn_res.rds")

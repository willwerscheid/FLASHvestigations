library(flashier)
source("extrapolate_util.R")

options(cl.cores = 8L, cl.type = "FORK")

fl.greedy <- flashier(data, greedy.Kmax = greedy.K, prior.family = prior.family)

# Flashier usually returns the difference in the ELBO after each iteration.
#   Get the raw ELBO instead.
get.elbo <- function(new, old, k) sprintf("%.1f", new$obj)

# 1. Sequential (flashr-style) backfit. --------------------------------------
message("Performing sequential backfit...")
# To track objective, write flashier output to a temporary file.
zz <- file("tmp.txt", open = "wt")
sink(zz)
t <- system.time({
  res <- flashier(init = fl.greedy, fit = "backfit.only", tol = tol,
                  backfit.order = "sequential", backfit.maxiter = maxiter,
                  verbose.fns = get.elbo, verbose.colwidths = 16,
                  verbose.colnames = "elbo", verbose.lvl = -1)
})
sink()
close(zz)
res <- read.table("tmp.txt", header = TRUE)
zz <- file.remove("tmp.txt")
all.res <- data.frame(fit = "SerialFlashr",
                      elbo = res$elbo,
                      elapsed.time = t[3] * 1:nrow(res) / nrow(res))


# 2. Dropout (flashier-style) backfit. ----------------------------------------
message("Performing dropout backfit...")
zz <- file("tmp.txt", open = "wt")
sink(zz)
t <- system.time({
  res <- flashier(init = fl.greedy, fit = "backfit.only", tol = tol,
                  backfit.order = "dropout", backfit.maxiter = maxiter,
                  verbose.fns = get.elbo, verbose.colwidths = 16,
                  verbose.colnames = "elbo", verbose.lvl = -1)
})
sink()
close(zz)
res <- read.table("tmp.txt", header = TRUE)
zz <- file.remove("tmp.txt")
next.res <- data.frame(fit = "SerialFlashier",
                       elbo = res$elbo,
                       elapsed.time = t[3] * 1:nrow(res) / nrow(res))
all.res <- rbind(all.res, next.res)

# 3. Parallel backfit. --------------------------------------------------------
message("Performing parallel backfit...")
zz <- file("tmp.txt", open = "wt")
sink(zz)
t <- system.time({
  res <- flashier(init = fl.greedy, fit = "backfit.only", tol = tol,
                  backfit.order = "parallel", backfit.maxiter = maxiter,
                  verbose.fns = get.elbo, verbose.colwidths = 16,
                  verbose.colnames = "elbo", verbose.lvl = -1)
})
sink()
close(zz)
res <- read.table("tmp.txt", header = TRUE)
zz <- file.remove("tmp.txt")
next.res <- data.frame(fit = "Parallel",
                       elbo = res$elbo,
                       elapsed.time = t[3] * 1:nrow(res) / nrow(res))
all.res <- rbind(all.res, next.res)

# 4. Serial extrapolation scheme. ---------------------------------------------
#   Update all factors serially; extrapolate EF, EF2, and tau; repeat.
message("Using serial extrapolation scheme...")
t <- system.time({
  res <- extrapolate.flash(fl.greedy, maxiter = maxiter, tol = tol)
})
next.res <- data.frame(fit = "Extrapolate",
                       elbo = res$obj,
                       elapsed.time = t[3] * 1:length(res$obj) / length(res$obj))
all.res <- rbind(all.res, next.res)

# 5. Parallel extrapolation scheme. -------------------------------------------
#   Same as #4, but updates are done in parallel. Monotonicity is not
#   guaranteed.
message("Using parallel extrapolation scheme...")
t <- system.time({
  res <- extrapolate.flash(fl.greedy, maxiter = maxiter, tol = tol,
                           parallel = TRUE)
})
next.res <- data.frame(fit = "Parapolate",
                       elbo = res$obj,
                       elapsed.time = t[3] * 1:length(res$obj) / length(res$obj))
all.res <- rbind(all.res, next.res)

# 6. Serial stutter-step extrapolation scheme. --------------------------------
#   Update twice before extrapolating: update, update, extrapolate, repeat.
message("Serial stutter-step extrapolation scheme...")
t <- system.time({
  res <- extrapolate.flash(fl.greedy, maxiter = maxiter / 2, tol = tol,
                           add.step = TRUE)
})
next.res <- data.frame(fit = "Extrastutter",
                       elbo = res$obj,
                       elapsed.time = t[3] * 1:length(res$obj) / length(res$obj))
all.res <- rbind(all.res, next.res)

# 7. Parallel stutter-step extrapolation scheme. ------------------------------
t <- system.time({
  res <- extrapolate.flash(fl.greedy, maxiter = maxiter / 2, tol = tol,
                           add.step = TRUE, parallel = TRUE)
})
next.res <- data.frame(fit = "Parastutter",
                       elbo = res$obj,
                       elapsed.time = t[3] * 1:length(res$obj) / length(res$obj))
all.res <- rbind(all.res, next.res)

saveRDS(all.res, out)

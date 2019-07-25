data <- readRDS("../../data/flashier_bench/sp_trachea.rds")
out  <- "../../output/extrapolate/trachea_pnres.rds"

maxiter <- 300
tol <- 2.0
greedy.K <- 40
prior.family <- flashier::prior.point.normal()

source("extrapolate.R")

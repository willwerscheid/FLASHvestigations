data <- readRDS("../../data/flashier_bench/sp_trachea.rds")
out  <- "../../output/extrapolate/trachea_ashres.rds"

maxiter <- 300
tol <- 2.0
greedy.K <- 100
prior.family <- flashier::prior.normal.scale.mix()

source("extrapolate.R")

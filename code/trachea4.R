# Download and pre-process pathways -----------------------------------------

# I want the pathways matrix to be a matrix of ones and zeros, with rows
#   corresponding to genes and columns corresponding to pathways.
# I download the following file from reactome.org:
#   wget https://reactome.org/download/current/Ensembl2Reactome.txt
setwd("/project2/mstephens/willwerscheid/")
pathways <- data.table::fread("Ensembl2Reactome.txt",
                              header = FALSE,
                              stringsAsFactors = TRUE)
pathways <- subset(pathways, V6 == "Mus musculus")
pw.mat <- reshape2::acast(pathways, V1 ~ V4, length)
pw.mat <- Matrix::Matrix(pw.mat)
saveRDS(pw.mat, "mus_pathways.rds")


# Download and pre-process data ---------------------------------------------

# File GSE103354_Trachea_droplet_UMIcounts.txt.gz was downloaded at
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103354
trachea <- read.table("droplet.txt")
trachea <- Matrix::Matrix(as.matrix(trachea))
# scale.constants <- 1e6 / Matrix::colSums(trachea)
# saveRDS(scale.constants, "trachea_scale.rds")

# Convert the MGI gene symbols to Ensembl IDs.
gene.symbols <- rownames(trachea)
ensembl <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
bm <- biomaRt::getBM(filters = "mgi_symbol",
                     attributes = c("mgi_symbol", "ensembl_gene_id"),
                     values = gene.symbols,
                     mart = ensembl)
saveRDS(bm, "pw_bm.rds")

trachea <- trachea[which(gene.symbols %in% bm$mgi_symbol), ]
gene.symbols <- rownames(trachea)
ensembl.ids <- bm$ensembl_gene_id[match(gene.symbols, bm$mgi_symbol)]
rownames(trachea) <- ensembl.ids

# I drop genes that have nonzero counts in five or fewer cells.
gene.cts <- Matrix::rowSums(trachea > 0)
trachea <- trachea[which(gene.cts > 5), ]
saveRDS(trachea, "trachea.rds")


# Fit the pathways to the data ----------------------------------------------

setwd("/project2/mstephens/willwerscheid/")
trachea <- readRDS("trachea.rds")
pw.mat <- readRDS("mus_pathways.rds")

library(ashr)
library(flashier)
library(ebnm)
library(Matrix)
library(susieR)

common.genes <- rownames(trachea)[rownames(trachea) %in% rownames(pw.mat)]
trachea <- trachea[common.genes, ]
pw.mat <- pw.mat[common.genes, ]

# Remove pathways with five or fewer genes.
pw.mat <- pw.mat[, colSums(pw.mat) > 5]
# Some genes no longer belong to any pathway. Remove them.
to.remove <- which(rowSums(pw.mat) == 0)
trachea <- trachea[-to.remove, ]
pw.mat <- pw.mat[-to.remove, ]
saveRDS(trachea, "pw_trachea.rds")
saveRDS(pw.mat, "pw_mat.rds")

# Log-transform the data and set it for flash.
data <- log1p(trachea)
S <- sqrt(data) / (data + 1)
S[S == 0] <- 0.5
fl.data <- set.flash.data(data, S, var.type = NULL)
rm(S)

# Subset ionocytes and set for flash.
ion.idx <- (sapply(strsplit(colnames(trachea), "_"), `[`, 3) == "Ionocyte")
ionocytes <- trachea[, ion.idx]
saveRDS(ionocytes, "ionocytes.rds")
# all-zero rows need to be removed (for flash only, not SuSiE):
nz.idx <- (rowSums(ionocytes) > 0)
ionocytes <- log1p(ionocytes[nz.idx, ])
S <- sqrt(ionocytes) / (ionocytes + 1)
S[S == 0] <- 0.5
ion.data <- set.flash.data(ionocytes, S, var.type = NULL)
rm(S)

# Time to flash. First, fit the mean factors. For genes, a "nonnegative" prior
#   (which is unimodal at zero) works best. I don't expect cell sizes to be
#   unimodal at zero, so I use a "nonzero.mode" prior for cells.
fl <- flashier(fl.data,
               var.type = NULL,
               prior.type = c("nonnegative", "nonzero.mode"),
               fixed.factors = c(ones.factor(2), ones.factor(1)),
               backfit = "only",
               verbose.lvl = 3)
saveRDS(fl, "/scratch/midway2/willwerscheid/pw_tmp.rds")

# For the cell type subsets, I only want a mean factor for cells.
ion.fl <- flashier(ion.data,
                   var.type = NULL,
                   prior.type = c("nonnegative", "nonzero.mode"),
                   fixed.factors = ones.factor(2),
                   greedy.Kmax = 0,
                   verbose.lvl = 3)
saveRDS(ion.fl, "/scratch/midway2/willwerscheid/ion_tmp.rds")

# I add a pathway as follows: First, I add a factor greedily with a nonnegative
#   prior on genes (since one expects genes in a pathway to be regulated in the
#   same direction).
get.next.greedy <- function(flash.init) {
  next.greedy <- flashier(flash.init = flash.init,
                          prior.type = c("nonnegative", "normal.mix"),
                          ash.param = list(method = "fdr"),
                          greedy.Kmax = 1,
                          verbose.lvl = 3)
  return(next.greedy$loadings$normalized.loadings[[1]][, next.greedy$n.factors])
}

# I then use SuSiE to find the pathways that best match the new factor.
get.candidate.pathways <- function(next.greedy, pw.mat, idx = NULL, level = 0.95) {
  if (!is.null(idx)) {
    Y <- rep(0, nrow(pw.mat))
    Y[idx] <- next.greedy
  }
  suz <- susie(pw.mat, Y, L = 1)
  suz.order <- order(suz$alpha, decreasing = TRUE)
  levels <- cumsum(suz$alpha[suz.order])
  n.candidates <- sum(levels < level) + 1
  return(suz.order[1:n.candidates])
}

# For each candidate pathway, I run flash to add the pathway as a fixed sparse
#   factor (i.e., all genes not in the pathway are fixed at zero). I keep the
#   pathway that gives the best flash objective.
fit.pathway <- function(fl, pw.mat, next.pathway, idx = NULL) {
  if (is.null(idx))
    idx <- 1:nrow(pw.mat)
  pathway.idx <- which(pw.mat[idx, next.pathway] == 1)
  fl <- flashier(flash.init = fl,
                 prior.type = c("nonnegative", "normal.mix"),
                 fixed.factors = sparse.factors(1, pathway.idx),
                 greedy.Kmax = 0,
                 ash.param = list(method = "fdr"),
                 verbose.lvl = 3)
  return(fl)
}

add.pathway <- function(fl, pw.mat, nz.idx = NULL, level = 0.95, save.fl = TRUE) {
  cat("FITTING NEXT GREEDY FACTOR \n")
  Y <- get.next.greedy(fl)

  cat("FINDING CANDIDATE PATHWAYS \n")
  candidates <- get.candidate.pathways(Y, pw.mat, nz.idx, level)
  n.candidates <- length(candidates)

  for (i in 1:n.candidates)
    cat("CANDIDATE PATHWAY:", colnames(pw.mat)[candidates[i]], "\n")

  best.obj <- -Inf
  for (i in 1:n.candidates) {
    cat("FITTING PATHWAY:", colnames(pw.mat)[candidates[i]], "\n")
    next.fl <- fit.pathway(fl, pw.mat, candidates[i], nz.idx)
    if (next.fl$objective > best.obj) {
      best.obj <- next.fl$objective
      best.fl <- next.fl
      best.pathway <- candidates[i]
      cat("NEW BEST PATHWAY:", colnames(pw.mat)[best.pathway], "\n")
    }
  }

  # Check that the pathway improves the objective before adding.
  if (best.fl$objective <= fl$objective)
    return(fl)

  cat("ADDING PATHWAY:", colnames(pw.mat)[best.pathway], "\n")
  if (save.fl)
    saveRDS(best.fl, "/scratch/midway2/willwerscheid/fl_tmp.rds")
  return(best.fl)
}

add.n.pathways <- function(fl, pw.mat, n, nz.idx = NULL, level = 0.95, save.fl = TRUE) {
  for (i in 1:n) {
    new.fl <- add.pathway(fl, pw.mat, nz.idx, level, save.fl)
    if (identical(fl, new.fl))
      break
    fl <- new.fl
  }
  return(fl)
}

backfit.fl <- function(fl) {
  return(flashier(flash.init = fl,
                  backfit = "only", backfit.reltol = 5,
                  verbose.lvl = 3, output.lvl = 4))
}

# Full data.
fl <- add.n.pathways(fl, pw.mat, 10)
fl <- backfit.fl(fl)
saveRDS(fl, "pw_fl10.rds")

fl <- add.n.pathways(fl, pw.mat, 10)
fl <- backfit.fl(fl)
saveRDS(fl, "pw_fl20.rds")

fl <- add.n.pathways(fl, pw.mat, 10)
fl <- backfit.fl(fl)
saveRDS(fl, "pw_fl23.rds") # didn't find a pathway that increased the objective

fl <- add.n.pathways(fl, pw.mat, 10)
fl <- backfit.fl(fl)
saveRDS(fl, "pw_fl33.rds")

fl <- add.n.pathways(fl, pw.mat, 10)
fl <- backfit.fl(fl)
saveRDS(fl, "pw_fl41.rds") # too many candidate factors

fl <- add.n.pathways(fl, pw.mat, 9)
fl <- backfit.fl(fl)
saveRDS(fl, "pw_fl50.rds")

# Ionocytes.
ion.fl <- add.n.pathways(ion.fl, pw.mat, 10, nz.idx, save.fl = FALSE)
ion.fl <- backfit.fl(ion.fl)
saveRDS(ion.fl, "ion_fl10.rds")

ion.fl <- add.n.pathways(ion.fl, pw.mat, 10, nz.idx, save.fl = FALSE)
ion.fl <- backfit.fl(ion.fl)
saveRDS(ion.fl, "ion_fl12.rds") # didn't find a pathway that increased obj

# This subsequent fit is unable to add any other pathways.
ion.fl <- add.n.pathways(ion.fl, pw.mat, 10, nz.idx, save.fl = FALSE)

# Recover pathways from flash object.
pathway.idx <- lapply(fl$fit$fix.idx[1:fl$n.factors], function(f) {
  setdiff(1:nrow(fl$fit$EF[[1]]), f)
})
pathways <- apply(pw.mat, 2, function(x) which(x == 1))
added.idx <- match(pathway.idx, pathways)
pathways.added <- names(pathways)[added.idx]

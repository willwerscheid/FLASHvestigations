# Download GTEx data from GTEx Portal and load into R using Peter's code ------

samples.file <- "https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt"
read.counts.file <- "~/Downloads/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct.gz"

source("https://raw.githubusercontent.com/stephenslab/topics/master/code/gtex.R?token=AWTMG_VIxH9P52pyv6O3a3tRQEFn0F9Zks5cN431wA%3D%3D")
gtex <- read.gtex.data(samples.file, read.counts.file)


# Pre-process the data as described in Wang, Fischer, and Song (2017) ---------

# "Normalization was performed using the size factors produced by the
#   estimateSizeFactors function of DESeq2."

#     Calculate the (logged) geometric means of the counts for each gene.
gene.geom.means <- apply(gtex$counts, 2, function(x) {
  if (all(x == 0)) {-Inf} else {sum(log(x[x > 0])) / length(x)}
})
#     Take the median of the normalized counts for each sample:
size.factors <- apply(gtex$counts, 1, function(x) {
  exp(median((log(x) - gene.geom.means)[x > 0]))
})
#     Normalize the samples:
gtex$counts <- apply(gtex$counts, 2, `/`, size.factors)

# "We required at least 15 samples to include a given tissue and an average of
#   at least 500 normalized reads in one or more tissues to retain a gene."

tissues.to.include <- names(which(table(gtex$samples$specific) > 14))
samples.to.retain <- samples$specific %in% tissues.to.include
gtex$counts <- gtex$counts[samples.to.retain, ]
gtex$samples <- gtex$samples[samples.to.retain, ]

gene.max.avg.reads <- apply(gtex$counts, 2, function(x) {
  max(aggregate(x, by = list(gtex$samples$specific), FUN = mean)$x)
})
genes.to.retain <- gene.max.avg.reads > 500
gtex$counts <- gtex$counts[, genes.to.retain]

# Log-transform data:

gtex$counts <- log1p(gtex$counts)


# Convert the data from a matrix to a tissue x individual x gene array --------

individual <- as.factor(sapply(strsplit(rownames(gtex$counts), "-"), function(x) {
  paste(x[1:2], collapse = "-")
}))
gtex.df <- data.frame(individual = individual,
                      tissue = droplevels(gtex$samples$specific))
gtex.df <- cbind(gtex.df, gtex$counts)
gtex <- reshape2::melt(gtex.df)
colnames(gtex)[3] <- "gene"
rm(gtex.df)

gtex <- reshape2::acast(gtex, tissue ~ individual ~ gene)
# object size: 4.6 Gb (a 49 x 714 x 17792 array)

saveRDS(gtex, "~/Downloads/gtex_v7_array.rds") # temporary


# Create smaller array using only brain tissues -------------------------------

brain.tissues <- which(substr(dimnames(gtex)[[1]], 1, 5) == "brain")
brain <- gtex[brain.tissues, , ]
# Remove individuals with no brains:
brain <- brain[, apply(brain, 2, function(x) sum(!is.na(x))) > 0, ]
# object size: 450 Mb (a 13 x 254 x 17792 array)

saveRDS(brain, "~/Downloads/gtex_v7_brain.rds") # temporary


# Fit flash object ------------------------------------------------------------

devtools::load_all("~/Github/ashr")
devtools::load_all("~/Github/flashier")
devtools::load_all("~/Github/ebnm")

brain <- set.flash.data(brain) # object size: 675 Mb
brain.flash <- flashier(brain,
                        var.type = 3,
                        prior.type = c("nonnegative", "nonzero.mode", "normal.mix"),
                        conv.crit.fn = function(new, old, k) {
                          flashier:::calc.max.abs.chg.EF(new, old, k, n = 1)
                        },
                        greedy.Kmax = 10,
                        greedy.tol = 5e-4,
                        backfit.after = 2,
                        backfit.every = 1,
                        inner.backfit.maxiter = 1,
                        ash.param = list(optmethod = "mixSQP"),
                        nonmissing.thresh = c(0, 0.05, 0),
                        verbose = "O L1 E2")

saveRDS(brain.flash, "~/Downloads/brain_flash.rds") # temporary

brain.flash <- flashier(brain,
                        flash.init = brain.flash,
                        backfit = "only",
                        backfit.order = "dropout",
                        conv.crit.fn = function(new, old, k) {
                          flashier:::calc.max.abs.chg.EF(new, old, k, n = 1)
                        },
                        backfit.tol = 5e-4,
                        verbose = "O L1 E2")

saveRDS(brain.flash, "~/Downloads/brain_flash_bf.rds") # temporary

# Remove large fields from flash object.
brain.flash$fit$Y <- NULL
brain.flash$fit$Z <- NULL
brain.flash$sampler <- NULL

saveRDS(brain.flash, "./data/brain/brain05.rds")

# Other "brain" objects were created using different settings for
#   nonmissing.thresh.


# Create data frame containing demographics and technical factors -------------

phenotypes <- readRDS("~/Downloads/GTEx_v7_Subjects.rds")
class(phenotypes) <- "data.frame"
rownames(phenotypes) <- phenotypes$SUBJID
phenotypes <- phenotypes[, -1]
# Get subset corresponding to individuals with brains.
phenotypes <- phenotypes[rownames(brain.flash$loadings$normalized.loadings[[2]]), ]
# Remove COHORT because they're all postmortem.
phenotypes <- phenotypes[, -1]
# Remove ETHNCTY because there are no hispanics.
phenotypes <- phenotypes[, -4]

samples <- readr::read_delim("https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt",
                             delim = "\t")
samples <- samples[, c("SAMPID", "SMGEBTCHD")]
# Extract individuals from sample IDs.
samples$SAMPID <- sapply(lapply(strsplit(samples$SAMPID, "-"), `[`, 1:2),
                         paste, collapse = "-")
# Convert sequencing date to number of days after first sequencing date.
samples$SMGEBTCHD <- as.numeric(as.POSIXct(samples$SMGEBTCHD, format = "%m/%d/%Y"))
samples$SMGEBTCHD <- samples$SMGEBTCHD - min(samples$SMGEBTCHD, na.rm = TRUE)
samples$SMGEBTCHD <- samples$SMGEBTCHD / (60 * 60 * 24)

SEQDATE <- aggregate(SMGEBTCHD ~ SAMPID, data = samples,
                       FUN = mean, na.rm = TRUE)
rownames(SEQDATE) <- SEQDATE$SAMPID
SEQDATE <- SEQDATE[, 2, drop = FALSE]
SEQDATE <- SEQDATE[rownames(brain.flash$loadings$normalized.loadings[[2]]), ]

all.covar <- cbind(phenotypes, SEQDATE)

# Do not upload to GitHub (data includes protected attributes!)
saveRDS(all.covar, "~/Downloads/GTEx_v7_brain_covar.rds")

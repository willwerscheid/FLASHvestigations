library(Matrix)
library(mgcv)
library(ggplot2)
library(mixsqp)

# Functions from https://github.com/adamh-broad/single_cell_airway:

extract.field <- function(string, field = 1, delim = "_", fixed = TRUE) {
  return(strsplit(string,delim, fixed=fixed)[[1]][field])
}

info <- function(text, ...) {
  cat(sprintf(paste(Sys.time(),"INFO:", text,"\n")))
}

tpm <- function(counts, mult = 10000) {
  info("Running TPM normalisation")
  total.counts = Matrix::colSums(counts)
  scaled.counts = t(t(counts) / total.counts)
  scaled.counts * mult
}

get.variable.genes.umis <- function(umi.cts, residual.threshold = -0.25,
                                    use.spline = FALSE, batch = NULL,
                                    ret.plot = FALSE) {
  if(!is.null(batch)){
    v = as.vector(table(batch))
    total_transcripts = data.frame(as.matrix(t(Matrix.utils::aggregate.Matrix(t(umi.cts), groupings = batch, fun="sum"))))
    detection_frac = Matrix.utils::aggregate.Matrix(t(umi.cts > 0), groupings = batch, fun="sum")
    detection_frac = data.frame(as.matrix(t(detection_frac / v)))
    test_genes = rownames(detection_frac)[Matrix::rowSums(detection_frac > 0) == length(unique(batch))]
    detection_frac = detection_frac[test_genes, ]
    total_transcripts = total_transcripts[test_genes, ]
    detection_frac$gene = rownames(detection_frac)
    total_transcripts$gene = rownames(total_transcripts)
    detection_frac = melt(detection_frac, id.vars="gene")
    colnames(detection_frac) = c("gene", "batch", "alpha")
    total_transcripts = melt(total_transcripts, id.vars="gene")
    colnames(total_transcripts) = c("gene", "batch", "UMIs")
    z = cbind(total_transcripts, detection_frac)[, c("gene", "batch", "alpha", "UMIs")]
    info("Fitting logistic GLM (controlling for batch covariate)")
    model.logit = glm(data = z, formula = alpha ~ log10(UMIs) + batch, family = binomial)
    #model.logit = robust::glmRob(data = z, formula = alpha ~ log10(UMIs), family = binomial)
    info("Fitting spline quantile regression (controlling for batch covariate)")
    model.gam = quantreg::rq(data = z, formula = alpha ~ splines::ns(log10(UMIs), df=15) + batch, tau=0.8)
    #model.gam = mgcv::gam(data = z, formula = alpha ~ s(log10(UMIs)), method="REML")
  }else{
    info("Computing gene dection rates (alphas)..")
    z = data.frame(UMIs = Matrix::rowSums(umi.cts), alpha= Matrix::rowSums(umi.cts>0) / ncol(umi.cts))
    z = subset(z, UMIs > 0 | alpha > 0)
    info("Fitting GLMs..")
    model.logit = glm(data = z, formula = alpha ~ log10(UMIs), family = binomial)
    #model.logit = robust::glmRob(data = z, formula = alpha ~ log10(UMIs), family = binomial)
    model.gam = quantreg::rq(data = z, formula = alpha ~ splines::ns(log10(UMIs), df=15), tau=0.8)
    #model.gam = mgcv::gam(data = z, formula = alpha ~ s(log10(UMIs)), method="REML")
  }

  if(use.spline){
    info("use.spline is ON. Using GAM fit (blue), logit in red")
    z$predicted = predict(object = model.gam, z, type="response")
    z$predicted.alternate = predict(object = model.logit, z, type="response")
    z$residual = model.gam$residuals
  }else{
    info("use.spline is OFF. Using logit fit (blue), GAM in red")
    z$predicted = predict(model.logit, type="response") #predict(object = model.logit, z, type="response")
    z$predicted.alternate = predict(object = model.gam, z, type="response")
    z$residual = model.logit$residuals
  }
  if(is.null(batch)) {z$gene = rownames(z)}
  outliers = subset(z, residual < residual.threshold)
  g = ggplot(z, aes(x=log10(UMIs), y=alpha, label=gene)) + geom_point(color="grey50", size=0.5, stroke=0) +
    ylim(c(0,1)) + geom_line(aes(y=predicted), size=0.5, color="blue", linetype="dotted") +
    geom_line(aes(y=predicted.alternate), size=0.5, color="red", linetype="dotted") +
    geom_text(data=outliers, color="black", size=1.5, vjust=2)
  if(!is.null(batch)){g = g + facet_wrap(~batch)}
  if(!ret.plot){print(g)}
  rv = unique(unlist(lapply(rownames(outliers), extract.field, 1, delim="_")))
  if(ret.plot){return(list("var.genes"=rv, "plot"=g, "fit.data"=z, "logit"=model.logit))}
  rv
}

# Load and preprocess data as in Montoro et al.
counts <- read.table("~/Downloads/GSE103354_Trachea_droplet_UMIcounts.txt")
counts <- Matrix::Matrix(as.matrix(counts))
var.genes <- get.variable.genes.umis(counts)
counts <- tpm(counts)
counts <- log1p(counts)
counts <- counts[var.genes, ]

# Add 30 factors with rough backfits after every 5 factors.
fl <- flashier(counts, greedy.Kmax = 30, var.type = 1,
               prior.type = c("nonnegative", "normal.mix"),
               ash.param = list(optmethod = "mixSQP"),
               backfit.every = 5, final.backfit = TRUE,
               backfit.order = "montaigne", warmstart.backfits = FALSE,
               verbose.lvl = 3)

saveRDS(fl, "~/Downloads/ppdrop30.rds")

# Refine by backfitting.
fl <- flashier(counts, flash.init = fl, backfit = "only",
               backfit.order = "dropout", backfit.maxiter = 200,
               verbose.lvl = 3)

saveRDS(fl, "~/Downloads/ppdrop30_backfit.rds")

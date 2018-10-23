library(mixsqp)
devtools::load_all("~/GitHub/flashrtools/")

t_fl30 <- system.time({
  fl30 <- flashier(trachea,
                   var_type = "by_row",
                   method = "nnloadings",
                   greedy_Kmax = 30,
                   backfit_maxiter = 100)
})

# Fit a hierarchical FLASH object to ensure there are no serious problems:
hier_fl <- flashier(fl30$ldf$f,
                    greedy_Kmax = 30,
                    var_type = "by_row",
                    method = "fastest")
plot(hier_fl, plot_factors = TRUE, plot_scree = FALSE)

# Extract factors and append cell type (as established by Montoro et al.):
factors_df <- data.frame(fl30$ldf$f)

cell_types <- sapply(strsplit(colnames(trachea), "_"), function(x) x[3])
cell_types <- as.factor(cell_types)
lvls <- levels(cell_types)
lvls[lvls == "Neuroendocrine"] <- "NEC"
levels(cell_types) <- lvls

factors_df$cell_type <- cell_types
saveRDS(factors_df, "~/GitHub/FLASHvestigations/data/trachea/factors_df.rds")

top_genes <- get_top_loading_elements(fl30)
top_genes <- add_desc_to_top_genes(top_genes,
                                   dataset = "mmusculus_gene_ensembl")
saveRDS(top_genes, "~/GitHub/FLASHvestigations/data/trachea/top_genes.rds")

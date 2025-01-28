# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                     ---- DIFFERENTIAL CO-EXPRESSION ANALYSIS ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
library(Seurat)
library(dplyr)
library(tidyverse)
library(viridisLite)
library(viridis)
library(ggplot2)
library(openxlsx)
library(future)
library(slingshot)
library(MAST)
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(Hmisc)
library(igraph)
library(corrplot)
library(ggplot2)
library(gridExtra)

library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(reshape2)
library(igraph)



rm(list = ls())
set.seed(18)
writeExcel = TRUE
options(future.globals.maxSize=10000*1024^2)



setwd(getwd())
output_folder <- "./Processed"


S <- LoadSeuratRds(file = './Processed/SeuratFinal.rds')


n_cores <- 16
# Enable parallel computing
enableWGCNAThreads(nThreads = n_cores)

# 1. Prepare pseudobulk expression data
# ====================================================
# Use imputed data for co-expression analysis
DefaultAssay(S) <- "RNA_imputed"

# Aggregate expression by Condition and Day
Idents(S) <- "ConditionDay"
pseudobulk <- AggregateExpression(S, 
                                  assays = "RNA_imputed", 
                                  return.seurat = TRUE, 
                                  group.by = "ConditionDay",
                                  slot = "data")  # Use "data" instead of "counts"

expr_matrix <- GetAssayData(S, assay = "RNA", slot = "data")
hvg <- VariableFeatures(S)
expr_matrix <- expr_matrix[hvg, ]

# Extract expression matrices
wt_expr <- GetAssayData(pseudobulk[, grepl("WT", colnames(pseudobulk))], 
                        assay = "RNA_imputed", 
                        slot = "counts")
pd_expr <- GetAssayData(pseudobulk[, grepl("ND", colnames(pseudobulk))], 
                        assay = "RNA_imputed", 
                        slot = "counts")

# Transpose for WGCNA (genes as columns, samples as rows)
wt_matrix <- t(as.matrix(wt_expr))
pd_matrix <- t(as.matrix(pd_expr))

# 2. Network construction with WGCNA
# ====================================================

construct_network <- function(datExpr, condition_name) {
  # Determine soft-thresholding power
  powers <- c(1:20)
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  power <- sft$powerEstimate
  
  if(is.na(power)) {
    power <- 6 # Fallback if no power is found
    warning(paste("No soft threshold found for", condition_name, "using power = 6"))
  }
  
  # Construct network
  net <- blockwiseModules(datExpr,
                          power = power,
                          TOMType = "unsigned",
                          minModuleSize = 30,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          numericLabels = TRUE,
                          pamRespectsDendro = FALSE,
                          saveTOMs = TRUE,
                          saveTOMFileBase = paste0(condition_name, "_TOM"),
                          verbose = 3)
  
  return(net)
}

# Build networks for both conditions
wt_net <- construct_network(wt_matrix, "WT")
pd_net <- construct_network(pd_matrix, "PD")

# 3. Module preservation analysis
# ====================================================

multiExpr <- list(WT = list(data = wt_matrix), 
                  PD = list(data = pd_matrix))

mp <- modulePreservation(multiExpr,
                         wt_net$colors,
                         referenceNetworks = 1,
                         nPermutations = 100,
                         randomSeed = 1,
                         verbose = 3)

# Analyze preservation statistics
preservation_stats <- mp$preservation$Z$ref.WT$inColumnsAlsoPresentIn.PD
print(head(preservation_stats[order(-preservation_stats$Zsummary.pres),], 20))

# 4. Differential co-expression analysis
# ====================================================

# Calculate correlation matrices
wt_cor <- cor(wt_matrix, method = "spearman")
pd_cor <- cor(pd_matrix, method = "spearman")

# Identify significantly different correlations
diff_cor <- wt_cor - pd_cor
cor_pvals <- corPvalueFisher(diff_cor, nSamples = nrow(wt_matrix))

# Threshold for significant differences
sig_threshold <- 0.01
sig_diff <- which(abs(diff_cor) > 0.6 & cor_pvals < sig_threshold, arr.ind = TRUE)

# Create differential co-expression network
diff_network <- graph_from_adjacency_matrix(abs(diff_cor) > 0.6 & cor_pvals < sig_threshold,
                                            mode = "undirected",
                                            diag = FALSE)

# 5. Pathway enrichment analysis
# ====================================================

perform_enrichment <- function(genes, background) {
  entrez_ids <- mapIds(org.Hs.eg.db, 
                       keys = genes,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  ego <- enrichGO(gene = entrez_ids,
                  universe = background,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.1)
  
  return(ego)
}

# Get background genes (all expressed genes)
background_genes <- rownames(wt_expr)

# Analyze WT modules
wt_modules <- labels2colors(wt_net$colors)
module_enrichment <- list()

for(module in unique(wt_modules)) {
  module_genes <- names(wt_modules)[wt_modules == module]
  ego <- perform_enrichment(module_genes, background_genes)
  module_enrichment[[paste0("WT_", module)]] <- ego
}

# Analyze PD modules
pd_modules <- labels2colors(pd_net$colors)
for(module in unique(pd_modules)) {
  module_genes <- names(pd_modules)[pd_modules == module]
  ego <- perform_enrichment(module_genes, background_genes)
  module_enrichment[[paste0("PD_", module)]] <- ego
}

# 6. Visualization
# ====================================================

# Plot module preservation
plot_preservation <- ggplot(preservation_stats, aes(x=moduleSize, y=Zsummary.pres)) +
  geom_point(aes(color=Zsummary.pres)) +
  scale_color_viridis_c() +
  geom_hline(yintercept = 5, linetype="dashed", color="red") +
  labs(title="Module Preservation Between Conditions",
       x="Module Size", y="Preservation Z-score") +
  theme_minimal()

print(plot_preservation)

# Plot differential co-expression network
plot(diff_network, 
     vertex.size=3,
     vertex.label=NA,
     layout=layout_with_fr,
     main="Differential Co-expression Network")

# Save results
if(writeExcel) {
  write.xlsx(module_enrichment, 
             paste0(output_folder, "/pathway_enrichment_results.xlsx"),
             colNames = TRUE, rowNames = TRUE)
}

saveRDS(list(wt_net = wt_net, pd_net = pd_net, diff_network = diff_network),
        paste0(output_folder, "/coexpression_networks.rds"))



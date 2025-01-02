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

rm(list = ls())
set.seed(18)
options(future.globals.maxSize=10000*1024^2)
writeExcel = TRUE

setwd(getwd())
output_folder <- "./Output"
S <- LoadSeuratRds(file = './Processed/ClusteredUmappedSeurat.rds')

#--------------FIND IDEAL UMAP PARAMETERS--------------

#NO NEED TO RERUN THIS IF YOU IMPORT THE PREPROCESSED SEURAT OBJECT, AS IT IS ALREADY STORED IN THE OBJECT
S = RunUMAP(S,dims = 1:100,assay = "integrated",n.neighbors = 200,seed.use = NULL,min.dist=0.5)

print(
  DimPlot(object = S, reduction = "umap", group.by = "ConditionDay",shuffle = TRUE, label.size = 5,label.box = TRUE,repel=TRUE) + 
    theme_minimal() +ggplot2::theme(legend.position = "none") +
    theme_minimal() +
    labs(x = "UMAP 1", y = "UMAP 2", title = "Umap by Condition and Day") +
    theme(plot.title = element_text(hjust = 0.5))
)

print(
  DimPlot(object = S, reduction = "umap", group.by = "Condition", split.by = "Day",shuffle = TRUE,label.size = 5,label.box = TRUE,repel=TRUE) + 
    theme_minimal() +ggplot2::theme(legend.position = "none") +
    theme_minimal() +
    labs(x = "UMAP 1", y = "UMAP 2", title = "Umap by Condition and Day") +
    theme(plot.title = element_text(hjust = 0.5))
)

FeaturePlot(S,features=c("VIM","TH"))   #find most interesting genes by deg

#--------------FIND BEST CLUSTER RESOLUTION--------------

#NO NEED TO RERUN THIS IF YOU IMPORT THE PREPROCESSED SEURAT OBJECT, AS IT IS ALREADY STORED IN THE OBJECT
DefaultAssay(S) = "integrated"

S = FindNeighbors(S, reduction = "pca", dims = 1:100)
S = FindClusters(S, resolution = 0.17)

DimPlot(S,split.by="ConditionDay",group.by="seurat_clusters")
DimPlot(S,split.by="Day",group.by="seurat_clusters")
DimPlot(S,group.by="seurat_clusters")

#SaveSeuratRds(S,"./Processed/ClusteredUmappedSeurat_copy.rds")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                            ---- DEGs DAYWISE ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# IMPORT

DEG_daywise <- read_excel("./Processed/DEG_daywise/DEG_DETAILED_MADS_new.xlsx")

add_up_down_reg <- function(df1,minpval,minFC){     #add a column to specify if the gene expression is up or down regulated
  df2 <- data.table::copy(df1)
  df2$diffexpressed <- "NOTSIG"
  df2$diffexpressed[df2$avg_log2FC > minFC & df2$p_val_adj < minpval] <- "UP"
  df2$diffexpressed[df2$avg_log2FC < -minFC & df2$p_val_adj < minpval] <- "DOWN"
  return(df2)
}

DEG_daywise <- add_up_down_reg(DEG_daywise,0.05,0.6)  #for daywise DEG´s exp1 is ctl, exp2 is parkin

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                         ---- VOLCANO PLOT DAYWISE ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

vulcano_plot_daywise <- function(df1,minpval,minFC){
  df2 <- data.table::copy(df1)
  df2$Day <- paste("Day",df2$Day,sep = " ")
  df2 <- df2 |> group_by(Day)
  ggplot(data = df2, aes(x = avg_log2FC, y = -log10(p_val_adj),col = diffexpressed)) +
    geom_vline(xintercept = c(-minFC, minFC), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(minpval), col = "gray", linetype = 'dashed') +
    guides(col = guide_legend(title = "exp1 = CTL\nexp2 = PARK")) +
    labs(y = "-log(p_val)", x = "log(FC)",title = "Vulcano plot daywise") +
    geom_point(size = 2) +
    scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"),
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    facet_grid(Day ~ .)
}


vulcano_plot_daywise(DEG_daywise,0.05,0.6)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                         ---- DEGs CLUSTERWISE ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# IMPORT

success = FALSE
try({
  iter = 0
  for (i in unique(S$seurat_clusters)){
    a <- read_excel(paste("./Processed/DEG_clusterwise/",i,".xlsx",sep=""))
    if (iter == 0){
      DEG_clusterwise <- a
    }else{
      DEG_clusterwise <- rbind(a,DEG_clusterwise)
    }
    iter = iter + 1
  }
  success = TRUE}
)
if (success == FALSE){
  iter = 0
  while (iter <= 9){
    a <- read_excel(paste("./Processed/DEG_clusterwise/",iter,".xlsx",sep=""))
    if (iter == 0){
      DEG_clusterwise <- a
    }else{
      DEG_clusterwise <- rbind(a,DEG_clusterwise)
    }
    iter = iter + 1
  }
}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                         ---- VOLCANO PLOT CLUSTERWISE ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

DEG_clusterwise <- add_up_down_reg(DEG_clusterwise,0.05,0.6)

vulcano_plot_clusterwise <- function(df1,minpval,minFC){
  df2 <- data.table::copy(df1)
  df2$ident1 <- paste("Cluster",df2$ident1,sep = " ")
  df2 <- df2|> group_by(ident1)
  ggplot(data = df2, aes(x = avg_log2FC, y = -log10(p_val_adj),col = diffexpressed)) +
    geom_vline(xintercept = c(-minFC, minFC), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(minpval), col = "gray", linetype = 'dashed') +
    geom_point(size = 1) +
    scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"),
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    guides(col = guide_legend(title = "exp1 = indicated cluster\nexp2 = other clusters")) +
    labs(y = "-log(p_val)", x = "log(FC)",title = "Vulcano plot clusterwise") +
    facet_grid(ident1 ~ .)
}


vulcano_plot_clusterwise(DEG_clusterwise[DEG_clusterwise$ident1 < 5,],0.05,0.6)
vulcano_plot_clusterwise(DEG_clusterwise[DEG_clusterwise$ident1 >= 5,],0.05,0.6)

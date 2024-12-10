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
library(readxl)
set.seed(18)
options(future.globals.maxSize=10000*1024^2)
writeExcel = TRUE

setwd(getwd())
output_folder <- "./Output"
S <- LoadSeuratRds(file = './Processed/ClusteredUmappedSeurat.rds')

#--------------FIND IDEAL UMAP PARAMETERS--------------

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

DefaultAssay(S) = "integrated"

S = FindNeighbors(S, reduction = "pca", dims = 1:100)
S = FindClusters(S, resolution = 0.17)

DimPlot(S,split.by="ConditionDay",group.by="seurat_clusters")
DimPlot(S,split.by="Day",group.by="seurat_clusters")
DimPlot(S,group.by="seurat_clusters")

SaveSeuratRds(S,"./Processed/ClusteredUmappedSeurat_copy.rds")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                            ---- DEGs DAYWISE ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# IMPORT

DEG_daywise <- read_excel("./Processed/DEG_daywise/DEG_DETAILED_MADS_new.xlsx")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                         ---- DEGs CLUSTERWISE ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# IMPORT

iter = 0
for (i in unique(S$seurat_clusters)){
  a <- read_excel(paste("./Processed/DEG_clusterwise/",i,".xlsx",sep=""))
  #assign(paste("DEG_cluster_", i,sep=""),a)
  if (iter == 0){
    DEG_clusterwise <- a
  }else{
    DEG_clusterwise <- rbind(a,DEG_clusterwise)
  }
  iter = iter + 1
}

# Creating summarized version of DEG table

for(i in unique(S$seurat_clusters)){
  temp = DEG_clusterwise %>%
    filter(ident1 == i) %>%
    filter(p_val_adj <= 0.05) %>%
    spread(key = ident2,value = avg_log2FC) %>%
    select(gene,setdiff(as.character(c(0:14)),as.character(i))) %>% 
    replace(is.na(.), 0) 
  
  temp = aggregate(x = temp[2:15],               
                   by = list(temp$gene),            
                   FUN = sum)  %>% 
    rename("Group.1" = "Gene")
  
  temp["absFCpower"] = rowMeans(abs(as.matrix(temp[,2:15])))
  temp = temp[order(temp$absFCpower,decreasing = TRUE),]
  
  if(writeExcel){
    write.xlsx(temp,paste(output_folder,"/ClusterWise/Summarized/",i,".xlsx",sep=""), colNames = TRUE, rowNames = FALSE)
  }
}

# Creating list of most discriminating genes

clust_genes = c()

for(i in unique(S$seurat_clusters)){
  temp = Clust_DEG_rep %>%
    filter(ident1 == i) %>%
    filter(p_val_adj <= 0.05) %>%
    spread(key = ident2,value = avg_log2FC) %>%
    select(gene,setdiff(as.character(c(0:14)),as.character(i))) %>% 
    replace(is.na(.), 0) 
  
  temp = aggregate(x = temp[2:15],               
                   by = list(temp$gene),            
                   FUN = sum)  %>% 
    rename("Group.1" = "Gene")
  
  temp["absFCpower"] = rowMeans(abs(as.matrix(temp[,2:15])))
  clust_genes = c(clust_genes,as.vector(temp %>%
                                          filter(absFCpower >= 2) %>%
                                          select(Gene))$Gene)
  
}

clust_genes = unique(clust_genes)

# Correlating Cluster Genes

clust_cor = cor(as.matrix(t(S@assays$RNA_imputed@data[clust_genes,])))
library(corrplot)
order_clust_genes = corrplot(clust_cor,order = "hclust")$corrPos$yName[1:length(clust_genes)]
#order_clust_genes = c(order_clust_genes[60:length(order_clust_genes)],order_clust_genes[1:59])
DefaultAssay(S) = "RNA"

DoHeatmap(S,features = order_clust_genes,group.by = "seurat_clusters")





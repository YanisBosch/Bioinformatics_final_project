library(Seurat)
library(dplyr)
library(tidyverse)
library(viridisLite)
library(viridis)
library(ggplot2)
library(openxlsx)
library(future)
library(slingshot)
set.seed(18)
options(future.globals.maxSize=10000*1024^2)
writeExcel = TRUE

setwd('C:/Users/yanis/Desktop/Universite/Master/Semestre_3/Bioinformatics/Final Project')
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
S = FindClusters(S, resolution = 0.2006)

DimPlot(S,split.by="ConditionDay")
DimPlot(S,split.by="Day")

#SaveSeuratRds(S,"./Processed/ClusteredUmappedSeurat.rds")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                            ---- DEGs DAYWISE ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

DefaultAssay(S) = "RNA"
S = SetIdent(object = S, value = S$ConditionDay)

for(pink in c("PINK1")){
  for(day in c(0,18,25,37,57)) {
    print(paste(pink,day))
    temp = FindMarkers(S, 
                       assay = "RNA", 
                       ident.1 = paste("CTL",day,sep=" "), 
                       ident.2 = paste(pink,day,sep=" "),
                       logfc.threshold = 0.25,
                       min.pct = 0.05,
                       densify = TRUE,
                       test.use = "MAST",
                       latent.vars = "nFeature_RNA",
                       slot = "data")
    temp$Day = day
    temp$Gene = rownames(temp)
    temp$ident.1 = "CTL"
    temp$ident.2 = pink
    rownames(temp) = paste(pink,day,c(1:dim(temp)[1]),sep="_")
    if(day==0 & pink == "PINK1"){
      DEG = temp
    }else{
      DEG = rbind(DEG,temp)
    }
  }
}

DEG = DEG %>% 
  select(Gene,Day,pct.1,pct.2,avg_log2FC,p_val,p_val_adj,ident.1,ident.2)  %>%
  arrange(Day,Gene,ident.2)

# Export

#DEG_RNA <- DEG
#DEG_RNAimp <- DEG

if(writeExcel){
  write.xlsx(DEG,paste(output_folder,"DEG_DETAILED_MADS_new.xlsx",sep=""),colNames = TRUE,rowNames = FALSE)
}

# Test the effect of the corresponding Assay
DefaultAssay(S) = "RNA_imputed"
#DefaultAssay(S) = "RNA"

FeaturePlot(S,features=c("ABHD13", "VIM", "TH"))


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                         ---- DEGs CLUSTERWISE ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

DefaultAssay(S) = "RNA"
S = SetIdent(object = S, value = S$seurat_clusters)
pairwise_clust = combn(unique(S$seurat_clusters),2)

for(idx in 1:dim(pairwise_clust)[2]) {
  print(idx)
  pw = pairwise_clust[,idx]
  temp = FindMarkers(S, 
                     assay = "RNA", 
                     ident.1 = pw[1],
                     ident.2 = pw[2],
                     logfc.threshold = 0.25,
                     min.pct = 0.1,
                     densify = TRUE,
                     test.use = "MAST",
                     latent.vars = "nFeature_RNA"
  )
  temp$gene = rownames(temp)
  temp$pct_diff = abs(temp$pct.1-temp$pct.2)
  temp$ident1 = pw[1]
  temp$ident2 = pw[2]
  if(idx==1){
    Clust_DEG = temp
  }else{
    Clust_DEG = rbind(Clust_DEG,temp)
  }
}

# order cluster along the differentiation trajectory
Clust_DEG = Clust_DEG[,c(6,1,5,2,3,4,7,8,9)]     #ADAPT TO OUR CLUSTERS

# Duplicating tests to ease search on excel file

Clust_DEG_rep = Clust_DEG
Clust_DEG_rep$avg_log2FC = -Clust_DEG_rep$avg_log2FC
temp = Clust_DEG_rep$ident2
Clust_DEG_rep$ident2 = Clust_DEG_rep$ident1
Clust_DEG_rep$ident1 = temp
temp = Clust_DEG_rep$pct.2
Clust_DEG_rep$pct.2 = Clust_DEG_rep$pct.1
Clust_DEG_rep$pct.1 = temp
Clust_DEG_rep = rbind(Clust_DEG,Clust_DEG_rep)

# Export
for(i in unique(S$seurat_clusters)){
  temp = Clust_DEG_rep %>%
    filter(ident1 == i)
  if(writeExcel){
    write.xlsx(temp,paste(output_folder,"/ClusterWise/Detailed/",i,".xlsx",sep=""), colNames = TRUE, rowNames = FALSE)
  }
}

# Creating summarized version of DEG table

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




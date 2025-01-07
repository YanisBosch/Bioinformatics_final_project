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
options(future.globals.maxSize=10000*1024^2)        #required for loading large seurat into ram
writeExcel = TRUE

setwd(getwd())
output_folder <- "./Processed"

#--------------PLEASE RUN THE PREPROCESSING SCRIPT IF YOU DO NOT HAVE THIS FILE ALREADY--------------

S <- LoadSeuratRds(file = './Processed/SeuratFinal.rds')

#--------------FIND IDEAL UMAP PARAMETERS--------------

#NO NEED TO RERUN THIS IF YOU IMPORT THE PREPROCESSED SEURAT OBJECT, AS IT IS ALREADY STORED IN THE OBJECT
S = RunUMAP(S,dims = 1:100,assay = "integrated",n.neighbors = 200,seed.use = NULL,min.dist=0.5)
#--------------

#--------------Plots--------------

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
#--------------

#--------------Plots--------------

S = FindNeighbors(S, reduction = "pca", dims = 1:100)
S = FindClusters(S, resolution = 0.17)        #Resolution was chosen to yield 10 clusters
                                              #as we have ten combinations of day/condition

DimPlot(S,split.by="ConditionDay",group.by="seurat_clusters")
DimPlot(S,split.by="Day",group.by="seurat_clusters")
DimPlot(S,group.by="seurat_clusters")

#ONLY NEEDS TO BE RUN ONCE AFTER COMPUTING THE CLUSTERS AND UMAP REPRESENTATION
SaveSeuratRds(S,"./Processed/SeuratFinal.rds")
#--------------

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                            ---- DEGs DAYWISE ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

#--------------RUN THIS CODE IF YOU DO NOT HAVE THE PREPROCESSED DATA ALREADY--------------

#COMPUTATION

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

if(writeExcel){
  write.xlsx(DEG,paste(output_folder,"/DEG_daywise/DEG_DETAILED_MADS_new.xlsx",sep=""),colNames = TRUE,rowNames = FALSE)
}

DEG_daywise <- DEG

#--------------RUN THIS CODE IF YOU ALREADY HAVE THE PREPROCESSED DATA--------------

# IMPORT

DEG_daywise <- read_excel("./Processed/DEG_daywise/DEG_DETAILED_MADS_new.xlsx")

#--------------ALWAYS RUN THIS CODE--------------

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
  #df1 = dataframe containing daywise DEG data
  #minpval = p value threshold
  #minFC = fold change threshold
  df2 <- data.table::copy(df1)      #create a copy of the supplied dataframe
                                    #to avoid mistakenly modifying it
  df2$Day <- paste("Day",df2$Day,sep = " ")
                                    #replace the numbers by "Day " + number for cleaner plotting
  df2 <- df2 |> group_by(Day)       #group by day for faceting
  ggplot(data = df2, aes(x = avg_log2FC, y = -log10(p_val_adj),col = diffexpressed)) +
    geom_vline(xintercept = c(-minFC, minFC), col = "gray", linetype = 'dashed') +
    #add vertical lines indicating the fold change threshold
    geom_hline(yintercept = -log10(minpval), col = "gray", linetype = 'dashed') +
    #add horizontal lines indicating the p-value threshold
    guides(col = guide_legend(title = "exp1 = CTL\nexp2 = PARK")) +
    #change title of legend
    labs(y = "-log(p_val)", x = "log(FC)",title = "Vulcano plot daywise") +
    #label axes properly and set title
    geom_point(size = 2) +
    scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"),
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    #colour points based the column in df2 that indicates up or down regulation
    facet_grid(Day ~ .)
    #facet by day
}


vulcano_plot_daywise(DEG_daywise,0.05,0.6)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                         ---- DEGs CLUSTERWISE ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

#--------------RUN THIS CODE IF YOU DO NOT HAVE THE PREPROCESSED DATA ALREADY--------------

#COMPUTATION

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
Clust_DEG = Clust_DEG[,c(6,1,8,2,5,3,9,0,7,4)] 

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
    write.xlsx(temp,paste(output_folder,"/DEG_clusterwise/",i,".xlsx",sep=""), colNames = TRUE, rowNames = FALSE)
  }
}

#--------------RUN THIS CODE IF YOU ALREADY HAVE THE PREPROCESSED DATA--------------

# IMPORT

success = FALSE
try({                         #try to get cluster numbers from Seurat object 
  iter = 0
  for (i in unique(S$seurat_clusters)){   #read all .xlsx files in succession
    a <- read_excel(paste("./Processed/DEG_clusterwise/",i,".xlsx",sep=""))
    if (iter == 0){
      DEG_clusterwise <- a    
    }else{
      DEG_clusterwise <- rbind(a,DEG_clusterwise)
      #build dataframe iteratively by joining the individual parts
    }
    iter = iter + 1
  }
  success = TRUE}
)
if (success == FALSE){        #if the Seurat object is not loaded in memory just
                              #manually input the cluster numbers
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
  #df1 = dataframe containing daywise DEG data
  #minpval = p value threshold
  #minFC = fold change threshold
  df2 <- data.table::copy(df1)      #create a copy of the supplied dataframe
                                    #to avoid mistakenly modifying it
  df2$ident1 <- paste("Cluster",df2$ident1,sep = " ")
                                    #replace the numbers by "Cluster " + number for cleaner plotting
  df2 <- df2|> group_by(ident1)     ##group by cluster for faceting
  ggplot(data = df2, aes(x = avg_log2FC, y = -log10(p_val_adj),col = diffexpressed)) +
    geom_vline(xintercept = c(-minFC, minFC), col = "gray", linetype = 'dashed') +
    #add vertical lines indicating the fold change threshold
    geom_hline(yintercept = -log10(minpval), col = "gray", linetype = 'dashed') +
    #add horizontal lines indicating the p-value threshold
    geom_point(size = 2) +
    scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"),
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    #colour points based the column in df2 that indicates up or down regulation
    guides(col = guide_legend(title = "exp1 = indicated cluster\nexp2 = other clusters")) +
    #change title of legend
    labs(y = "-log(p_val)", x = "log(FC)",title = "Vulcano plot clusterwise") +
    #label axes properly and set title
    facet_grid(ident1 ~ .)
    #facet by day
}


vulcano_plot_clusterwise(DEG_clusterwise[DEG_clusterwise$ident1 < 5,],0.05,0.6)
vulcano_plot_clusterwise(DEG_clusterwise[DEG_clusterwise$ident1 >= 5,],0.05,0.6)

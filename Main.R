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

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                         ---- VOLCANO PLOT DAYWISE ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

DEG_clusterwise <- add_up_down_reg(DEG_clusterwise,0.05,0.6)

vulcano_plot_clusterwise <- function(df1,minpval,minFC){
  df2 <- data.table::copy(df1)
  df2$ident1 <- paste("Cluster",df2$ident1,sep = " ")
  df2 <- df2|> group_by(ident1)
  ggplot(data = df2, aes(x = avg_log2FC, y = -log10(p_val_adj),col = diffexpressed)) +
    geom_vline(xintercept = c(-minFC, minFC), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(minpval), col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"),
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    facet_grid(ident1 ~ .)
}


vulcano_plot_clusterwise(DEG_clusterwise[DEG_clusterwise$ident1 < 5,],0.05,0.6)
vulcano_plot_clusterwise(DEG_clusterwise[DEG_clusterwise$ident1 >= 5,],0.05,0.6)


# order cluster along the differentiation trajectory
Clust_DEG = DEG_clusterwise[,c(1,6,4,0,9,7,5,3,2,8)] 

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

# Creating list of most discriminating genes

clust_genes = c()

for(i in unique(S$seurat_clusters)){
  temp = Clust_DEG_rep %>%
    filter(ident1 == i) %>%
    filter(p_val_adj <= 0.05) %>%
    spread(key = ident2,value = avg_log2FC) %>%
    select(gene,setdiff(as.character(c(0:9)),as.character(i))) %>% 
    replace(is.na(.), 0) 
  
  temp = aggregate(x = temp[2:10],               
                   by = list(temp$gene),            
                   FUN = sum)  %>% 
    rename("Group.1" = "Gene")
  
  temp["absFCpower"] = rowMeans(abs(as.matrix(temp[,2:10])))
  clust_genes = c(clust_genes,as.vector(temp %>%
                                          filter(absFCpower >= 2) %>%
                                          select(Gene))$Gene)
}

clust_genes = unique(clust_genes)

# Correlating Cluster Genes

clust_cor = cor(as.matrix(t(S@assays$RNA_imputed@data[clust_genes,])))
order_clust_genes = corrplot(clust_cor,order = "hclust")$corrPos$yName[1:length(clust_genes)]
#order_clust_genes = c(order_clust_genes[60:length(order_clust_genes)],order_clust_genes[1:59])
DefaultAssay(S) = "RNA"

DoHeatmap(S,features = order_clust_genes,group.by = "seurat_clusters")


#--------------CORRELATION--------------

#--------------Filtering for high correlation, split by ctrl and pink

tempa <- subset(S, ConditionDay %in% "CTL 0")
expression_matrix_CTL_0 <- GetAssayData(tempa, slot = "data")

tempa <- subset(S, ConditionDay %in% "PINK1 0")
expression_matrix_PINK_0 <- GetAssayData(tempa, slot = "data")

#for control condition at day 0
rcorr_result_ctr_0 <- rcorr(t(as.matrix(expression_matrix_CTL_0)))
cor_matrix_ctr_0 <- rcorr_result_ctr_0$r
#for pink1 condition day 0
rcorr_result_pink_0 <- rcorr(t(as.matrix(expression_matrix_PINK_0)))
cor_matrix_pink_0 <- rcorr_result_pink_0$r

#replaces NaN with 0
cor_matrix_ctr_0[is.na(cor_matrix_ctr_0)] <- 0
cor_matrix_pink_0[is.na(cor_matrix_pink_0)] <- 0

hist(cor_matrix_ctr_0,breaks=200)
hist(cor_matrix_pink_0,breaks=200)

# to reduce correlation matrix by identifying rows and columns with at least one correlation > cthresh
# play with treshhold
ctresh <- 0.6
# for filtering have to set trivial diagonal to 0
diag(cor_matrix_ctr_0) <- 0
rows_to_keep_ctr_0 <- apply(cor_matrix_ctr_0, 1, function(row) any(abs(row) > ctresh))
cols_to_keep_ctr_0 <- apply(cor_matrix_ctr_0, 2, function(col) any(abs(col) > ctresh))

diag(cor_matrix_pink_0) <- 0
rows_to_keep_pink_0 <- apply(cor_matrix_pink_0, 1, function(row) any(abs(row) > ctresh))
cols_to_keep_pink_0 <- apply(cor_matrix_pink_0, 2, function(col) any(abs(col) > ctresh))

# Subset the correlation matrix
reduced_cor_matrix_ctr_0 <- cor_matrix_ctr_0[rows_to_keep_ctr_0, cols_to_keep_ctr_0]
reduced_cor_matrix_pink_0 <- cor_matrix_pink_0[rows_to_keep_pink_0, cols_to_keep_pink_0]

#--------------Plotting

#for CTRL condition
filtered_cor_matrix_ctr_0 <- ifelse(abs(reduced_cor_matrix_ctr_0) >= ctresh, reduced_cor_matrix_ctr_0, 0)

heatmap(filtered_cor_matrix_ctr_0, na.rm=T, Rowv = NA, Colv = NA,
        col = colorRampPalette(c("blue", "white", "red"))(50),
        scale = "none")

dim(filtered_cor_matrix_ctr_0 )

#for PINK1 condition
filtered_cor_matrix_pink_0 <- ifelse(abs(reduced_cor_matrix_pink_0) >= ctresh_2, reduced_cor_matrix_pink_0, 0)

heatmap(filtered_cor_matrix_pink_0, na.rm=T, Rowv = NA, Colv = NA,
        col = colorRampPalette(c("blue", "white", "red"))(50),
        scale = "none")

dim(filtered_cor_matrix_pink_0 )

# adjacency matrix 
#ctr
adj_matrix_ctr_0 <- ifelse(abs(reduced_cor_matrix_ctr_0) > ctresh, 1, 0) #how to break this down into pos and negative
diag(adj_matrix_ctr_0) <- 0

heatmap(adj_matrix_ctr_0, na.rm=T, Rowv = NA, Colv = NA,
        col = colorRampPalette(c("white", "black"))(3),
        scale = "none")

nodes_ctr_0 <-  rowSums(adj_matrix_ctr_0)
hist(nodes_ctr_0, breaks=100)

#pink1
adj_matrix_pink_0 <- ifelse(abs(reduced_cor_matrix_pink_0) > ctresh, 1, 0) #how to break this down into pos and negative
diag(adj_matrix_pink_0) <- 0

heatmap(adj_matrix_pink_0, na.rm=T, Rowv = NA, Colv = NA,
        col = colorRampPalette(c("white", "black"))(3),
        scale = "none")

nodes_pink_0 <-  rowSums(adj_matrix_pink_0)
hist(nodes_pink_0, breaks=100)

#--------------NETWORKS--------------

# ctr
network_graph_ctr_0 <- graph_from_adjacency_matrix(adj_matrix_ctr_0, mode = "undirected")

plot(network_graph_ctr_0, 
     vertex.size = 5,        # Size of nodes
     vertex.label.cex = 0.5,  # Font size for labels
     vertex.label.color = "black",  # Label color
     edge.width = 2,          # Width of edges
     edge.color = "gray")     # Color of edges


V(network_graph_ctr_0)$color <- c(colorRampPalette(c("black", "grey", "orange", "red"))(dim(adj_matrix_ctr_0)[1]))
E(network_graph_ctr_0)$width <- 2

plot(network_graph_ctr_0,
     vertex.size = 5,        # Size of nodes
     vertex.label.cex = .6,  # Font size for labels
     vertex.label.color = "black",  # Label color
     edge.width = 4,          # Width of edges
     edge.color = "gray"
)

# pink
network_graph_pink_0 <- graph_from_adjacency_matrix(adj_matrix_pink_0, mode = "undirected")

plot(network_graph_pink_0, 
     vertex.size = 5,        # Size of nodes
     vertex.label.cex = 0.5,  # Font size for labels
     vertex.label.color = "black",  # Label color
     edge.width = 2,          # Width of edges
     edge.color = "gray")     # Color of edges


V(network_graph_pink_0)$color <- c(colorRampPalette(c("black", "grey", "orange", "red"))(dim(adj_matrix_pink_0)[1]))
E(network_graph_pink_0)$width <- 2

plot(network_graph_pink_0,
     vertex.size = 5,        # Size of nodes
     vertex.label.cex = .6,  # Font size for labels
     vertex.label.color = "black",  # Label color
     edge.width = 4,          # Width of edges
     edge.color = "gray"
)


layout(matrix(c(1, 2), nrow = 1, ncol = 2, byrow = TRUE))

# reduce adjacency matrix to only highly connected nodes
ntres = 5
nodes_to_keep_ctr_0 <- which(nodes_ctr_0 >= ntres)
red_adj_matrix_ctr_0 <- adj_matrix_ctr_0[nodes_to_keep_ctr_0, nodes_to_keep_ctr_0]

nodes_to_keep_pink_0 <- which(nodes_pink_0 >= ntres)
red_adj_matrix_pink_0 <- adj_matrix_pink_0[nodes_to_keep_pink_0, nodes_to_keep_pink_0]


heatmap(red_adj_matrix_ctr_0, na.rm=T, Rowv = NA, Colv = NA,
        col = colorRampPalette(c("white", "black"))(3),
        scale = "none")

heatmap(red_adj_matrix_pink_0, na.rm=T, Rowv = NA, Colv = NA,
        col = colorRampPalette(c("white", "black"))(3),
        scale = "none")

# for network representation

#ctr
red_network_graph_ctr_0 <- graph_from_adjacency_matrix(red_adj_matrix_ctr_0, mode = "undirected")

plot(red_network_graph_ctr_0, 
     vertex.size = 5,        # Size of nodes
     vertex.label.cex = 0.5,  # Font size for labels
     vertex.label.color = "black",  # Label color
     edge.width = 2,          # Width of edges
     edge.color = "gray")     # Color of edges

#pink1
red_network_graph_pink_0 <- graph_from_adjacency_matrix(red_adj_matrix_pink_0, mode = "undirected")

plot(red_network_graph_pink_0, 
     vertex.size = 5,        # Size of nodes
     vertex.label.cex = 0.5,  # Font size for labels
     vertex.label.color = "black",  # Label color
     edge.width = 2,          # Width of edges
     edge.color = "gray")     # Color of edges


# to compare condition by adjacency pattern

# Subset the correlation matrix
com_ind_0 <- c(which(rows_to_keep_ctr_0),which(rows_to_keep_pink_0)) 
shared_red_cor_matrix_ctr_0 <- cor_matrix_ctr_0[com_ind_0, com_ind_0  ]
shared_red_cor_matrix_pink_0 <- cor_matrix_pink_0[com_ind_0, com_ind_0  ]


h1 <- Heatmap(shared_red_cor_matrix_ctr_0, border=T,
              name = "Correlation", cluster_rows = FALSE, cluster_columns = FALSE, 
              show_row_names = FALSE, show_column_names = FALSE,
              col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

h2 <- Heatmap(shared_red_cor_matrix_pink_0,  border=T,
              name = "Correlation", cluster_rows = FALSE, cluster_columns = FALSE, 
              show_row_names = FALSE, show_column_names = FALSE,
              col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))


draw(h1 + h2)

# for common adjacency as above

com_adj_matrix_ctr_0 <- ifelse(abs(shared_red_cor_matrix_ctr_0) > ctresh, 1, 0) #how to break this down into pos and negative
diag(com_adj_matrix_ctr_0) <- 0

heatmap(com_adj_matrix_ctr_0, na.rm=T, Rowv = NA, Colv = NA,
        col = colorRampPalette(c("white", "black"))(3),
        scale = "none")

#pink1
com_adj_matrix_pink_0 <- ifelse(abs(shared_red_cor_matrix_pink_0) > ctresh, 1, 0) #how to break this down into pos and negative
diag(com_adj_matrix_pink_0) <- 0

heatmap(com_adj_matrix_pink_0, na.rm=T, Rowv = NA, Colv = NA,
        col = colorRampPalette(c("white", "black"))(3),
        scale = "none")


h3 <- Heatmap(com_adj_matrix_ctr_0,  border=T,
              name = "Adjacency CTR", cluster_rows = FALSE, cluster_columns = FALSE, 
              show_row_names = FALSE, show_column_names = FALSE,
              col = colorRamp2(c( 0, 1), c( "white", "black")))

h4 <- Heatmap(com_adj_matrix_pink_0,  border=T,
              name = "Adjecency PINK1", cluster_rows = FALSE, cluster_columns = FALSE, 
              show_row_names = FALSE, show_column_names = FALSE,
              col = colorRamp2(c(0, 1), c( "white", "black")))


draw(h3 + h4)

# can be done for other time points - separate?


##### for big correlation calculations:

# Install gputools package
#install.packages("gputools")

# Load the package
#library(gputools)

# Compute correlation matrix
#cor_matrix <- gpuCor(data_matrix, method = "pearson")

# or alternative
# Install package
install.packages("bigstatsr")
# Load the package
library(bigstatsr)
# to use binary format
cm_0_c <- as_FBM(as.matrix(expression_matrix_CTL_0)) # binary expression matrix for cell correlations
cor_matrix_cells <- big_cor(cm_0_c) # correlation for cells


########## For String analysis

degs <- read.xlsx('./Results/DEG_DETAILED_MADS_new_RNA.xlsx')

degs_0_all <- degs[degs$Day==0,]

degs_0 <- degs_0_all[abs(degs_0_all$avg_log2FC)>1,1] 

write_csv(as.data.frame(degs_0), file='./degs_day0.csv')

#this csv file can be then uploaded into string 








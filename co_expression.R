# ================================= LIBRARIES ==================================

library(Seurat)
library(dplyr)
library(tidyverse)
library(viridisLite)
library(viridis)
library(ggplot2)
library(openxlsx)
library(slingshot)

library(Hmisc)
library(ComplexHeatmap)
library(circlize)
library(igraph)

# to get clean environment
rm(list = ls())


# ================================ PARAMETERS ==================================

# set your working directory accordingly
setwd(getwd())

input_file = "./Preprocessed/SeuratFinal.rds"
output_folder = "./Results/co_expression"

#to export things later
writeExcel = TRUE
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                             ---- Gene co expression Analysis  ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
gene_co_exp_day <- function(day){
  
  
  # Load mapped data-sets
  
  S = readRDS(input_file)
  
  
  # Convert meta-data to factors
  
  S$Day = factor(S$Day, levels = c(0,18,25,37,57))
  S$Condition = str_replace_all(S$Condition, "ND","PINK1")
  S$Condition = str_replace_all(S$Condition,"WT","CTL")
  S$ConditionDay = paste(S$Condition,S$Day,sep = " ")
  
  S$ConditionDay = factor(S$ConditionDay, levels = c("CTL 0","PINK1 0",
                                                     "CTL 18","PINK1 18",
                                                     "CTL 25","PINK1 25",
                                                     "CTL 37","PINK1 37",
                                                     "CTL 57","PINK1 57"))
  
  S$ConditionDay = factor(S$ConditionDay, levels = c("CTL 0","CTL 18","CTL 25","CTL 37","CTL 57",
                                                     "PINK1 0","PINK1 18","PINK1 25","PINK1 37","PINK1 57"))
  
  S$Condition = factor(S$Condition, levels = c("CTL","PINK1"))
  
  
  # Create output directory
  output_dir <- file.path(output_folder, paste0("day_", as.character(day)))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  calculate_cor_mat<- function(condition) {
    
    condition_day <- paste(condition, day)
    print(paste("Filtering for:", condition_day))
    
    ## Step 1: Extract expression matrix for each condition and day
    print(dim(S))
    
    tempa <- subset(S, ConditionDay %in% condition_day)
    expression_matrix <- GetAssayData(tempa, layer = "data")
    
    ## Step 2: Compute Correlation matrix and convert NaN to 0
    
    ###for control condition at day 
    rcorr_result <- rcorr(t(as.matrix(expression_matrix)))
    cor_matrix <- rcorr_result$r
    
    ###replaces NaN with 0
    cor_matrix[is.na(cor_matrix)] <- 0
    cat("The correlation matrix for", condition, "has dimensions:", dim(cor_matrix), "\n")
    
    ## Step 3: Filter the correlation matrix by threshold and diagonal to 0
    # reduce correlation matrix by identifying rows and columns with at least one correlation > cthresh
    ctresh <- 0.8 # treshhold
    # for filtering have to set trivial diagonal to 0
    diag(cor_matrix) <- 0
    rows_to_keep <- apply(cor_matrix, 1, function(row) any(abs(row) > ctresh))
    cols_to_keep <- apply(cor_matrix, 2, function(col) any(abs(col) > ctresh))
    
    # Subset the correlation matrix
    reduced_cor_matrix <- cor_matrix[rows_to_keep, cols_to_keep]
    cat("The reduced correlation matrix for", condition, "has dimensions:", dim(reduced_cor_matrix), "\n")
    
    
    
    ### Step 4: Filter the correlation matrix by threshold and diagonal to 0
    
    ctresh_2 = 0.8
    #for CTRL condition
    filtered_cor_matrix <- ifelse(abs(reduced_cor_matrix) >= ctresh_2, reduced_cor_matrix, 0)
    
    cat("The filtered cor matrix for", condition, "has dimensions:", dim(filtered_cor_matrix), "\n")
    
    ### Step 5: adjacency matrix 
    
    adj_matrix <- ifelse(abs(reduced_cor_matrix) > ctresh, 1, 0)
    diag(adj_matrix) <- 0
    cat("The adj matrix for", condition, "has dimensions:", dim(adj_matrix), "\n")
    ### reduce adjacency matrix to only highly connected nodes
    nodes <- rowSums(adj_matrix)
    ntres = 5
    nodes_to_keep <- which(nodes >= ntres)
    red_adj_matrix <- adj_matrix[nodes_to_keep, nodes_to_keep]
    
    cat("The high cor matrix for", condition, "has dimensions:", dim(red_adj_matrix), "\n")
    return(
      list(
            cor_matrix = cor_matrix, 
            rows_to_keep = rows_to_keep,
            cols_to_keep = cols_to_keep,
            reduced_cor_matrix = reduced_cor_matrix, 
            filtered_cor_matrix=filtered_cor_matrix,  
            adj_matrix = adj_matrix, 
            nodes= nodes, 
            red_adj_matrix = red_adj_matrix))
  }
  ctl <- calculate_cor_mat("CTL")
  pink1 <- calculate_cor_mat("PINK1")
  dim(pink1$cor_matrix)
  gc()
  
  #############################correlation matrix##############################################
  # Specify the folder path
  folder_path <- file.path(output_dir, "cor_mat")
  dir_of_folder_path <- dirname(folder_path)
  print(dir_of_folder_path)
  
  # Create the folder
  dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  
  cor_matrix_ctr<- ctl$cor_matrix
  cor_matrix_pink1<- pink1$cor_matrix
  
  # Ensure the output directory exists
  output_subdir <- paste0(output_dir, "/cor_mat/")
  dir.create(output_subdir, showWarnings = FALSE, recursive = TRUE)
  

  
  # Save histogram with larger PNG dimensions
  png(
    filename = paste0(output_subdir, "/hist_cor_mat_ctr_", day, "_.png"), 
    res = 300, 
    width = 2400,   # 8 inches at 300 dpi
    height = 1800   # 6 inches at 300 dpi
  )
  hist(cor_matrix_ctr, 
       breaks = 200, 
       main = paste0("Distribution of Correlation Values for CTL", day),
       xlab = "Correlation", 
       ylab = "Frequency")
  dev.off()
  
  # Save PINK1 histogram with larger dimensions
  png(
    filename = paste0(output_subdir, "/hist_cor_mat_pink1_", day, "_.png"), 
    res = 300, 
    width = 2400,   # 8 inches at 300 dpi
    height = 1800   # 6 inches at 300 dpi
  )
  hist(cor_matrix_pink1, 
       breaks = 200, 
       main = paste0("Distribution of Correlation Values for PINK1 ", day),  
       xlab = "Correlation", 
       ylab = "Frequency")
  dev.off()
  
  # Common parameters for all PNG devices
  png_width <- 2400   # 8 inches at 300 dpi
  png_height <- 2400  # 8 inches at 300 dpi (square for symmetric matrices)
  png_res <- 300
  
  ##-----------------------------------------
  # Base R Heatmaps (with Clustering)
  #-----------------------------------------
  png(
    paste0(output_subdir, "/heatmap_cor_mat_ctl_", day, "_.png"),
    width = png_width, height = png_height, res = png_res
  )
  heatmap(
    cor_matrix_ctr,
    na.rm = TRUE,
    col = colorRampPalette(c("blue", "white", "red"))(50),
    scale = "none",
    main = paste0("Heatmap of Correlation Matrix for CTL at day ", day)
  )
  dev.off()
  
  png(
    paste0(output_subdir, "/heatmap_cor_mat_pink1_", day, "_.png"),
    width = png_width, height = png_height, res = png_res
  )
  heatmap(
    cor_matrix_pink1,
    na.rm = TRUE,
    col = colorRampPalette(c("blue", "white", "red"))(50),
    scale = "none",
    main = paste0("Heatmap of Correlation Matrix for PINK1 at day ", day)
  )
  dev.off()
  
  
  #-----------------------------------------
  # Base R Heatmaps (Without Clustering)
  #-----------------------------------------
  png(
    paste0(output_subdir, "/heatmap_cor_mat_no_clustering_ctr_", day, "_.png"),
    width = png_width, height = png_height, res = png_res
  )
  heatmap(
    cor_matrix_ctr,
    na.rm = TRUE,
    Rowv = NA, Colv = NA,  # Disable clustering
    col = colorRampPalette(c("blue", "white", "red"))(50),
    scale = "none",
    main = paste0("Heatmap without Clustering for CTL at day ", day)
  )
  dev.off()
  
  
  png(
    paste0(output_subdir, "/heatmap_cor_mat_no_clustering_pink_", day, "_.png"),
    width = png_width, height = png_height, res = png_res
  )
  heatmap(
    cor_matrix_pink1,
    na.rm = TRUE,
    Rowv = NA, Colv = NA,  # Disable clustering
    col = colorRampPalette(c("blue", "white", "red"))(50),
    scale = "none",
    main = paste0("Heatmap without Clustering for PINK1 at day ", day)
  )
  dev.off()
  
  
  #-----------------------------------------
  # ComplexHeatmap Plots (Clustered)
  #-----------------------------------------
  library(ComplexHeatmap)  # Ensure package is loaded
  
  png(
    paste0(output_subdir, "/complex_cor_mat_ctl_", day, "_.png"),
    width = png_width, height = png_height, res = png_res
  )
  draw(  # Explicitly draw the Heatmap object
    Heatmap(
      cor_matrix_ctr,
      name = "Correlation",
      col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    )
  )
  dev.off()
  
  png(
    paste0(output_subdir, "/complex_cor_mat_pink_", day, "_.png"),
    width = png_width, height = png_height, res = png_res
  )
  draw(
    Heatmap(
      cor_matrix_pink1,
      name = "Correlation",
      col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    )
  )
  dev.off()
  
  #-----------------------------------------
  # ComplexHeatmap Plots (Without Clustering)
  #-----------------------------------------
  png(
    paste0(output_subdir, "/complex_cor_mat_no_clustering_ctr_", day, "_.png"),
    width = png_width, height = png_height, res = png_res
  )
  draw(
    Heatmap(
      cor_matrix_ctr,
      name = "Correlation",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    )
  )
  dev.off()
  
  png(
    paste0(output_subdir, "/complex_cor_mat_no_clustering_pink_", day, "_.png"),
    width = png_width, height = png_height, res = png_res
  )
  draw(
    Heatmap(
      cor_matrix_pink1,
      name = "Correlation",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    )
  )
  dev.off()
  
  
  #############################reduced correlation matrix##############################################
 
  folder_path <- file.path(output_dir, "red_cor_mat")
  
  # Create the folder
  dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  
  
  reduced_cor_matrix_ctr <- ctl$reduced_cor_matrix
  reduced_cor_matrix_pink1 <- pink1$reduced_cor_matrix
  
  #----------------------------------------------------------
  # ComplexHeatmap (No Clustering, No Row/Column Names) 
  #----------------------------------------------------------
  png(
    paste0(output_dir, "/red_cor_mat/heatmap_reduced_cor_ctr_", day, "_no_names_clustering.png"),
    width = png_width, height = png_height, res = png_res
  )
  draw(  # Use draw() to render ComplexHeatmap objects
    Heatmap(
      reduced_cor_matrix_ctr,
      name = "Correlation",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    )
  )
  dev.off()
  
  
  png(
    paste0(output_dir, "/red_cor_mat/heatmap_reduced_cor_pink1_", day, "_no_names_clustering.png"),
    width = png_width, height = png_height, res = png_res
  )
  draw(
    Heatmap(
      reduced_cor_matrix_pink1,
      name = "Correlation",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    )
  )
  dev.off()
  
  
  
  #----------------------------------------------------------
  # Base R Heatmap (No Clustering) 
  #----------------------------------------------------------
  
  png(
    paste0(output_dir, "/red_cor_mat/heatmap_reduced_cor_ctr_", day, "_no_clustering.png"),
    width = png_width, height = png_height, res = png_res
  )
  heatmap(
    reduced_cor_matrix_ctr,
    na.rm = TRUE,
    Rowv = NA,  # Disable row clustering
    Colv = NA,  # Disable column clustering
    col = colorRampPalette(c("blue", "white", "red"))(50),
    scale = "none"
  )
  dev.off()
  
  
  png(
    paste0(output_dir, "/red_cor_mat/heatmap_reduced_cor_pink1_", day, "_no_clustering.png"),
    width = png_width, height = png_height, res = png_res
  )
  heatmap(
    reduced_cor_matrix_pink1,
    na.rm = TRUE,
    Rowv = NA,  # Disable row clustering
    Colv = NA,  # Disable column clustering
    col = colorRampPalette(c("blue", "white", "red"))(50),
    scale = "none"
  )
  dev.off()
  

  
  ##############################Filtering the Correlation Matrices by threshold############################
  
  # Specify the folder path
  folder_path <- file.path(output_dir, "fil_cor_mat")
  
  # Create the folder
  dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  
  filtered_cor_matrix_ctr <- ctl$filtered_cor_matrix
  filtered_cor_matrix_pink1 <- pink1$filtered_cor_matrix
  
  #----------------------------------------------------------
  # Filtered Heatmap (Base R)
  #----------------------------------------------------------
  png(
    paste0(output_dir, "/fil_cor_mat/heatmap_filtered_cor_ctr_", day, "_.png"),
    width = png_width, 
    height = png_height, 
    res = png_res
  )
  heatmap(
    filtered_cor_matrix_ctr,
    na.rm = TRUE,
    Rowv = NA,  # Disable row clustering
    Colv = NA,  # Disable column clustering
    col = colorRampPalette(c("blue", "white", "red"))(50),
    scale = "none"
  )
  dev.off()
  
  png(
    paste0(output_dir, "/fil_cor_mat/heatmap_filtered_cor_pink1_", day, "_.png"),
    width = png_width, 
    height = png_height, 
    res = png_res
  )
  heatmap(
    filtered_cor_matrix_pink1,
    na.rm = TRUE,
    Rowv = NA,
    Colv = NA,
    col = colorRampPalette(c("blue", "white", "red"))(50),
    scale = "none"
  )
  dev.off()
  
  ##################### Adjacency Matrix #################################################################
  # Specify the folder path
  folder_path <- file.path(output_dir, "adj_cor_mat")
  
  # Create the folder
  dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  
  adj_matrix_ctr<-ctl$adj_matrix
  adj_matrix_pink1<-pink1$adj_matrix
  
  #----------------------------------------------------------
  # Adjacency Matrix Heatmap 
  #----------------------------------------------------------
  png(
    paste0(output_dir, "/adj_cor_mat/heatmap_adj_ctr_", day, "_.png"),
    width = png_width, 
    height = png_height, 
    res = png_res
  )
  heatmap(
    adj_matrix_ctr,
    na.rm = TRUE,
    Rowv = NA,  # No row clustering
    Colv = NA,  # No column clustering
    col = colorRampPalette(c("white", "black"))(3),  # Binary-like colors
    scale = "none"
  )
  dev.off()
  
  png(
    paste0(output_dir, "/adj_cor_mat/heatmap_adj_pink1_", day, "_.png"),
    width = png_width, 
    height = png_height, 
    res = png_res
  )
  heatmap(
    adj_matrix_pink1,
    na.rm = TRUE,
    Rowv = NA,
    Colv = NA,
    col = colorRampPalette(c("white", "black"))(3),
    scale = "none"
  )
  dev.off()
  #----------------------------------------------------------
  # Adjacency Matrix Network Visualization 
  #----------------------------------------------------------
  
  # Create the network graph
  network_graph_ctr <- graph_from_adjacency_matrix(adj_matrix_ctr, mode = "undirected")
  

  # Open the PNG device
  png(
    paste0(output_dir, "/adj_cor_mat/network_adj_ctr_basic_", day, "_.png"),
    width = 2400, 
    height = 1800, 
    res = png_res
  )
  
  # Set up the layout for two plots in one row
  layout(matrix(c(1, 2), nrow = 1, ncol = 2, byrow = TRUE))
  
  # First plot
  plot(network_graph_ctr, 
       vertex.size = 5,        # Size of nodes
       vertex.label.cex = 0.5,  # Font size for labels
       vertex.label.color = "black",  # Label color
       edge.width = 2,          # Width of edges
       edge.color = "gray")     # Color of edges
  
  # Modify node colors and edge width
  V(network_graph_ctr)$color <- c(colorRampPalette(c("black", "grey", "orange", "red"))(dim(adj_matrix_ctr)[1]))
  E(network_graph_ctr)$width <- 1
  
  # Second plot
  plot(network_graph_ctr,
       vertex.size = 5,        # Size of nodes
       vertex.label.cex = .6,  # Font size for labels
       vertex.label.color = "black",  # Label color
       edge.width = 4,          # Width of edges
       edge.color = "gray"
  )
  
  # Close the PNG device
  dev.off()
  
  # pink
  # Create the network graph
  network_graph_pink <- graph_from_adjacency_matrix(adj_matrix_pink1, mode = "undirected")
  
  # Open the PNG device first
  png(
    paste0(output_dir, "/adj_cor_mat/network_adj_ctr_basic_", day, "_.png"),
    width = 2400, 
    height = 1800, 
    res = png_res
  )
  
  # Set up the layout for two plots in one row after opening the device
  layout(matrix(c(1, 2), nrow = 1, ncol = 2, byrow = TRUE))
  
  # First plot with default settings
  plot(network_graph_pink, 
       vertex.size = 5,
       vertex.label.cex = 0.5,
       vertex.label.color = "black",
       edge.width = 2,
       edge.color = "gray")
  
  # Customize vertex colors based on the number of nodes in the graph
  num_nodes <- vcount(network_graph_pink)
  V(network_graph_pink)$color <- colorRampPalette(c("black", "grey", "orange", "red"))(num_nodes)
  E(network_graph_pink)$width <- 1
  
  # Second plot with updated attributes
  plot(network_graph_pink,
       vertex.size = 5,
       vertex.label.cex = 0.6,
       vertex.label.color = "black",
       edge.width = 4,
       edge.color = "gray")
  
  # Close the PNG device
  dev.off()
  
  
  ########################################################################################################  
  ##################### Highly Connected Nodes's Adjacency Matrix ########################################
  ######################################################################################################## 
  # Specify the folder path
  folder_path <- file.path(output_dir, "high_cor_mat")
  
  # Create the folder
  dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
  
  #############################Nodes #####################################################################
  nodes_ctr <- ctl$nodes
  nodes_pink1 <- pink1$nodes
  
  #----------------------------------------------------------
  # Histogram for CTL Nodes
  #----------------------------------------------------------
  png(
    paste0(output_dir, "/high_cor_mat/hist_nodes_ctr_", day, "_.png"),
    width = 2400,
    height = 1800,
    res = png_res
  )
  hist(
    nodes_ctr,
    breaks = 100,
    main = paste0("Node Distribution - CTL at day", day),  # Added title for clarity
    xlab = "Node Values",
    ylab = "Frequency"
  )
  dev.off()
  
  
  png(
    paste0(output_dir, "/high_cor_mat/hist_nodes_pink1_", day, "_.png"),
    width = png_width,
    height = png_height,
    res = png_res
  )
  hist(
    nodes_pink1,
    breaks = 100,
    main = paste0("Node Distribution - PINK1 at day", day),
    xlab = "Node Values",
    ylab = "Frequency"
  )
  dev.off()
  
  #############################high related Adjacency Matrix####################################################
  
  red_adj_matrix_ctl <- ctl$red_adj_matrix
  red_adj_matrix_pink1 <- pink1$red_adj_matrix
  
  #----------------------------------------------------------
  # Highly Connected Adjacency Matrix Heatmap
  #----------------------------------------------------------
  if (nrow(red_adj_matrix_ctl)<2   && ncol(red_adj_matrix_ctl) < 2){
    message("The matrix does not meet the requirements. Exiting the function...")
    return(NULL)
  }
  
  png(
    paste0(output_dir, "/high_cor_mat/heatmap_high_cor_ctl_", day, "_.png"),
    width = png_width,
    height = png_height,
    res = png_res
  )
  heatmap(
    red_adj_matrix_ctl,
    na.rm = TRUE,
    Rowv = NA,  # Disable row clustering
    Colv = NA,  # Disable column clustering
    col = colorRampPalette(c("white", "black"))(3),  # Binary-like colors
    scale = "none",
    main = paste0("Highly Connected Adjacency Matrix - CTL at day ", day)  # Add title
  )
  dev.off()
  
  png(
    paste0(output_dir, "/high_cor_mat/heatmap_high_cor_pink1_", day, "_.png"),
    width = png_width,
    height = png_height,
    res = png_res
  )
  heatmap(
    red_adj_matrix_pink1,
    na.rm = TRUE,
    Rowv = NA,
    Colv = NA,
    col = colorRampPalette(c("white", "black"))(3),
    scale = "none",
    main = paste0("Highly Connected Adjacency Matrix - PINK1 at day ", day)
  )
  dev.off()
  
  ###############network representation of the highly connected nodes######################################
  # Define improved parameters
  png_width <- 4000  # Larger canvas for better spacing (13.3 inches at 300 dpi)
  png_height <- 4000
  png_res <- 300
  vertex_size <- 15  # Increased node size
  label_cex <- 1.3  # Larger labels
  layout_zoom <- 1.5       # Spread out nodes
  
  #----------------------------------------------------------
  # CTL Network Plot
  #----------------------------------------------------------
  # Create output directory if needed
  dir.create(paste0(output_dir, "/high_cor_mat/"), showWarnings = FALSE, recursive = TRUE)
  
  # Create graph object
  red_network_graph_ctl <- graph_from_adjacency_matrix(
    red_adj_matrix_ctl, 
    mode = "undirected"
  )
  
  # Save plot
  png(
    paste0(output_dir, "/high_cor_mat/network_high_cor_ctl_", day, "_.png"),
    width = png_width, 
    height = png_height, 
    res = png_res
  )
  
  # Set layout and margins
  par(mar = c(0,0,1,0))  # Remove margins, keep space for title
  set.seed(123)  # For reproducible layout
  
  plot(red_network_graph_ctl,
       layout = layout_with_kk(red_network_graph_ctl, dim = 2) * layout_zoom,
       vertex.size = vertex_size,
       vertex.color = "#2c3e50",  # Dark blue nodes
       vertex.frame.color = NA,
       vertex.label.cex = label_cex,
       vertex.label.color = "red",
       vertex.label.family = "sans",  # Clean font
       vertex.label.dist = 2.,     # Distance from node center (in plot units)
       vertex.label.degree = pi/2,
       edge.width = 1.5,  # Thinner edges
       edge.color = adjustcolor("#7f8c8d", alpha.f = 0.3),  # Semi-transparent grey
       main = paste0("High Connectivity Network - CTL at day ", day),
       main.color = "#2c3e50",
       main.cex = 3,
       main.font = 2
  )
  
  dev.off()
  
  #----------------------------------------------------------
  # PINK1 Network Plot
  #----------------------------------------------------------
  red_network_graph_pink1 <- graph_from_adjacency_matrix(
    red_adj_matrix_pink1, 
    mode = "undirected"
  )
  
  png(
    paste0(output_dir, "/high_cor_mat/network_high_cor_pink1_", day, "_.png"),
    width = png_width, 
    height = png_height, 
    res = png_res
  )
  
  par(mar = c(0,0,2,0))
  set.seed(123)  # Same seed for comparable layouts
  
  plot(red_network_graph_pink1,
       layout = layout_with_kk,  # Kamada-Kawai layout (good for small networks)
       vertex.size = 8,
       vertex.color = "#e74c3c",  # Red nodes for contrast
       vertex.frame.color = NA,
       vertex.label.cex = 0.5,
       vertex.label.color = "#2c3e50",
       vertex.label.family = "sans",
       vertex.label.dist = 0.5,
       vertex.label.degree = pi,
       edge.width = 1.5,
       edge.color = adjustcolor("#7f8c8d", alpha.f = 0.3),
       main = paste0("High Connectivity Network - PINK1 at day ", day),
       main.color = "#2c3e50",
       main.cex = 2
  )
  
  dev.off()
  ############################ Comparing Conditions by Adjacency Pattern####################################

  
  # Subset the correlation matrix
  rows_to_keep_ctl <- ctl$rows_to_keep
  rows_to_keep_pink1 <- pink1$rows_to_keep
  com_ind<- c(which(rows_to_keep_ctl),which(rows_to_keep_pink1)) 
  
  shared_red_cor_matrix_ctr <- cor_matrix_ctr[com_ind, com_ind]
  shared_red_cor_matrix_pink <- cor_matrix_pink1[com_ind, com_ind]
  
  # Define consistent PNG parameters
  png_width <- 3600  # 12 inches at 300 dpi (wider for side-by-side plots)
  png_height <- 2400 # 8 inches at 300 dpi
  png_res <- 300
  
  h1 <- Heatmap(shared_red_cor_matrix_ctr, border=T,
                name = "Correlation", cluster_rows = FALSE, cluster_columns = FALSE, 
                show_row_names = FALSE, show_column_names = FALSE,
                col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
  
  h2 <- Heatmap(shared_red_cor_matrix_pink,  border=T,
                name = "Correlation", cluster_rows = FALSE, cluster_columns = FALSE, 
                show_row_names = FALSE, show_column_names = FALSE,
                col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))
  
  
  # Combined correlation heatmaps
  png(paste0(output_dir, "/high_cor_mat/heatmap_shared_high_cor_", day, "_.png"), 
      width = png_width, height = png_height, res = png_res)
  draw(h1 + h2, 
       ht_gap = unit(1, "cm"),  # Space between heatmaps
       padding = unit(c(2, 2, 4, 2), "cm"))  # Top, right, bottom, left margins
  dev.off()
  
  ##########################common adjacency as above####################################################
  
  ctresh <- 0.8
  com_adj_matrix_ctr <- ifelse(abs(shared_red_cor_matrix_ctr) > ctresh, 1, 0) #how to break this down into pos and negative
  diag(com_adj_matrix_ctr) <- 0
  
  # CTL adjacency matrix
  png(paste0(output_dir, "/high_cor_mat/heatmap_com_adj_matrix_ctr_",day,"_.png"), 
      width = png_width, height = png_height, res = png_res)
  heatmap(com_adj_matrix_ctr, na.rm = TRUE, Rowv = NA, Colv = NA,
          col = colorRampPalette(c("white", "black"))(3),
          scale = "none",
          margins = c(1,1))  # Reduce axis margins
  dev.off()
  
  #pink1
  com_adj_matrix_pink <- ifelse(abs(shared_red_cor_matrix_pink) > ctresh, 1, 0) #how to break this down into pos and negative
  diag(com_adj_matrix_pink) <- 0
  
  # PINK1 adjacency matrix
  png(paste0(output_dir, "/high_cor_mat/heatmap_com_adj_matrix_pink1_",day,"_.png"), 
      width = png_width, height = png_height, res = png_res)
  heatmap(com_adj_matrix_pink, na.rm = TRUE, Rowv = NA, Colv = NA,
          col = colorRampPalette(c("white", "black"))(3),
          scale = "none",
          margins = c(1,1))
  dev.off()
  
  h3 <- Heatmap(com_adj_matrix_ctr,  border=T,
                name = "Adjacency CTR", cluster_rows = FALSE, cluster_columns = FALSE, 
                show_row_names = FALSE, show_column_names = FALSE,
                col = colorRamp2(c( 0, 1), c( "white", "black")))
  
  h4 <- Heatmap(com_adj_matrix_pink,  border=T,
                name = "Adjecency PINK1", cluster_rows = FALSE, cluster_columns = FALSE, 
                show_row_names = FALSE, show_column_names = FALSE,
                col = colorRamp2(c(0, 1), c( "white", "black")))
  
  # Combined adjacency heatmaps
  png(paste0(output_dir, "/high_cor_mat/heatmap_common_high_adj_matrix_",day,".png"), 
      width = png_width, height = png_height, res = png_res)
  draw(h3 + h4,
       ht_gap = unit(1, "cm"),
       padding = unit(c(2, 2, 4, 2), "cm"),
       auto_adjust = FALSE)  # Prevent automatic size reduction
  dev.off()
}


gc()
gene_co_exp_day(0)

gc()
gene_co_exp_day(18)


gc()
gene_co_exp_day(25)


gc()
gene_co_exp_day(37)


gc()
gene_co_exp_day(57)


# ================================= LIBRARIES ==================================

library(glmGamPoi)
library(Seurat)
library(stringr)
library(ggplot2)
library(openxlsx)
library(Matrix)
library(genefilter) 
library(tidyverse)
library(SAVER)
library(modeest)
library(sctransform)

#ipsc used, after 21days develop into neurons, check what is the difference between wt and pd

rm(list=ls())
set.seed(420)
options(future.globals.maxSize=10000*1024^2)

# ================================ PARAMETERS ==================================

# assumes that working directory is project1 folder

setwd('C:/Users/yanis/Desktop/Universite/Master/Semestre_3/Bioinformatics/Final Project')

raw_folder = "./Data/"
SampleInfo_path = "./Data/SampleInfo.xlsx"
output_folder = "./Processed/"

export_plots = FALSE
impute = TRUE
size.factor = 1
n_cores = 12
remove_ribosomal = TRUE

subset_condition = c("WT","ND")
subset_days = c(0,18,25,37,57)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                            ---- PREPROCESSING ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# ---- Load metadata ----

SampleInfo = read.xlsx(SampleInfo_path)

# ---- Load files ----

S_raw = list()
counter = 1

for(file in list.files(raw_folder,full.names = TRUE)){
  
  condition_name = SampleInfo$Condition[counter]
  date = SampleInfo$Date[counter]
  sample_number = counter
  sample = paste(sample_number,condition_name,date,sep="_")
  condition_day = paste(condition_name,date,sep="_")
  output_file = paste(output_folder,sample,sep="")
  counter = counter + 1
  
  if(condition_name %in% subset_condition & date %in% subset_days){
    
    print(condition_day)
    
    # Loading raw matrix
    
    M = as.matrix(read.csv(file, header = T, row.names = 1, sep = "\t"))
    
    # Renaming cells and genes
    
    colnames(M) =  paste(condition_day,as.character(1:ncol(M)),sep="_")  
    rownames(M) = str_replace_all(toupper(rownames(M)), "-",".")
    rownames(M) = str_replace_all(toupper(rownames(M)), "_",".")
    
    # Creating Seurat Object
    
    S_raw[condition_day] = CreateSeuratObject(counts = M, min.cells = 0, min.features = 0)
    
  }
  
}

Samples = names(S_raw)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                                ---- QC ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

S_list_filtered = list()
S_list = list()

nArray = data.frame(Name = character(),nCount = numeric(),nFeature = numeric())

counter = 1

QC_df = data.frame(Sample = Samples,mode_mt = numeric(10),mode_rb = numeric(10),mode_nFeat = numeric(10),mode_nCount = numeric(10))

for(SeuratObj in S_raw){
  
  print(Samples[counter])
  
  # Mitochondrial RNA
  SeuratObj[["percent.mt"]] = PercentageFeatureSet(SeuratObj, pattern = "^MT-|^MT\\.")
  
  # Ribosomal RNA
  SeuratObj[["percent.rb"]] = PercentageFeatureSet(SeuratObj, pattern = "^RPL|^RPS|^RPSA|^RPP")
  
  SeuratObj = subset(SeuratObj, subset = percent.mt > 0 & percent.rb > 0)
  
  # Number of features (standardized)
  mode_nFeat = as.numeric(mlv(SeuratObj$nFeature_RNA,method = "asselin", bw = 1,na.rm = TRUE)[1])
  SeuratObj[["nFeature_RNA_standardized"]] = scale(SeuratObj$nFeature_RNA - mode_nFeat, center = FALSE, scale = TRUE)
  
  # Number of counts (standardized)
  mode_nCount = as.numeric(mlv(SeuratObj$nCount_RNA,method = "asselin", bw = 1,na.rm = TRUE)[1])
  SeuratObj[["nCount_RNA_standardized"]] = scale(SeuratObj$nCount_RNA - mode_nCount, center = FALSE, scale = TRUE)
  
  # Mitochondrial RNA (standardized)
  mode_mt = as.numeric(mlv(SeuratObj$percent.mt,method = "asselin", bw = 1,na.rm = TRUE)[1])
  SeuratObj[["percent.mt_standardized"]] = scale(SeuratObj$percent.mt - mode_mt, center = FALSE, scale = TRUE)
  
  # Ribosomal RNA (standardized)
  mode_rb = as.numeric(mlv(SeuratObj$percent.rb,method = "asselin", bw = 1,na.rm = TRUE)[1])
  SeuratObj[["percent.rb_standardized"]] = scale(SeuratObj$percent.rb - mode_rb, center = FALSE, scale = TRUE)
  
  #Filter
  S_list[[Samples[counter]]] = SeuratObj
  S_list_filtered[[Samples[counter]]] = subset(SeuratObj, subset = percent.mt_standardized < 2 & nCount_RNA_standardized < 3)
  
  if(Samples[counter] == "PINK_0"){
    S_list_filtered[[Samples[counter]]] = subset(SeuratObj, subset = percent.mt_standardized < 2 & nCount_RNA_standardized < 3 & nCount_RNA > 1600)
  }
  
  # Fill dataframe
  QC_df$mode_mt[counter] = mode_mt
  QC_df$mode_rb[counter] = mode_rb
  QC_df$mode_nFeat[counter] = mode_nFeat
  QC_df$mode_nCount[counter] = mode_nCount
  
  nArray = rbind(nArray,data.frame(Name = Samples[counter],nCount = S_list[[Samples[counter]]]$nCount_RNA,nFeature = S_list[[Samples[counter]]]$nFeature_RNA))
  
  counter = counter + 1
  
}

# Removing ribosomal genes

if(remove_ribosomal){
  for(s in Samples){
    non_rib_genes <- grep(pattern = "^RPL|^RPS|^RPSA|^RPP",
                          rownames(S_list_filtered[[s]]),
                          value=TRUE, invert=TRUE)
    S_list_filtered[[s]] = subset(S_list_filtered[[s]],features = non_rib_genes)
  }
}


saveRDS(S_list_filtered,paste(output_folder,"SeuratFiltered.rds",sep=""))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                           ---- QC PLOTS ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# to safe plots u may have to generate the corresponding folder
output_folder = "./Preprocessed_Data/QC_Plots/"

if(export_plots){
  for(S in Samples){
    
    df_unfilt = data.frame(
      cell_id = names(S_list[[S]]$percent.mt),
      percent.mt = as.numeric(S_list[[S]]$percent.mt),
      percent.rb = as.numeric(S_list[[S]]$percent.rb),
      nFeature = as.numeric(S_list[[S]]$nFeature_RNA),
      nCount = as.numeric(S_list[[S]]$nCount_RNA),
      percent.mt_standardized = as.numeric(S_list[[S]]$percent.mt_standardized),
      percent.rb_standardized = as.numeric(S_list[[S]]$percent.rb_standardized),
      nFeature_standardized = as.numeric(S_list[[S]]$nFeature_RNA_standardized),
      nCount_standardized = as.numeric(S_list[[S]]$nCount_RNA_standardized))
    
    df_unfilt = gather(df_unfilt,"Type","Value",2:9)
    
    df_filt = data.frame(
      cell_id = names(S_list_filtered[[S]]$percent.mt),
      percent.mt = S_list_filtered[[S]]$percent.mt,
      percent.rb = S_list_filtered[[S]]$percent.rb,
      nFeature = S_list_filtered[[S]]$nFeature_RNA,
      nCount = S_list_filtered[[S]]$nCount_RNA,
      percent.mt_standardized = S_list_filtered[[S]]$percent.mt_standardized,
      percent.rb_standardized = S_list_filtered[[S]]$percent.rb_standardized,
      nFeature_standardized = S_list_filtered[[S]]$nFeature_RNA_standardized,
      nCount_standardized = S_list_filtered[[S]]$nCount_RNA_standardized)
    
    df_filt = gather(df_filt,"Type","Value",2:9)
    
    pdf(file = paste(output_folder,S,".pdf",sep=""), width=30, height=12)
    
    print(df_unfilt %>% ggplot(aes(x = "",y = Value)) + geom_violin(adjust = 0.8) +
            geom_jitter(alpha = 0.05,size = 0.6,fill = "black") +
            facet_wrap(~Type,scales = "free",ncol = 8) +
            labs(title = paste("Sample ",S," before QC filtering")))
    
    print(df_filt %>% ggplot(aes(x = "",y = Value)) + geom_violin(adjust = 0.8) +
            geom_jitter(alpha = 0.05,size = 0.6,fill = "black") +
            facet_wrap(~Type,scales = "free",ncol = 8) +
            labs(title = paste("Sample ",S," after QC filtering")))
    
    dev.off()
    
  }
}

output_folder = "./Preprocessed/"

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                           ---- ADD METADATA ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

for(s in Samples){
  
  S_list_filtered[[s]] = AddMetaData(object = S_list_filtered[[s]], 
                                     metadata = str_split(s,"_")[[1]][1], 
                                     col.name = "Condition")
  
  S_list_filtered[[s]] = AddMetaData(object = S_list_filtered[[s]], 
                                     metadata = as.numeric(str_split(s,"_")[[1]][2]), 
                                     col.name = "Day")
  
  S_list_filtered[[s]] = AddMetaData(object = S_list_filtered[[s]], 
                                     metadata = s, 
                                     col.name = "ConditionDay")
  
}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                   ---- VARIANCE STABILIZATION TRANSFORM ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

Condition_per_sample = unlist(str_split(Samples,"_"))[seq(1,20,2)]
Day_per_sample = unlist(str_split(Samples,"_"))[seq(2,20,2)]

# Merge all samples in one Seurat object (S)

S = list()
is.first = TRUE

#S = merge(S_list_filtered[[Samples[1]]], c(S_list_filtered[[Samples[2]]],S_list_filtered[[Samples[3]]],S_list_filtered[[Samples[4]]],S_list_filtered[[Samples[5]]],S_list_filtered[[Samples[6]]],S_list_filtered[[Samples[7]]],S_list_filtered[[Samples[8]]],S_list_filtered[[Samples[9]]],S_list_filtered[[Samples[10]]]))

for(s in Samples){
  print(s)
  if(is.first){
    S = S_list_filtered[[s]]
    is.first = FALSE
  }
  else{
    S = merge(S,S_list_filtered[[s]])
  }
}

S <- JoinLayers(S) # to be done for Seurate 5
# Run VST with the batch as dependent variable

S$orig.ident = S$ConditionDay
S$log_umi_per_gene = log10(S$nCount_RNA/S$nFeature_RNA)
S$batch = S$ConditionDay

# Attention - this will take quite some time (1h or so)
S_vst = vst(S@assays$RNA$counts, cell_attr = S@meta.data, 
            latent_var = "log_umi_per_gene", 
            batch_var = "batch", 
            method = "glmGamPoi",
            return_corrected_umi = TRUE)

# Create new Seurat object with corrected counts and pearson residuals

S_corrected = CreateSeuratObject(counts = S_vst$umi_corrected)
#S_corrected@assays$RNA@layers$scale.data = S_vst$y
# or better ? by 
#S_corrected[["RNA"]]@scale.data
S_corrected$Condition = S$Condition
S_corrected$Day = S$Day
S_corrected$ConditionDay = S$ConditionDay

saveRDS(S_corrected,paste(output_folder,"SeuratCorrected.rds",sep=""))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                           ---- DATA IMPUTATION ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 


# NOTE: will take long time - we go through the steps but will use final object later
S_list_corrected = list()

S_corrected = subset(S_corrected,cells = which(S_corrected$nCount_RNA!=0))

cells_to_keep = colnames(S_corrected)

for(s in Samples){
  print(s)
  S_list_corrected[[s]] = subset(S_corrected,subset = ConditionDay == s)
}

S_list_imputed = list()

counter = 1

for(SeuratObj in S_list_corrected){
  
  print(Samples[counter])
  
  if(impute){
    Saver_object = saver(SeuratObj@assays$RNA$counts, ncores = n_cores, size.factor = 1)
    SeuratObj[["RNA_imputed"]] = CreateAssayObject(counts = Saver_object$estimate)
    S_list_imputed[[Samples[counter]]] = SeuratObj
    counter = counter + 1
  }
  
}

saveRDS(S_list_imputed,paste(output_folder,"SeuratImputed.rds",sep=""))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                         ---- DATA INTEGRATION ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Removing ribosomal genes from corrected assay 

non_rib_genes = grep(pattern = "^RPL|^RPS|^RPSA|^RPP",
                     rownames(S_corrected),
                     value=TRUE, invert=TRUE)

S_corrected_nonrib = subset(S_corrected,features = non_rib_genes)

S_merged_Cond = list()

S_merged_Cond[["WT"]] = subset(S_corrected_nonrib, subset = Condition=="WT")
S_merged_Cond[["ND"]] = subset(S_corrected_nonrib, subset = Condition=="ND")

S_merged_Cond[["WT"]] = NormalizeData(S_merged_Cond[["WT"]])
S_merged_Cond[["ND"]] = NormalizeData(S_merged_Cond[["ND"]])

features = SelectIntegrationFeatures(object.list = S_merged_Cond)
anchors = FindIntegrationAnchors(object.list = S_merged_Cond, anchor.features = features, scale = FALSE, normalization.method = "LogNormalize")                            
S_integrated = IntegrateData(anchorset = anchors)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                     ---- DIMENSIONALITY REDUCTION ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

S_integrated = ScaleData(S_integrated)
S_integrated = RunPCA(S_integrated,npcs = 100)
S_integrated = RunUMAP(S_integrated,dims = 1:100,assay = "integrated",n.neighbors = 30)

DimPlot(object = S_integrated, reduction = "umap", group.by = "Condition", split.by = "Day",shuffle = TRUE) 

saveRDS(S_integrated,paste(output_folder,"SeuratIntegrated.rds",sep=""))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                         ---- UNIFY DATASETS ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Merge all samples in one Seurat object (S)

is.first = TRUE

for(s in Samples){
  print(s)
  if(is.first){
    S_imputed_merged = S_imputed[[s]]
    is.first = FALSE
  }
  else{
    S_imputed_merged = merge(S_imputed_merged,S_imputed[[s]])
  }
}
f
S_integrated[["RNA_imputed"]] = CreateAssayObject(counts = S_imputed_merged@assays$RNA_imputed@counts)

saveRDS(S_integrated,paste(output_folder,"SeuratFinal.rds",sep=""))

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#                         ---- WORK ON PREPROCESSED DATA ----                  
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

preprocessed_data <- readRDS("./Preprocessed/SeuratFinal.rds")

summary(preprocessed_data)

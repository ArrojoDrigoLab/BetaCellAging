##
##      __________  _____
##     / ____/ __ \/ ___/
##    / /   / / / /\__ \ 
##   / /___/ /_/ /___/ / 
##   \____/_____//____/  
##
##  Creative Data Solutions
##  Vanderbilt University
##  https://cds.vanderbilt.edu
##
## Authors: Drs. Shristi Shrestha & JP Cartailler
## Date Created: 2022-08-10
## 

#############################################

# scRNA-seq meta-analysis of following published studies
# 1. Segerstolpe et al.2016
# 2. Xin et al.2018
# 3. Human Pancreas Analysis Program(HPAP)
# 4. Camunas et al. 2020
# 5. Enge et al.2017
# 6. Xin et al.2016
# 7. Tosti et al. 2020


############################################

#load libraries

library(Seurat) #v3.2.3
library(tidyverse)
library(ggplot2)


# load data---------------------------------------------

# seurat object for Tosti et al 2020,Xin et al.2018, Camunas et al. 2020 were created with
# pre-processed expression count table availabe in respective data deposits

# Xin et al.2016, Enge et al.2017 and Segerstolpe et al. 2016 were processed from raw reads provided in respective data deposits.
# See Methods for pre-processing of these 3 datasets 

# Final merged/integrated and processed Rds provided at the end of this script

# Tosti et al. 2020 dataset
gastro <- readRDS("data/Seurat_objects/Seurat_Gastro_January2021.rds")

# Xin et al.2018 dataset
gromada <- readRDS("data/Seurat_objects/Seurat_gromada_January2021.rds")

# Human Pancreas Analysis Program(HPAP)
hpap <- readRDS("data/Seurat_objects/Seurat_hpap_January2021.rds")

# Camunas et al. 2020 dataset
joan <- readRDS("data/Seurat_objects/Seurat_Joan_January2021.rds")

# Combined Xin et al.2016, Enge et al.2017, Segerstolpe et al. 2016
meta <- readRDS("data/Seurat_objects/Seurat_meta_January2021.rds")

#dimension of all datasets
dim(gastro) #2000 genes 13182 cells
dim(gromada) # 19167 genes 20623 cells
dim(hpap) #21812 genes  4456 cells
dim(joan)#48898 genes  5359 cells
dim(meta)#11362 genes  6559 cells

# Tosti et al. 2020 dataset(young donors) left out of all preprocess filtering
objects <- list(gromada,hpap,joan,meta)

# 1. QC, pre-processing  --------------------------------------------------


for (object in objects) {
    
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
    object[["percent.ribo"]] <- PercentageFeatureSet(object, pattern = "^RPS|^RPL")
    
    #Exploratory data analysis(gene,UMI, mito count distribution) was visualized
}

# 2. Based on EDA , filtered cells with genes <200 and transcript <3000----
# gastro(young donors) left out of all preprocess filtering

filtered_gromada <- subset(gromada, subset = nFeature_RNA >200 & nCount_RNA >3000) #32 cells removed
filtered_hpap <- subset(hpap, subset = nFeature_RNA >200 & nCount_RNA >3000) #24 cells removed
filtered_joan <- subset(joan, subset = nFeature_RNA >200 & nCount_RNA >3000) #144 cells removed
filtered_meta <- subset(meta, subset = nFeature_RNA >200 & nCount_RNA >3000) #178 cells removed

# 3. Normalize ---------------------------------------------

filtered_gromada <- NormalizeData(filtered_gromada, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_hpap <- NormalizeData(filtered_hpap, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_joan <- NormalizeData(filtered_joan, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_meta <- NormalizeData(filtered_meta, normalization.method = "LogNormalize", scale.factor = 10000)

# 4. Supervised filtering to remove double positive hormonal cells in each dataset ----

# gromada:
FeatureScatter(filtered_gromada, feature1 = "INS",feature2 = "GCG", group.by = "agegroup")
filtered_gromada <- subset(x = filtered_gromada, subset = INS >  6 & GCG >5.5, invert=TRUE)
dim(filtered_gromada) #19167 20157 , 434 cells removed

# hpap:
FeatureScatter(filtered_hpap, feature1 = "INS",feature2 = "GCG", group.by = "agegroup")
filtered_hpap <- subset(x = filtered_hpap, subset = INS >  4.5 & GCG >3.5, invert=TRUE)
dim(filtered_hpap)#21812  4191, 241 cells removed

# joan:
FeatureScatter(filtered_joan, feature1 = "INS",feature2 = "GCG", group.by = "agegroup")
filtered_joan <- subset(x = filtered_joan, subset = INS >  5.5 & GCG >4.5, invert=TRUE)
dim(filtered_joan) #48898  5095, 120 cells removed

# meta:
FeatureScatter(filtered_meta, feature1 = "INS",feature2 = "GCG", group.by = "agegroup")
filtered_meta <- subset(x = filtered_meta, subset = INS >  4 & GCG >4.5, invert=TRUE)
dim(filtered_meta)#11362  6329 #52 cells removed

# total 847 cells removed

# 5. Combining datasets (gastro, gromada, hpap, joan,meta)----
DefaultAssay(object = gastro) <- "RNA"
panc.combined <- merge(gastro, y = c(filtered_gromada, filtered_hpap,filtered_joan,filtered_meta), add.cell.ids = c("gastro", "gromada", "hpap", "joan", "meta"), project = "AllPanc")
panc.combined #66628 48954

# remove(gastro, gromada,hpap,joan,meta, filtered_gromada, filtered_hpap,filtered_joan,filtered_meta)
# saveRDS(panc.combined, "unintegrated_data.rds")
# panc.combined<-readRDS("unintegrated_data.rds")

# 6. Finding anchors and integrate datasets-----
pancreas.list <- SplitObject(panc.combined, split.by = "tech")
pancreas.list <- pancreas.list[c("meta", "gromada", "hpap", "facs", "cryo", "gastro")]

for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
    pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
                                               nfeatures = 8000, verbose = FALSE)
}

# produce 3000 anchor genes which will be used to build anchors between datasets.
features <- SelectIntegrationFeatures( object.list = pancreas.list, nfeatures = 3000 , verbose = TRUE )

reference.list <- pancreas.list[c("meta", "gromada", "hpap", "facs", "cryo", "gastro")]

pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list,
                                           anchor.features = features,
                                           dims = 1:30)

all_features <- lapply(pancreas.list, row.names) %>% Reduce(intersect, .) 
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, 
                                     dims = 1:30, 
                                     features.to.integrate= all_features)

nrow(pancreas.integrated[["integrated"]]@scale.data)
length(all_features)




# 7. Peek at highly variable genes-------------------
pancreas.integrated <- FindVariableFeatures(pancreas.integrated, 
                                            selection.method = "vst", 
                                            nfeatures = 2000)
top10 <- head(VariableFeatures(pancreas.integrated), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pancreas.integrated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)


# 8. Scaling------------------------------------------------------ 
all.genes <- rownames(pancreas.integrated)
pancreas.integrated <- ScaleData(pancreas.integrated, features = all.genes)

saveRDS(pancreas.integrated, "data/cds/pancreas.integrated.Rds")
#pancreas.integrated<-readRDS("data/cds/pancreas.integrated.Rds")


# 9. Perform PCA--------------------------------------------------
#All cells
pancreas.integrated <- RunPCA(pancreas.integrated, features = VariableFeatures(object = pancreas.integrated))

# PCA Heatmaps
DimHeatmap(pancreas.integrated, dims = 1:5, cells = 500, balanced = TRUE)

# Determining PCs to use downstream
ElbowPlot(pancreas.integrated)


# 10. Clustering cells-------------------------------------------
pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:20)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.7)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:20)

# Clustering visuals
DimPlot(pancreas.integrated, reduction = "umap", label = TRUE)
DimPlot(pancreas.integrated, reduction = "umap", label = FALSE)

# Stash current cluster ident
pancreas.integrated[["cluster_res0.7"]] <- Idents(object = pancreas.integrated)
Idents(object = pancreas.integrated) <- "cluster_res0.7"

dim(pancreas.integrated) #66628 genes 48954 cells

# 11. Fix annotation in metadata--------------------------------------
metadata<-pancreas.integrated@meta.data
table(metadata$sex)
table(metadata$agegroup)
table(metadata$study)

# combine <6 and 0-6 age group
metadata$agegroup[metadata$agegroup == "<6"] <- "0-6"

# capitalize the first letter of sex
library(tidyverse)
str_sub(metadata$sex, 1, 1) <- str_sub(metadata$sex, 1, 1) %>% str_to_upper()
table(metadata$sex)

# added name gastro in study metadata column
metadata$study[is.na(metadata$study)] <- "gastro"

# replace donor info for gastro dataset with patient ID since donor metadata was empty for gastro----
metadata$donor <- ifelse(is.na(metadata$donor), metadata$sample_ID, metadata$donor)

# add the metadata back to seurat obj
pancreas.integrated<-AddMetaData(pancreas.integrated, metadata = metadata)



# 12. CellType Annotation-------------------------------------------

# Plotting known cell type marker genes
Pancreas.markers <- c("GCG","IRX2","PCSK2","ARX",
                    "INS","IAPP","PDX1","MAFA", "MAFB","PDX1",
                    "SST","HHEX","LEPR","PPY","GHRL",
                    "PHGR1", "PRSS1","REG1A","CPA2", 
                    "KRT19","PROM1","CFTR",
                    "PECAM1","PLVAP","VWF", "PDGFRB","TIMP1","COL1A1", 
                    "HLA-DRA","CD68","CD74",
                    "VWF", "VIM", "CD3", "COL1A1", "MUC6", "CHGA", "FEV", "GFAP")

DotPlot(pancreas.integrated, features =Pancreas.markers)+ 
    scale_color_gradient2(low = "blue", high = "red",mid = "white")+
    labs(y=" ", x="")+RotatedAxis()+
    theme(text = element_text(size = 18),
          axis.text = element_text(size = 18))


# 13. Assigning cell type identity to clusters when resolution=1.6----

current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,13,14,15,16,17,18,19,20,21,22)
new.cluster.ids <- c("Beta","Alpha","Ductal","Alpha","Acinar","Delta-PP-Ep",6,"Ductal","Stellate","Alpha","Endothelial","Beta","Beta",13,"Stellate","Immune","Ductal",17,18,"Delta-PP-Ep","Stellate",21,"Ductal")                          
pancreas.integrated@active.ident <- plyr::mapvalues(x = pancreas.integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)

DimPlot(object = pancreas.integrated, pt.size = 0.8, label = TRUE)+
    theme(text = element_text(size = 18,face="bold"),
          axis.text = element_text(size = 18,face="bold"))

# stash cell type identities for later

pancreas.integrated[["CellTypes"]] <- Idents(object = pancreas.integrated)
Idents(object = pancreas.integrated) <- "CellTypes"



# 14. Cell type cells counts----------------------------------------
table(pancreas.integrated@meta.data$agegroup,pancreas.integrated@meta.data$CellTypes)
table(pancreas.integrated@meta.data$agegroup,pancreas.integrated@meta.data$tech)
prop.table(pancreas.integrated@meta.data$CellTypes)

# 15. Retain only protein coding genes-----------------------------

# Connect to AnnotationHub
ah <- AnnotationHub()
# Access the Ensembl database for human
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)
# Acquire the latest annotation files
id <- ahDb %>%
    mcols() %>%
    rownames() %>%
    tail(n = 1)
# Download the appropriate Ensembldb database
edb <- ah[[id]]
# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
    dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

protein_coding<-annotations %>% filter(gene_biotype=="protein_coding" )
panc.filtered <- pancreas.integrated[protein_coding$gene_name, ]
dim(panc.filtered) #19407 48954

# 17. Save Rds panc.filtered with non-protein coding genes removed ---------
saveRDS(panc.filtered, "data/panc.filtered.Rds")
dim(panc.filtered) #19407 genes 48954 cells
# This Rds available in zenodo deposit, see link in Readme



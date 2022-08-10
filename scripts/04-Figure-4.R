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


# Libraries
library(dplyr)
library(Seurat) #3.2.3
library(ggplot2)
library(stringr)
library(ggpubr)
library(viridis)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(readr)
library(dplyr)
library(ggcorrplot)
library(Matrix)
library(plotly)

# Figure 4B ----
# Heatmap of pyscenic regulons among all cell types


# Prepare input files for pySCENIC analysis
# Load pre-processed and cell type annotated data from scripts/01-Preprocessing.R

panc.filtered <- readRDS("data/cds/panc.filtered.Rds")

Idents(panc.filtered) <- "status"
levels(panc.filtered)

# Subset cells from only healthy donors
allcelltype_healthy <- subset(panc.filtered, idents = "healthy")
dim(allcelltype_healthy) #19407 38754

Idents(allcelltype_healthy) <- "CellTypes"
allcelltype_healthy <- subset(allcelltype_healthy, idents = c("Beta","Alpha","Delta-PP-Ep","Acinar","Ductal","Stellate","Endothelial","Immune"))
dim(allcelltype_healthy) #19407 35566

# Randomize and subset 13,000 cells (arbitrary limit to run pySCENIC analysis for this figure only)
set.seed(111)
allcelltype_healthy_13000 <- allcelltype_healthy[, sample(colnames(allcelltype_healthy), size =13000, replace=F)]
dim(allcelltype_healthy_13000)#19407 13000

# This Rds was used for pySCENIC analysis for figure 4B
saveRDS(allcelltype_healthy_13000, "allcelltype_healthy_13000.Rds")

# After pySCENIC analysis, pyscenic binary matrix and auc scores were added back to this same Rds (allcelltype_healthy_13000.Rds)
# auc scores are in assay "AUC" and binary matrix are added to "BIN" assay

obj <- readRDS("allcelltype_healthy_13000.pySCENIC.rds") #available in zenodo deposit, see link in Readme

#Prepare metadata to add to heatmap
agegroup_anno <- obj@meta.data$agegroup
sex_anno <- obj@meta.data$sex
study_anno <- obj@meta.data$study
celltype_anno <- obj@meta.data$CellTypes

# Add % ribo metadata
DefaultAssay(obj) <- "RNA"
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = '^RPL|^RPS' )

# Assay containing the pyscenic binary matrix, auc score are in assay "AUC"
DefaultAssay(obj) <- "BIN"

# Get matrix to be plotted in heatmap
mat <- obj[["BIN"]]@data %>% as.matrix()

# Choose heatmap colors
col_fun = circlize::colorRamp2(c(0, 1), c("white","black")) #color for heatmap
col_fun2 = circlize::colorRamp2(c(0, 15, 30), c("blue", "white", "red")) #color for percent.ribo

# Choose regulons to show and match row index
regulons_to_show_in_rows <- c("PDX1","MAFA","MAFB","NKX6-1","NKX2-2","SIX2","SIX3","PAX6","PTF1A","PAX4","ARX","IRX2","ISL1","NEUROD1","SOX9","ONECUT2","GATA6","RBPJL","NANOG","IRF2","TBX2","XBP1","ATF4")
row_index <- match(regulons_to_show_in_rows,rownames(mat))

# Heatmap plot
set.seed(100)
heat<-Heatmap(mat, name = "Expression",  
              column_split = factor(celltype_anno, levels = c("Beta","Alpha","Delta-PP-Ep","Acinar","Ductal","Stellate","Endothelial","Immune")),
              cluster_columns = TRUE,
              show_column_dend = TRUE,
              cluster_column_slices = TRUE,
              column_title_gp = gpar(fontsize = 8),
              column_gap = unit(0.5, "mm"),
              cluster_rows = TRUE,
              row_km = 10,
              show_row_dend = TRUE,
              col = col_fun,
              row_names_gp = gpar(fontsize = 4),
              column_title_rot = 90,
              show_row_names = FALSE,
              top_annotation = HeatmapAnnotation(celltype = celltype_anno,
                                                 percent.ribo = obj$percent.ribo,
                                                 agegroup = agegroup_anno,
                                                 sex =sex_anno,
                                                 study = study_anno,
                                                 col = list(celltype = c('Alpha'='#F8766D','Beta'='#39B600','Delta-PP-Ep'='#D89000','Acinar'='#00BFC4','Ductal'='#0035f6','Endothelial'='#f64600','Stellate'='#E76BF3','Immune'='#FF62BC'),
                                                            percent.ribo = col_fun2,
                                                            agegroup =  c("0-6"="white" , "14-18" = "yellow" , "20-49" = "#54A8DD", "50-76" =  "#D83089"),
                                                            sex = c("Male" = "red", "Female" = "blue"),
                                                            study = c("Segerstolpe" = "#FF62BC","HPAP" = "#00B0F6","gromada" = "#00BFC4", "joan" = "#F8766D","Quake" = "#A3A500","gastro" = "#39B600", "Xin"="#FFFF00"))),
              right_annotation = rowAnnotation(foo = anno_mark(at = row_index, labels = regulons_to_show_in_rows)),
              show_column_names = FALSE,
              raster_quality = 4
)


heat=draw(heat)


#Figure 4C ----

# Extract kmeans cluster 8 regulon identities to plot in beta cells only
row_order_list <- row_order(heat) 

# convert row_order_list to dataframe
kmean_clusters_df <- enframe(row_order_list) %>%
    unnest(cols = c(value)) %>% 
    mutate(regulons= rownames(mat[value,]))%>%
    mutate(kmean = paste0("cluster",name))

# Subset only beta cells from from the object

Idents(obj)
Beta_BIN<-subset(obj, idents = ("Beta"))

# Add ribo metadata
DefaultAssay(Beta_BIN) <- "RNA"
Beta_BIN[["percent.ribo"]] <- PercentageFeatureSet(Beta_BIN, pattern = '^RPL|^RPS' )

# Get matrix to be plotted in heatmap for beta cells only
DefaultAssay(Beta_BIN) <- "BIN"
mat <- Beta_BIN[["BIN"]]@data %>% as.matrix()

#regulons to display
cluster_subset <- kmean_clusters_df %>% filter(kmean == "cluster8") 
regulons_to_show <- cluster_subset$regulons
row_index <- match(regulons_to_show, rownames(mat))

#subset only the cluster8 regulons
mat2 <- mat[row_index,] 

# Heatmap
set.seed(100)
known_heat <- Heatmap(mat2, name = "Expression",  
                    cluster_columns = TRUE,
                    show_column_dend = TRUE,
                    cluster_column_slices = FALSE,
                    column_title_gp = gpar(fontsize =12),
                    cluster_rows = FALSE,
                    show_row_dend = FALSE,
                    col = col_fun,
                    row_names_gp = gpar(fontsize = 14),
                    top_annotation = HeatmapAnnotation(celltype = Beta_BIN$CellTypes,
                                                       percent.ribo = Beta_BIN$percent.ribo,
                                                       agegroup = Beta_BIN$agegroup,
                                                       sex =Beta_BIN$sex,
                                                       study = Beta_BIN$study,
                                                       col = list(celltype = c('Alpha'='#F8766D','Beta'='#39B600','Delta-PP-Ep'='#D89000','Acinar'='#00BFC4','Ductal'='#0035f6','Endothelial'='#f64600','Stellate'='#E76BF3','Immune'='#FF62BC'),
                                                                  percent.ribo = col_fun2,
                                                                  agegroup =  c("0-6"="white" , "14-18" = "yellow" , "20-49" = "#54A8DD", "50-76" =  "#D83089"),
                                                                  sex = c("Male" = "red", "Female" = "blue"),
                                                                  study = c("Segerstolpe" = "#FF62BC","HPAP" = "#00B0F6","gromada" = "#00BFC4", "joan" = "#F8766D","Quake" = "#A3A500","gastro" = "#39B600", "Xin"="#FFFF00"))),
                    show_column_names = FALSE,
                    show_row_names = TRUE,
                    raster_quality = 4
)
known_heat=draw(known_heat)


# Figure 4F ----
# Read data
# this worksheet contains top 15 % regulons ranked by their pyscenic "importance" metric. 
# The matrix are binary counts output from pyscenic analysis
# Cells are only included from healthy beta cells

datafile <- 'data/top15pct_regulon_matrix_beta_healthy2.csv'
data <- read.table(file = datafile, header = TRUE, row.names = 1, sep = ",")
data[1:3, 1:3]


TF_all <- rownames(data)
TF_all_unique <- unique(TF_all)
as.data.frame(TF_all) %>% filter(TF_all == '')

# If rows are TFs, then tranpose the data 
data_t <- t(data)
data_t[1:3, 1:3]

# Correlation analysis
# method = c("pearson", "kendall", "spearman"))
data_t.cor <- data_t %>% cor(use="pairwise.complete.obs") # ~1min, 
dim(data_t.cor) #264 264

# Compute a correlation matrix p-values.
pvalue.mat <- cor_pmat(data_t.cor)
head(pvalue.mat[, 1:4])

# Heatmap
g <- ggcorrplot(data_t.cor, p.mat = pvalue.mat, hc.order = TRUE, type = "full", insig = "blank", sig.level = 1e-100)
g
ggsave("corr_plot.png", g, scale=5, height = 8, width = 8.5, dpi=300)

ggplotly(g) # interactive plotly

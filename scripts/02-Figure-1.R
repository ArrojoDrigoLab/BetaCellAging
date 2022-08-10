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

# Load libraries
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)

# Load pre-processed and cell type annotated data from scripts/01-Preprocessing.R----
panc.filtered <- readRDS("data/cds/panc.filtered.Rds")

# Figure1B ----
# UMAP clustering of pooled cells from seven single cell RNAseq human pancreatic islet 
# Donors of all age and diseases are included.
# Only beta cell cluster(from healthy donors) was considered for further analysis beyond this figure

# Remove donor R257(duplicate donor id for T2D and elevated HBA1c) for this umap.
Idents(panc.filtered)<-"donor"
panc.filtered<-subset(panc.filtered, idents = "R257", invert=TRUE )
#check if removal was ok
check<-as.data.frame(table(panc.filtered$donor))
"R257" %in% check$Var1


DimPlot(object = panc.filtered, pt.size = 0.8, group.by = "CellTypes",
        cols = c('Alpha'='#F8766D','Beta'='#39B600','Delta-PP-Ep'='#D89000','Acinar'='grey','Ductal'='grey','Endothelial'='grey','Stellate'='grey','Immune'='grey',"13"="grey","6"="grey","17"="grey","18"="grey","21"="grey"))+
    theme(text = element_text(size = 14),
          axis.text = element_text(size = 14))



# Supplemental Figure 1A----
p1<-DimPlot(object = panc.filtered, pt.size = 0.8,group.by = "cluster_res0.7")+
    theme(text = element_text(size = 14,face="bold"),
          axis.text = element_text(size = 14,face="bold"))+
    labs(title = "Clusters")

p2<-DimPlot(panc.filtered, reduction = "umap", group.by = "study")+
    theme(text = element_text(size = 14,face="bold"),
          axis.text = element_text(size = 14,face="bold"))+
    labs(title = "Data Source")

p3<-DimPlot(panc.filtered, reduction = "umap", group.by = "status")+
    theme(text = element_text(size = 14,face="bold"),
          axis.text = element_text(size = 14,face="bold"))+
    labs(title = "Disease status")

p4<-DimPlot(panc.filtered, reduction = "umap", group.by = "sex")+
    theme(text = element_text(size = 14,face="bold"),
          axis.text = element_text(size = 14,face="bold"))+
    labs(title = "Sex")

p<-grid.arrange(arrangeGrob( p2,p3,p4,p1, ncol=4))

# Supplemental Figure 1B----

#First address ambient expression of marker genes before plotting

DefaultAssay(panc.filtered)<-"RNA"
Idents(panc.filtered)<-'CellTypes'
VlnPlot(panc.filtered, features ="INS", pt.size=-1)
VlnPlot(panc.filtered, features ="GCG", pt.size=-1)
VlnPlot(panc.filtered, features ="SST", pt.size=-1)
VlnPlot(panc.filtered, features ="PRSS1", pt.size=-1)
VlnPlot(panc.filtered, features ="KRT19", pt.size=-1)
VlnPlot(panc.filtered, features ="PECAM1", pt.size=-1)
VlnPlot(panc.filtered, features ="PDGFRB", pt.size=-1)
VlnPlot(panc.filtered, features ="HLA-DRA", pt.size=-1)

# Based on bimodality of cell type markers, 
# set threshold on log exp to remove counts below the threshold 
# This is only for the purpose of cell annotation on featurePlots. Counts were not removed in the original gene expression matrix
amb<-GetAssayData(object = panc.filtered, slot = "data")
amb["INS",]<-ifelse(amb["INS",]<5,0,amb["INS",])
amb["GCG",]<-ifelse(amb["GCG",]<5,0,amb["GCG",])
amb["SST",]<-ifelse(amb["SST",]<5,0,amb["SST",])
amb["PRSS1",]<-ifelse(amb["PRSS1",]<4,0,amb["PRSS1",])
panc.filtered<-SetAssayData(panc.filtered, slot = "data", new.data = amb)

FeaturePlot(panc.filtered, features = "INS")+
    theme(text = element_text(size = 14,face="bold"),
          axis.text = element_text(size = 14,face="bold"))

FeaturePlot(panc.filtered, features = "GCG")+
    theme(text = element_text(size = 14,face="bold"),
          axis.text = element_text(size = 14,face="bold"))

FeaturePlot(panc.filtered, features = "SST")+
    theme(text = element_text(size = 14,face="bold"),
          axis.text = element_text(size = 14,face="bold"))

# similar plots for rest of the cell type markers were considered for cell type annotations

# Supplemental Figure 1C----

DefaultAssay(panc.filtered) <- "integrated"
Idents(panc.filtered) <- "CellTypes"

# Reorder cell type identity 
levels(x = panc.filtered)
my_levels <- c("Alpha", "Beta", "Delta-PP-Ep","Acinar","Ductal","Endothelial","Stellate","Immune","6","13","17","18","21")
panc.filtered@active.ident <- factor(x = panc.filtered@active.ident, levels = my_levels)

Pancreas <- c("GCG","IRX2","PCSK2","INS","IAPP","PDX1","SST","HHEX","LEPR","PPY","GHRL", "PRSS1","REG1A","CPA2", "KRT19","PROM1","PECAM1","PLVAP","VWF", "PDGFRB","TIMP1","COL1A1", "HLA-DRA","CD68","CD74")
genes<-rev(Pancreas)

DotPlot(panc.filtered, feature=genes, dot.scale = 8)+ 
    scale_color_gradient2(low = "#8080AD", high = "#BB3C34", mid = "white")+
    labs(y=" ", x="")+RotatedAxis()+
    theme(text = element_text(size = 15),
          axis.text = element_text(size = 15))

# Supplemental Figure 1D----
TF <- c("RFX6","PDX1","PAX6","NKX6-1","NKX2-2","NEUROD1","MAFB","MAFA","ISL1","IRX2","FOXA2","ARX")

DotPlot(panc.filtered, features=TF,dot.scale = 8, col.min = -1, col.max = 1.0)+ 
    coord_flip()+ 
    scale_color_gradient2(low = "#8080AD", high = "#BB3C34", mid = "white")+
    labs(y=" ", x="")+RotatedAxis()+
    theme(text = element_text(size = 15),
          axis.text = element_text(size = 15))





# Supplemental Figure 1E----
# plot age distribution among data source

# Create unique donor ID in metadata column 'donorid2'
meta<-panc.filtered@meta.data
meta$cellID<-rownames(meta)
meta<- meta %>% group_by(donor)%>% mutate(id=cur_group_id())
#check<-meta %>% group_by(donor,id) %>% tally()
meta<-meta %>% mutate(donorid2=paste("Donor", id, sep = "_")) %>% select(-id)
rownames(meta)<-meta$cellID
panc.filtered<-AddMetaData(panc.filtered, metadata = meta)

# Donor/disease cell counts
disease_count<-meta %>% group_by(status, donorid2) %>% count()
disease_count<-disease_count[,1] %>% count()
#68ND, 35diseased

# Distribution of age per decade only in healthy donors

# subset healthy donors
Idents(panc.filtered)<-"status"
healthy<-subset(panc.filtered, idents = "healthy")

# classify age per decade and add to metadata
meta<-healthy@meta.data
meta$cellID<-rownames(meta)
meta$age<-as.numeric(meta$age)
meta<-meta %>% mutate(age_breakby_10 = case_when(age < 10  ~ "<10",
                                                 age >= 10 & age < 20 ~ "10_19",
                                                 age >= 20 & age < 30 ~ "20_29",
                                                 age >= 30 & age < 40 ~ "30_39",
                                                 age >= 40 & age < 50 ~ "40_49",
                                                 age >= 50 & age < 60 ~ "50_59",
                                                 age >= 60 & age <= 76 ~ "60_76"
                                                 
                                                 ))

table(meta$age_breakby_10)
healthy <- AddMetaData(healthy, metadata = meta)

# plot

tb<-table(healthy$age_breakby_10, healthy$study)%>%
    as.data.frame()
tb %>% 
    ggplot(aes(x=Freq, y=Var2,fill=Var1)) + 
    #geom_bar(position="stack", stat="identity") +
    geom_bar(position="fill", stat="identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("%age group per dataset ")+
    xlab("%")+ylab("source of data")+
    theme(text = element_text(size = 18,face="bold"),
          axis.text = element_text(size = 18,face="bold"),
          axis.line.x = element_line(color="black", size = 0.8),
          axis.line.y = element_line(color="black", size = 0.8))


# Supplemental Figure 1F----

# cell type frequency by donors
celltype_counts<-as.data.frame(table(panc.filtered$donorid2,panc.filtered$CellTypes))
names(celltype_counts) <- c("donors","celltypes","freq")

# cell type frequency by age
age<-as.data.frame(table(panc.filtered$donorid2,panc.filtered$age))
names(age)<-c("donors","age","cell_freq")

#join age info with cell type freq per donor
celltype_counts <- celltype_counts %>% left_join(age, by="donors")

#get the order of descending age

donor_info <- meta %>% 
    select(age, donorid2) %>%
    group_by(donorid2) %>%
    count(age)%>%
    mutate(age2 = as.numeric(age))%>%
    arrange(desc(age2))

level_order <- donor_info$donorid2

celltype_counts %>% 
    ggplot(aes(x=freq, y=factor(donors, level = level_order),fill=celltypes)) + 
    geom_bar(position="fill", stat="identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    xlab("% cell types")+ylab("sample_id")+labs(fill='cell types') +
    theme(text = element_text(size = 18,face="bold"),
          axis.text = element_text(size = 18,face="bold"),
          axis.line.x = element_line(color="black", size = 0.8),
          axis.line.y = element_line(color="black", size = 0.8))


# Subset healthy beta cells only----
panc.filtered<-readRDS("data/cds/panc.filtered.Rds")
Beta<-subset(panc.filtered,idents = "Beta")
Idents(Beta)<-"status"
dim(Beta)

# subset only healthy cells
Healthy_Beta<-subset(Beta, idents = "healthy")
dim(Healthy_Beta) #19407 11279

# Scaling
all.genes <- rownames(Healthy_Beta)
Healthy_Beta <- ScaleData(Healthy_Beta, features = all.genes)

# Perform PCA
Healthy_Beta <- RunPCA(Healthy_Beta)
DimHeatmap(Healthy_Beta, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Healthy_Beta, dims = 1:12, cells = 500, balanced = TRUE)

# Visualize PCA plot
DimPlot(Healthy_Beta, 
        reduction = "pca", 
        label = T, 
        #group.by = "agegroup",
        dims = c(9,5),
        pt.size = 1)


# Determining PCs
ElbowPlot(Healthy_Beta)

# Clustering
Healthy_Beta <- FindNeighbors(Healthy_Beta, dims = 1:20)
Healthy_Beta <- FindClusters(Healthy_Beta, resolution = 0.5)
Healthy_Beta <- RunUMAP(Healthy_Beta, dims = 1:20)


Healthy_Beta[["cluster_res0.5"]]<-Idents(Healthy_Beta)
Healthy_Beta[["percent.ribo"]] <- PercentageFeatureSet(Healthy_Beta, pattern = '^RPL|^RPS' )

# this Healthy_Beta seurat object was used as input for pyscenic

# Convert to ann data obj for scanpy visuals for rest of Figure 1----
DefaultAssay(Healthy_Beta)<-"RNA"

library(SeuratDisk)
SaveH5Seurat(Healthy_Beta, filename = "beta_healthy4.h5Seurat")
Convert("beta_healthy.h5Seurat", dest = "h5ad")




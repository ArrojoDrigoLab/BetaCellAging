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

import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as pl

# load ann data object (was created frin script/02-Figure-1.R)
adata_RNA = sc.read_h5ad("beta_healthy4.h5ad")

# When converted from Seurat to h5ad using SeuratDisk. 
# adata.raw.var_names is different than adata.var_names.
# As a result, plotting will be problematic.The problem is really just that you're not able to directly edit adata.raw.var_names, 
# but you can work around by convert .raw to a temporary AnnData, change the var_name to match with the adata.var_names, and convert back to .raw
# https://github.com/theislab/scanpy/issues/1406

adata_tmp = adata_RNA.raw.to_adata()
adata_tmp.var_names = adata_RNA.var_names
adata_RNA.raw = adata_tmp

adata_RNA.obs["age_breakby10"] = adata_RNA.obs["age_breakby_10"].astype("category")

# Figure 1C ----
# Correlation matrix among age groups per decade
sc.tl.dendrogram(adata_RNA, groupby='age_break10')
sc.pl.correlation_matrix(adata_RNA, 'age_break10')


# Figure 1D ----

TF_betaenriched = {'Transcriptional factors' : ['MAFA','MAFB','PAX6','PDX1','NKX2-2','NKX6-1','NEUROD1','ISL1','FOXA2'],
                   'Beta_enriched' : ['INS','IAPP','HADH','DLK1','RBP4','NPTX2','GAD2','EIF4A2','SAMD11','G6PC2','GSN','ABCC8','ADCYAP1','HOPX','PCSK1']
                     }
sc.pl.dotplot(adata_RNA, TF_betaenriched, groupby='age_break10', dendrogram=True, use_raw=False, standard_scale='var',var_group_rotation=0)

# Figure 1E ----
beta_func_genes = {'Ion channels': ["KCNK16","KCNJ8","KCNJ6","CACNB2","CACNA2D1","CACNA1D","CACNA1C","CACNA1A","ABCC8","SLC30A8","ACLY"],
                     'Insulin secretion': ["GPX3","GPI","G6PC2","ENO2","PKM","PDHA1"],
                     'Vessicle organization': ["VAMP3","SYT7","SYT3","CREB1","EXOC8","RAB1A"], #https://www.cell.com/cell-metabolism/pdf/S1550-4131(17)30731-3.pdf FIGURE4
                     'Protein folding': ["B2M","HSP90AB1","DNAJB1","ATF6","BMP5"],
                    'Ubiquitin processing' : ['BIRC2','UBE2G2','CANX','PSMA4','PTEN','STYX'],
                     'ER Stress': ["HERPUD1","DDIT3","HSPA5","ATF4","XBP1","JUND","CREB3L2","ERN1","ATF6"]
                   }
            
sc.pl.dotplot(adata_RNA, TF_betaenriched, groupby='age_break10', dendrogram=True, use_raw=False, standard_scale='var',var_group_rotation=0)            

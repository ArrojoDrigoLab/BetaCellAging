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
from matplotlib.pyplot import rc_context

# load ann data object 
# seurat object subsetted to only beta cells(healthy donors) from script/02-Figure-1.R was converted to ann data object for following figures)

# This h5ad available in zenodo deposit, see link in Readme
adata_RNA = sc.read_h5ad("beta_healthy4.h5ad")

# When converted from Seurat to h5ad using SeuratDisk. 
# adata.raw.var_names is different than adata.var_names.
# As a result, plotting will be problematic.The problem is really just that you're not able to directly edit adata.raw.var_names, 
# but you can work around by convert .raw to a temporary AnnData, change the var_name to match with the adata.var_names, and convert back to .raw
# https://github.com/theislab/scanpy/issues/1406

adata_tmp = adata_RNA.raw.to_adata()
adata_tmp.var_names = adata_RNA.var_names
adata_RNA.raw = adata_tmp

adata_RNA.obs.dtypes #seurat clusters are considered integers ,convert to category

adata_RNA.obs["age_breakby10"] = adata_RNA.obs["age_breakby_10"].astype("category")
adata_RNA.obs["cluster_res0.5"] = adata_RNA.obs["cluster_res0.5"].astype("category")

# DE genes per healthy beta cluster

# Figure 2A ----
with rc_context({'figure.figsize': (6, 5)}):
    sc.pl.umap(adata_RNA, color='cluster_res0.5', add_outline=True,
               legend_fontsize=12, legend_fontoutline=2,frameon=False,title='', legend_loc = 'on data')
               
               
color_dict = dict({'<10':'#E69F00',
                  '10_19':'#56B4E9',
                  '20_29': '#009E73',
                  '30_39': '#F0E442',
                   '40_49': '#0072B2',
                   '50_59': '#D55E00',
                   '60_76': '#CC79A7'})
adata_RNA.obs['age_break10'].cat.reorder_categories(['<10','10_19','20_29','30_39', '40_49', '50_59','60_76'], inplace = True)

# Figure 2B ----
# UMAP of healthy beta clusters
with rc_context({'figure.figsize': (6, 5)}):
    sc.pl.umap(adata_RNA, color='age_breakby_10', add_outline=True,
               legend_fontsize=12, legend_fontoutline=2,frameon=False,title='', palette=color_dict)

# Figure 2C----
marker_genes_dict = {'0': ["INS",	"PCSK1N",	"EIF3E",	"MT1X",	"MT-ND3",	"MTRNR2L1"],
                     '1': ["CD99",	"MAFB"],
                     '2': ["TUBA4A",	"PPP1R1A",	"TMEM14A",	"TAGLN2",	"YWHAH",	"PCBP1",	"PKIB",	"TUBB",	"TXNL4A",	"TUBB4B", "DLK1"],
                     '3': [	"TMSB4X",	"TFF3",	"CMPK1",	"CPB1",	"AQP3",	"TMSB10",	"IGFBP5",	"CDKN1A",	"ID3",	"ID1"],
                     '4': ["RPS4X",	"RPL8",	"MT-CO3",	"RPL5",	"RPL6",	"TPT1","NACA","FKBP2","SNRPN",	"RBP1"],
                     '5': ["SERF2",	"PCP4",	"FTH1",	"RPS2",	"NENF",	"S100A6",	"PFN1",	"TMSB10",	"PPDPF"],
                     '6': ["FTH1",	"RBP4",	"CD63",	"PPDPF",	"MYL6",	"RPS29",	"CLU",	"HLA-C",	"PARK7"],
                     '7': ["PDZD4",	"PPP2R3B",	"CORO2A",	"SPRTN",	"ADCY1",	"TTC7A",	"NAF1",	"LRRC7"],
                     '8': ["MT-CO1",	"MT-CO3",	"MT-CYB",	"MT-ND4",	"CHGB",	"HLA-A",	"CLU",	"MT-ND1",	"PCSK2",	"FXYD3"],
                     '9': ["HSPA5",	"DUSP1",	"DNTTIP2",	"KLF6",	"WTAP",	"PPP1CB",	"TMED4",	"EZR",	"PSMG3",	"HSPB1"],
                     '10': ["DDIT3",	"EIF1",	"PPP1R15A",	"MAP1LC3B",	"SLC3A2",	"RHEB",	"GADD45A",	"SQSTM1",	"PFDN2",	"EIF5",],
                     '11': ["RBPJL",	"SPARCL1",	"SLC35D1",	"TM9SF2",	"IAPP",	"PCSK1",	"ENPP2",	"SYT4",	"DDHD1",	"G6PC2"],
                     '12': ["BCL6B",	"TRPV4",	"HLA-DRB1",	"PTN",	"PRTN3",	"AURKC",	"FAM83H",	"ATP10B",	"ZMYND12",	"EHD2"]
                     }
                     
mp=sc.pl.matrixplot(adata_RNA, marker_genes_dict, groupby='cluster_res0.5', dendrogram=False, use_raw=False,return_fig=True,standard_scale='var', cmap='viridis',figsize=(30,4))
mp.add_totals().style(edge_color='black').show()

#Figure 2D ----

marker_genes_dict2 = {'0': ["INS",	"MT1X"],
                     '1': ["CD99",	"MAFB"],
                     '2': ["TUBA4A",	"PKIB"],
                     '3': ["CDKN1A",	"ID3", "TFF3"],
                     '4': ["FKBP2",	"RBP1"],
                     '5': [	"PCP4",	"FTH1"],
                     '6': ["RBP4",	"CD63"],
                     '7': [	"CORO2A",	"ADCY1"],
                     '8': [	"CHGB",	"MT-CO1"],
                     '9': ["HSPA5",	"HSPB1"],
                     '10': ["DDIT3",	"GADD45A"],
                     '11': ["IAPP",	"PCSK1"],
                     '12': ["TRPV4",	"EHD2"]
                     }
mp=sc.pl.matrixplot(adata_RNA, marker_genes_dict2, groupby='age_break10', dendrogram=False, use_raw=False,return_fig=True,standard_scale='var', cmap='viridis',figsize=(20,4))
mp.add_totals().style(edge_color='black').show()


# Supplemental Figure 2C ----
# Dotplot for disallowed genes:

#list of disallowed genes for beta cells listed herhe-> (https://www.frontiersin.org/files/Articles/249155/fgene-08-00041-HTML/image_m/fgene-08-00041-g001.jpg)
disallowed_genes = ['LDHA','CDKN1C','IGF1R','HSD11B1','OLFML1','IGFBP4','ZCCHC24','SMOC2','NDRG2','FCGRT','LRIG3','MGLL','RASGRP2','FAM13A','DAPK2','HPGD','YAP1','TPM2','RARRES2','GAS1','IGF1','SLC16A1','GLP1R','GIPR','ADCY5','EIF2AK4','EIF2AK3']
sc.pl.dotplot(adata_RNA, disallowed_genes, groupby='age_break10', dendrogram=True, use_raw=False, standard_scale='var',var_group_rotation=0)





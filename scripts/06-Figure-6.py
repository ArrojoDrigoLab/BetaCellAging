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

import os
import pandas as pd # v 1.2.4
import scanpy as sc #v 1.7.2
import numpy as np #1.20.2
import scipy as scipy
import seaborn as sns
import matplotlib.pyplot as plt #v 3.4.1
from matplotlib.pyplot import rc_context

# below h5ad file is same as in 03-Figure-3.py and 02-Figure-1.py except pySCENIC regulon binary matrix is added instead of gene expression matrix
adata_BIN = sc.read_h5ad("beta_healthy2_BIN.h5ad")

adata_BIN

# Fix for seurat-to-anndata conversion screwup
adata_tmp = adata_BIN.raw.to_adata()
adata_tmp.var_names = adata_BIN.var_names
adata_BIN.raw = adata_tmp


# PAGA's follows proceeds in two steps: 
# 1. cluster the single-cell graph 
# 2. construct the abstracted graph based on the connectivity between the clusters.

# 1. cluster the single cell graph

# Compute PCA/neighbor graphs
import random
random.seed(30)
sc.tl.pca(adata_BIN, svd_solver='arpack')
sc.pp.neighbors(adata_BIN, n_neighbors=10, n_pcs=20)

#leiden clustering
sc.tl.leiden(adata_BIN, resolution= 0.4)

#reannotate
adata_BIN.obs['leiden_anno'] = adata_BIN.obs['leiden']

#change to right data types
adata_BIN.obs["age_breakby_10"] = adata_BIN.obs["age_breakby_10"].astype("category")
adata_BIN.obs["cluster_res0.5"] = adata_BIN.obs["cluster_res0.5"].astype("category") #this is seurat clusters from RNA

# 1. Dimension reduction visuals
# Draw_graph: An alternative to tSNE 
sc.tl.draw_graph(adata_BIN)
sc.pl.draw_graph(adata_BIN, color=['leiden_anno','age_breakby_10', 'cluster_res0.5'], legend_loc='on data')


#2. Change to right colors and orders of the metadata
adata_BIN.obs['age_breakby_10'].cat.reorder_categories(['<10','10_19','20_29','30_39', '40_49', '50_59', '60_76'], inplace = True)
cols = ["#E69F00" ,"#56B4E9", "#009E73", "#F0E442","#0072B2", "#D55E00", "#CC79A7"]

cols = np.array(cols)
new_colors = np.array(adata_BIN.uns['age_breakby_10_colors'])
new_colors[[0]] = cols[[0]]
new_colors[[1]] = cols[[1]]
new_colors[[2]] = cols[[2]]
new_colors[[3]] = cols[[3]]
new_colors[[4]] = cols[[4]]
new_colors[[5]] = cols[[5]]
new_colors[[6]] = cols[[6]]

adata_BIN.uns['age_breakby_10_colors'] = new_colors


# 3. PAGA
sc.tl.paga(adata_BIN, groups='leiden_anno')
sc.pl.paga(adata_BIN, color=['leiden_anno'],cmap='viridis')


# 4. Recomputing the embedding using PAGA-initialization
sc.tl.draw_graph(adata_BIN, init_pos='paga')


# Figure 6A and 6B----
#draw graph-paga initialized visuals
sc.pl.draw_graph(adata_BIN, color=['leiden_anno', 'age_breakby_10'],frameon=False)

# PAGA graphs
sc.pl.paga(adata_BIN,color=['leiden_anno','age_breakby_10'],title='PAGA grapgh',fontoutline=None,threshold=0,
           edge_width_scale=0.7, node_size_scale=1.8, node_size_power=1, fontsize=14, solid_edges='connectivities_tree',
           dashed_edges='connectivities', left_margin=0.01, frameon=False)

# Figure 6C ----
# PAGA tree layout , root 4 chosen(<10 agegroup)
sc.pl.paga(adata_BIN,color=['leiden_anno','age_breakby_10'],root=[4],title='',colorbar=True,fontoutline=None,threshold=0,
               edge_width_scale=0.7,node_size_scale=1.8,node_size_power=1,fontsize=14,solid_edges='connectivities_tree',
               dashed_edges='connectivities',left_margin=0.01,frameon=False,layout='rt')


# Supplemental Figure 6B ----

# Add pseudotime
adata_BIN.uns['iroot'] = np.flatnonzero(adata_BIN.obs['leiden_anno']  == '4')[0] #Choose a root cell for diffusion pseudotime.
sc.tl.dpt(adata_BIN)
sc.pl.draw_graph(adata_BIN, color='dpt_pseudotime', cmap="BuPu")

adata_BIN.obs['distance'] = adata_BIN.obs['dpt_pseudotime']

# Reconstructing gene changes along PAGA paths for a given set of genes
# Select some of the marker gene names to show in paga_path heatmap
pyscenic_kmean_clus = ['KLF4','NFATC2','FOXO6','NKX6-1','RXRG','MAFA','PDX1','MYCN','LHX4','JUND','XBP1','ATF6','ATF4','NEUROD1','FOXA2','RFX2','NPDC1'
                       ,'CEBPD','ETV2','YBX1','RFX3','FOSB','HIF1A','KLF3','FOXA1','TAF7','NELFE','CREB3L2','HNF4A']

# Need log conversion and compute var genes below to be able to run paga_path
# Use the full raw data for visualization.
sc.pp.log1p(adata_BIN)
sc.pp.highly_variable_genes(adata_BIN)


# Choose path based on nodes trajectory of sc.paga.pl
paths = [('Path1', [4,5,1,3,]),
         ('Path2', [4,5,1,0]),
         ('Path3', [4,5,9]),
         ('Path4', [4,5,7,2,8]),
         ('Path5', [4,5,7,2,6]),
         ('Path6', [1,5,2])]


paths = [('Path6', [1,5,2])]
#plot paga_path
_, axs = plt.subplots(ncols=6, figsize=(16, 12), gridspec_kw={'wspace': 0.1, 'left': 0.12})
plt.subplots_adjust(left=0.02, right=0.98, top=0.82, bottom=0.2)
for ipath, (descr, path) in enumerate(paths):
    _, data = sc.pl.paga_path(
        adata_BIN, path, pyscenic_kmean_clus,
        show_node_names=True,
        ax=axs[ipath],
        ytick_fontsize=12,
        left_margin=0.2,
        n_avg=50,
        show_colorbar=True,
        annotations=['distance'],
        title_fontsize= 18,
        normalize_to_zero_one = True,
        show_yticks=True if ipath==0 else False,
        color_map=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True),
        color_maps_annotations={'distance': 'BuPu'},
        title='{} path'.format(descr),
        return_data=True,
        show=False)
pl.show()





# Compute differential expression among leiden clusters

sc.tl.rank_genes_groups(adata_BIN, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata_BIN, n_genes=10, sharey=False)
pd.DataFrame(adata_BIN.uns['rank_genes_groups']['names']).head(10)

#Get a table with the scores and groups.
result = adata_BIN.uns['rank_genes_groups']
groups = result['names'].dtype.names
data=pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(10)

sc.tl.rank_genes_groups(adata_BIN, 'leiden', groups=['5'], reference='4', method='wilcoxon')
sc.pl.rank_genes_groups(adata_BIN, groups=['5'], n_genes=20)

sc.tl.rank_genes_groups(adata_BIN, 'leiden', groups=['5'], reference='0', method='wilcoxon')
sc.pl.rank_genes_groups(adata_BIN, groups=['5'], n_genes=20)

#export metadata
adata_BIN.obs.to_csv('metadata.csv')



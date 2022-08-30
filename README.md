# Aging compromises human islet beta cell function and identity by decreasing transcription factor activity and inducing ER stress

Shristi Shrestha<sup>1</sup> , Galina Erikson<sup>2</sup>, James Lyon<sup>4</sup>, Aliya F Spigelman<sup>4</sup>, Austin Bautista<sup>4</sup>, Jocelyn E Manning Fox<sup>4</sup>, Maxim Shokhirev<sup>2</sup>, Jean-Phillippe Cartailler<sup>1</sup>, Martin W. Hetzer<sup>3</sup>, Patrick E. MacDonald<sup>4</sup>, Rafael Arrojo e Drigo<sup>5</sup>

> <sup>1 – Creative Data Solutions, Vanderbilt Center for Stem Cell Biology, Nashville, Tennessee, USA \
2 – Integrative Genomics and Bioinformatics Core, Salk Institute of Biological Studies, La Jolla, CA, USA 92037 \
3 – Molecular and Cell Biology Laboratory, Salk Institute of Biological Studies, La Jolla, CA, USA 92037 \
4 – Department of Pharmacology and Alberta Diabetes Institute, University of Alberta, Edmonton, Canada, T6G2E1 \
5 – Department of Molecular Physiology and Biophysics, Vanderbilt University, Nashville, Tennessee, USA</sup>



![main_figure](main_figure.png)

This repository contains scripts and links to analysis data used to generate manuscript figures.

## Scripts

The `script` folder contain both R and python script organized by the figure numbers in the manuscript.

## Analysis Data

Seurat object containing annoted cell types, clusters, donor metadata, pySCENIC regulon binary matrix, pySCENIC AUC scores and scRNA-seq gene expression are posted to [Zenodo](https://zenodo.org/deposit/6491944), including additional tabular data in CSV format:

  1. `panc.filtered.Rds` - Seurat object for post-processed, merged and integrated dataset.
  2. `allcelltype_healthy_13000.pySCENIC` - Seurat object for subset of cells (13,000 randomly selected) from healthy donors only, containing pySCENIC-generated regulon AUC score and binary assignments separate assays.
  3. `beta_healthy4.h5ad` - Anndata object for beta cell subset from healthy donors only.
  4. `beta_healthy2_BIN.h5ad` - Anndata object for beta cell subset from healthy donors only that includes pySCENIC-generated regulon binary assignments.
  5. Tabular data used as inputs for several figures, generated from a pySCENIC run on various agegroup and cell type subsets. See scripts for respective figure script documentation for details on these CSV files.

## Raw Data

Raw files of data analyzed in this manuscript were originally retrieved from::

  1. GEO accession number [GSE81547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81547) (Enge et al 2017)
  2. GEO accession number [GSE81608](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81608) (Xin et al 2016)
  3. ArrayExpress accession [E-MTAB-5061](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/) (Segerstolpe et al 2016)
  4. GEO accession number [GSE114297](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114297) (Xin et al 2018) 
  5. GEO accession number [GSE124742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124742) (Camunas et al 2020)
  6. European Genome-Phenome Archive [EGAS00001004653](https://ega-archive.org/studies/EGAS00001004653) (Tosti et al 2020)
  7. [Human Pancreas Analysis Program](https://hpap.pmacs.upenn.edu/)



## Donor Metadata

| Data Source            | Donors | Disease                                                      | Cells |
| ---------------------- | ------ | ------------------------------------------------------------ | ----- |
| Tosti et al 2020       | 6      | 4 healthy, 2 pancreatitis                                    | 13182 |
| Xin et al 2018         | 12     | 12 healthy                                                   | 20157 |
| HPAP                   | 16     | 6 AA+, 4 healthy, 2 T1D, 4 T2D                               | 4191  |
| Camunas et al 2020     | 35     | 1 elevated HbA1c, 22 healthy, 1 prediabetic, 3 T1D, 7 T2D, 1 T2Dreversed | 5095  |
| Enge et al 2017        | 8      | 8 healthy                                                    | 2327  |
| Segerstolpe et al 2016 | 10     | 6 healthy, 4 T2D                                             | 2501  |
| Xin et al 2016         | 18     | 12 healthy, 6 T2D                                            | 1501  |






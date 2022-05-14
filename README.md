# MCDA

This repository contains code to cluster and generate figures for "Systematic identification of cell fate regulatory programs using a single-cell atlas of mouse development".

## MCDA_construction

**MCDA_construction.ipynb**: The script was used to generate the Mouse cell differentiation atlas.


## Entropy_calculation

**CCAT_calculation_using_homologous_genes.R**: The script was utilized to perform single-cell entropy estimation using CCAT method with a weighted matrix to leverage all the homology genes between human and other species. Hydra was used as an example. The gene homology relationship used for weighted homology matrixes between human and other species were given in "MCDA/homologous _genes/all_all" folder.

**Entropy_calculation_using_three_methods.R**: The script was utilized to perform single-cell entropy estimation using CCAT, StemID, and SLICE methods between Xbp1-/- and wild-type cells per cluster.


## commonly_regulated_genes

**function.R**: The script contains some functions that were needed in the following R script.

**Hydra_cell_type_pairs.ipynb**: The script is used for PAGA anaylsis per lineage among invertebrates and Hydra vulgaris (Hydra) was selected as an example. Edges in an abstracted graph with a probability higher than 0.0005 were considered as possible connections of cell-type hierarchies. 


**Hydra_differentially_expressed_genes.R**: The script is used for perform differential expression analysis for cells in each lineage-cell type separately according to the cell type hierarchy of invertebrates using FindMarkers function in Seurat. Hydra was used as an example.


**Human_cell_type_pairs.ipynb": The script is used for connecting cell states across time according to gene expression similarity of cell types in vertebrates. Human was used as an example.


**Human_differentially_expressed_genes.R**: The script is used for perform differential expression analysis for cells in each tissue-cell type separately according to the cell type hierarchy of vertebrates using FindMarkers function in Seurat. Human was used as an example.


**commonly_regulated_genes.R**: The script is used to generating the commonly regulated genes in seven species. The one to one gene homology relationship between human and other species were given in "MCDA/homologous _genes/one_one" folder.



## homologous _genes

**MCDA/homologous _genes/one_one**: This folder contains the gene homology relationship used for weighted homology matrixes between human and other seven species were given.
 
**MCDA/homologous _genes/all_all**: This folder contains the one to one gene homology relationship between human and other seven species.



############################################################################## Figure 1c
##################### 
import scanpy as sc
import pandas as pd
import os
os.chdir("./merge/scanpy")


#### input
adata = sc.read_h5ad("MCDA_raw.h5ad")

#### filter
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=3)
adata.obs['n_counts'] = adata.X.sum(axis=1)


#### Feature selection
filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=0.005, max_mean=15, min_disp=0.5)
import collections
collections.Counter(filter_result.gene_subset)
#Counter({False: 31173, True: 3039})
sc.pl.filter_genes_dispersion(filter_result)
adata = adata[:, filter_result.gene_subset]

####  Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.
sc.pp.regress_out(adata, ['n_counts'])

####  scale the data
sc.pp.scale(adata, max_value=10)

####  PCA
sc.tl.pca(adata, n_comps=100)
sc.pl.pca_loadings(adata)
adata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
sc.pl.pca_scatter(adata, color='Col1a1')
# PC
sc.pl.pca_variance_ratio(adata, log=True,  show=100,n_pcs=100)sc.pp.neighbors(adata, n_neighbors=10,n_pcs=50)
sc.tl.leiden(adata, resolution=2)

#### TSNE
sc.tl.tsne(adata,n_pcs=50)
sc.pl.tsne(adata, color=['tissue'])

#### Clutering
sc.tl.leiden(adata, resolution=8)
sc.pl.tsne(adata, color=['leiden'])
adata.obs["leiden8"]=adata.obs["leiden"]
sc.pl.tsne(adata, color=['leiden'],legend_loc="on data")
sc.pl.tsne(adata, color=['stage'],save="_stage_tsne.pdf")
sc.pl.tsne(adata, color=['tissue1'],save="_tissue_tsne_final.pdf")


new=range(1,96)
adata.rename_categories('leiden_final',new)
sc.set_figure_params(dpi=200)
sc.pl.tsne(adata, color=['leiden_final'],legend_loc="on data",palette=sc.pl.palettes.godsnot_102,legend_fontsize=6)


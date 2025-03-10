RNA velocity with scVelo
================
Aaron Mohammed

``` r
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(velocyto.R)
library(ggplot2)

seurat_input_dir <- file.path(dirname(getwd()),
                              "2_Annotation_and_integration")

loom_input_dir <- file.path(getwd(), "velocyto_output")
```

``` r
merged_s <- readRDS(file.path(seurat_input_dir, "final_cardiac_seurat.rds"))

myo_s <- readRDS(file.path(seurat_input_dir, "final_myo_seurat.rds"))

WT.ldat <- ReadVelocity(file = file.path(loom_input_dir, "WT", "Wildtype.loom"))
KO.ldat <- ReadVelocity(file = file.path(loom_input_dir, "KO", "Knockout.loom"))
```

``` r
insert.ldat <- function (ldat, seurat, patt, rep) {
  for (i in names(x = ldat)) {
    ### Store assay in a new variable
    assay <- ldat[[i]]
  
    ### Rename cell names in loom file to match cell names in Seurat object
    colnames(assay) <- gsub(patt, rep, colnames(assay))
    colnames(assay) <- gsub("x", "-1", colnames(assay))
  
    ### Subset to filtered cells in Seurat object
    assay <- assay[,colnames(seurat)]
  
    ### Add assay to Seurat object
    seurat[[i]] <- CreateAssayObject(counts = assay)
  }
  return(seurat)
}
```

``` r
DefaultAssay(merged_s) <- "RNA"
Idents(merged_s) <- "orig.ident"

WT <- subset(merged_s, idents = "WT")

WT <- insert.ldat(ldat = WT.ldat,
                      seurat = WT,
                      patt = "Wildtype:",
                      rep = "WT_")

SaveH5Seurat(WT, filename = "WT.h5Seurat", overwrite=TRUE)
Convert("WT.h5Seurat", dest = "h5ad", overwrite=TRUE)

################################################################################

KO <- subset(merged_s, idents = "KO")

KO <- insert.ldat(ldat = KO.ldat,
                      seurat = KO,
                      patt = "Knockout:",
                      rep = "KO_")

SaveH5Seurat(KO, filename = "KO.h5Seurat", overwrite=TRUE)
Convert("KO.h5Seurat", dest = "h5ad", overwrite=TRUE)

################################################################################

merged <- merge(x = WT, 
                y = KO,
                  add.cell.ids = c("WT",
                                   "KO"),
                  project = "scRNAseq_Ren")

SaveH5Seurat(merged, filename = "merged.h5Seurat", overwrite=TRUE)
Convert("merged.h5Seurat", dest = "h5ad", overwrite=TRUE)
```

``` r
DefaultAssay(myo_s) <- "RNA"
Idents(myo_s) <- "orig.ident"

WT_myo <- subset(myo_s, idents = "WT")

WT_myo <- insert.ldat(ldat = WT.ldat,
                      seurat = WT_myo,
                      patt = "Wildtype:",
                      rep = "WT_")

SaveH5Seurat(WT_myo, filename = "WT_myo.h5Seurat", overwrite=TRUE)
Convert("WT_myo.h5Seurat", dest = "h5ad", overwrite=TRUE)

################################################################################

KO_myo <- subset(myo_s, idents = "KO")

KO_myo <- insert.ldat(ldat = KO.ldat,
                      seurat = KO_myo,
                      patt = "Knockout:",
                      rep = "KO_")

SaveH5Seurat(KO_myo, filename = "KO_myo.h5Seurat", overwrite=TRUE)
Convert("KO_myo.h5Seurat", dest = "h5ad", overwrite=TRUE)

################################################################################

myo <- merge(x = WT_myo, 
             y = KO_myo,
             add.cell.ids = c("WT",
                              "KO"),
             project = "scRNAseq_Ren")

SaveH5Seurat(myo, filename = "myo.h5Seurat", overwrite=TRUE)
Convert("myo.h5Seurat", dest = "h5ad", overwrite=TRUE)
```

``` r
library(reticulate)
use_condaenv("RNAvelocity")  
```

``` python
import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from matplotlib import pyplot as plt
```

``` python
WT_adata = sc.read_h5ad('/Users/asm/Creighton/core_lab/scRNAseq_Bubr1_Ren_3/4_RNA_velocity/WT.h5ad')

scv.pp.filter_and_normalize(WT_adata)
scv.pp.moments(WT_adata, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(WT_adata,n_jobs=7)
scv.tl.velocity(WT_adata, mode='dynamical')
scv.tl.velocity_graph(WT_adata)

scv.pl.proportions(WT_adata, groupby='Broad_Annotation', save='WT_proportions.pdf')

scv.pl.velocity_embedding_grid(WT_adata, basis='umap', title='WT - Dynamical', scale=0.4,color='Broad_Annotation',palette = ("#00BF7D", "#A3A500", "#00B0F6", "#F8766D","#E76BF3"),groups=("Endocardial", "Epicardial", "Mesenchymal", "MP", "Myocardial"), save='WT_dynamical.svg')

scv.pl.velocity_embedding_stream(WT_adata,basis='umap', color='Broad_Annotation', title='WT - Dynamical', density=2.5,palette = ("#00BF7D", "#A3A500", "#00B0F6", "#F8766D","#E76BF3"),groups=("Endocardial", "Epicardial", "Mesenchymal", "MP", "Myocardial"), save='WT_dynamical_stream.svg', add_margin=0.05)

scv.tl.velocity_confidence(WT_adata)
scv.pl.scatter(WT_adata, c=('velocity_length', 'velocity_confidence'), cmap='coolwarm', perc=[5, 95],save= 'WT_velocity_length_confidence_dynamical.svg')

scv.tl.latent_time(WT_adata)
scv.pl.scatter(WT_adata, color='latent_time', color_map='gnuplot', size=80,title='WT - Latent Time', save= 'WT_latent_time.svg')

################################################################################

KO_adata = sc.read_h5ad('/Users/asm/Creighton/core_lab/scRNAseq_Bubr1_Ren_3/4_RNA_velocity/KO.h5ad')

scv.pp.filter_and_normalize(KO_adata)
scv.pp.moments(KO_adata, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(KO_adata,n_jobs=7)
scv.tl.velocity(KO_adata, mode='dynamical')
scv.tl.velocity_graph(KO_adata)

scv.pl.proportions(KO_adata, groupby='Broad_Annotation', save='KO_proportions.pdf')

scv.pl.velocity_embedding_grid(KO_adata, basis='umap', title='KO - Dynamical', scale=0.4,color='Broad_Annotation',palette = ("#00BF7D", "#A3A500", "#00B0F6", "#F8766D","#E76BF3"),groups=("Endocardial", "Epicardial", "Mesenchymal", "MP", "Myocardial"), save='KO_dynamical.svg')

scv.pl.velocity_embedding_stream(KO_adata,basis='umap', color='Broad_Annotation', title='KO - Dynamical', density=2.5,palette = ("#00BF7D", "#A3A500", "#00B0F6", "#F8766D","#E76BF3"),groups=("Endocardial", "Epicardial", "Mesenchymal", "MP", "Myocardial"), save='KO_dynamical_stream.svg', add_margin=0.05)

scv.tl.velocity_confidence(KO_adata)
scv.pl.scatter(KO_adata, c=('velocity_length', 'velocity_confidence'), cmap='coolwarm', perc=[5, 95],save= 'KO_velocity_length_confidence_dynamical.svg')

scv.tl.latent_time(KO_adata)
scv.pl.scatter(KO_adata, color='latent_time', color_map='gnuplot', size=80,title='KO - Latent Time', save= 'KO_latent_time.svg')

################################################################################

merged_adata = sc.read_h5ad('/Users/asm/Creighton/core_lab/scRNAseq_Bubr1_Ren_3/4_RNA_velocity/merged.h5ad')

scv.pp.filter_and_normalize(merged_adata)
scv.pp.moments(merged_adata, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(merged_adata,n_jobs=7)
scv.tl.velocity(merged_adata, mode='dynamical')
scv.tl.velocity_graph(merged_adata)

scv.pl.proportions(merged_adata, groupby='Broad_Annotation', save='merged_proportions.pdf')

scv.pl.velocity_embedding_grid(merged_adata, basis='umap', title='Merged - Dynamical', scale=0.4,color='Broad_Annotation',palette = ("#00BF7D", "#A3A500", "#00B0F6", "#F8766D","#E76BF3"),groups=("Endocardial", "Epicardial", "Mesenchymal", "MP", "Myocardial"), save='merged_dynamical.svg')

scv.pl.velocity_embedding_stream(merged_adata,basis='umap', color='Broad_Annotation', title='Merged - Dynamical', density=2.5,palette = ("#00BF7D", "#A3A500", "#00B0F6", "#F8766D","#E76BF3"),groups=("Endocardial", "Epicardial", "Mesenchymal", "MP", "Myocardial"), save='merged_dynamical_stream.svg', add_margin=0.05)

scv.tl.velocity_confidence(merged_adata)
scv.pl.scatter(merged_adata, c=('velocity_length', 'velocity_confidence'), cmap='coolwarm', perc=[5, 95],save= 'merged_velocity_length_confidence_dynamical.svg')

scv.tl.latent_time(merged_adata)
scv.pl.scatter(merged_adata, color='latent_time', color_map='gnuplot', size=80,title='Merged - Latent Time', save= 'merged_latent_time.svg')

################################################################################

adatas = [WT_adata, KO_adata]
broad_adata = ad.concat(adatas,  join="outer")

scv.tl.rank_velocity_genes(broad_adata, groupby='orig.ident', min_corr=.3)

df = scv.DataFrame(broad_adata.uns['rank_velocity_genes']['names'])
df.head()

df.to_csv(path_or_buf='broad_top_dynamical_genes.csv')

scv.pl.scatter(broad_adata, df['WT'][:10],nrows = 2, frameon=True, color='orig.ident', size=80, linewidth=1.5, figsize = (7,8),fontsize=20,legend_loc="right", save= 'WT_top10_dynamical_1.pdf')

scv.pl.scatter(broad_adata, df['WT'][:10],nrows = 2, frameon=True, color='Broad_Annotation', size=80, linewidth=1.5, figsize = (7,8),fontsize=20,legend_loc="right", save= 'WT_top10_dynamical_2.pdf')

scv.pl.scatter(broad_adata, df['KO'][:10],nrows = 2, frameon=True, color='orig.ident', size=80, linewidth=1.5, figsize = (7,8),fontsize=20,legend_loc="right", save= 'KO_top10_dynamical_1.pdf')

scv.pl.scatter(broad_adata, df['KO'][:10],nrows = 2, frameon=True, color='Broad_Annotation', size=80, linewidth=1.5, figsize = (7,8),fontsize=20,legend_loc="right", save= 'KO_top10_dynamical_2.pdf')
```

``` python
WT_myo_adata = sc.read_h5ad('/Users/asm/Creighton/core_lab/scRNAseq_Bubr1_Ren_3/4_RNA_velocity/WT_myo.h5ad')

scv.pp.filter_and_normalize(WT_myo_adata)
scv.pp.moments(WT_myo_adata, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(WT_myo_adata,n_jobs=7)
scv.tl.velocity(WT_myo_adata, mode='dynamical')
scv.tl.velocity_graph(WT_myo_adata)

scv.pl.proportions(WT_myo_adata, groupby='myo_annot', save='WT_myo_proportions.pdf')

scv.pl.velocity_embedding_grid(WT_myo_adata, basis='umap', title='WT - Dynamical', scale=0.4,color='myo_annot',palette = ("#A3A500","#F8766D","#00BF7D","#00B0F6", "#E76BF3"),groups=("Atrial","AVC","OFT","SAN", "Ventricular"), save='WT_myo_dynamical.svg')

scv.pl.velocity_embedding_stream(WT_myo_adata,basis='umap', color='myo_annot', title='WT - Dynamical', density=2.5,palette = ("#A3A500","#F8766D","#00BF7D","#00B0F6", "#E76BF3"),groups=("Atrial","AVC","OFT","SAN", "Ventricular"), save='WT_myo_dynamical_stream.svg', add_margin=0.05)

scv.tl.velocity_confidence(WT_myo_adata)
scv.pl.scatter(WT_myo_adata, c=('velocity_length', 'velocity_confidence'), cmap='coolwarm', perc=[5, 95],save= 'WT_myo_velocity_length_confidence_dynamical.svg')

scv.tl.latent_time(WT_myo_adata)
scv.pl.scatter(WT_myo_adata, color='latent_time', color_map='gnuplot', size=80,title='WT - Latent Time', save= 'WT_myo_latent_time.svg')

################################################################################

KO_myo_adata = sc.read_h5ad('/Users/asm/Creighton/core_lab/scRNAseq_Bubr1_Ren_3/4_RNA_velocity/KO_myo.h5ad')

scv.pp.filter_and_normalize(KO_myo_adata)
scv.pp.moments(KO_myo_adata, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(KO_myo_adata,n_jobs=7)
scv.tl.velocity(KO_myo_adata, mode='dynamical')
scv.tl.velocity_graph(KO_myo_adata)

scv.pl.proportions(KO_myo_adata, groupby='myo_annot', save='KO_myo_proportions.pdf')

scv.pl.velocity_embedding_grid(KO_myo_adata, basis='umap', title='KO - Dynamical', scale=0.4,color='myo_annot',palette = ("#A3A500","#F8766D","#00BF7D","#00B0F6", "#E76BF3"),groups=("Atrial","AVC","OFT","SAN", "Ventricular"), save='KO_myo_dynamical.svg')

scv.pl.velocity_embedding_stream(KO_myo_adata,basis='umap', color='myo_annot', title='KO - Dynamical', density=2.5,palette = ("#A3A500","#F8766D","#00BF7D","#00B0F6", "#E76BF3"),groups=("Atrial","AVC","OFT","SAN", "Ventricular"), save='KO_myo_dynamical_stream.svg', add_margin=0.05)

scv.tl.velocity_confidence(KO_myo_adata)
scv.pl.scatter(KO_myo_adata, c=('velocity_length', 'velocity_confidence'), cmap='coolwarm', perc=[5, 95],save= 'KO_myo_velocity_length_confidence_dynamical.svg')

scv.tl.latent_time(KO_myo_adata)
scv.pl.scatter(KO_myo_adata, color='latent_time', color_map='gnuplot', size=80,title='KO - Latent Time', save= 'KO_myo_latent_time.svg')

################################################################################

merged_myo_adata = sc.read_h5ad('/Users/asm/Creighton/core_lab/scRNAseq_Bubr1_Ren_3/4_RNA_velocity/myo.h5ad')

scv.pp.filter_and_normalize(merged_myo_adata)
scv.pp.moments(merged_myo_adata, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(merged_myo_adata,n_jobs=9)
scv.tl.velocity(merged_myo_adata, mode='dynamical')
scv.tl.velocity_graph(merged_myo_adata)

scv.pl.proportions(merged_myo_adata, groupby='myo_annot', save='Merged_proportions.pdf')

scv.pl.velocity_embedding_grid(merged_myo_adata, basis='umap', title='Merged - Dynamical', scale=0.4,color='myo_annot',palette = ("#A3A500","#F8766D","#00BF7D","#00B0F6", "#E76BF3"),groups=("Atrial","AVC","OFT","SAN", "Ventricular"), save='merged_myo_dynamical.svg')

scv.pl.velocity_embedding_stream(merged_myo_adata,basis='umap', color='myo_annot', title='Merged - Dynamical', density=2.5,palette = ("#A3A500","#F8766D","#00BF7D","#00B0F6", "#E76BF3"),groups=("Atrial","AVC","OFT","SAN", "Ventricular"), save='merged_myo_dynamical_stream.svg', add_margin=0.05)

scv.tl.velocity_confidence(merged_myo_adata)
scv.pl.scatter(merged_myo_adata, c=('velocity_length', 'velocity_confidence'), cmap='coolwarm', perc=[5, 95],save= 'merged_myo_velocity_length_conf_dyn.svg')

scv.tl.latent_time(merged_myo_adata)
scv.pl.scatter(merged_myo_adata, color='latent_time', color_map='gnuplot', size=80,title='Merged - Latent Time', save= 'merged_myo_latent_time.svg')

################################################################################

adatas = [WT_myo_adata, KO_myo_adata]
myo_adata = ad.concat(adatas,  join="outer")

scv.tl.rank_velocity_genes(myo_adata, groupby='orig.ident', min_corr=.3)

df = scv.DataFrame(myo_adata.uns['rank_velocity_genes']['names'])
df.head()

df.to_csv(path_or_buf='myo_top_dynamical_genes.csv')

scv.pl.scatter(myo_adata, df['WT'][:10],nrows = 2, frameon=True, color='orig.ident', size=80, linewidth=1.5, figsize = (7,8),fontsize=20,legend_loc="right", save= 'WT_myo_top10_dynamical_1.pdf')

scv.pl.scatter(myo_adata, df['WT'][:10],nrows = 2, frameon=True, color='myo_annot', size=80, linewidth=1.5, figsize = (7,8),fontsize=20,legend_loc="right", save= 'WT_myo_top10_dynamical_2.pdf')

scv.pl.scatter(myo_adata, df['KO'][:10],nrows = 2, frameon=True, color='orig.ident', size=80, linewidth=1.5, figsize = (7,8),fontsize=20,legend_loc="right", save= 'KO_myo_top10_dynamical_1.pdf')

scv.pl.scatter(myo_adata, df['KO'][:10],nrows = 2, frameon=True, color='myo_annot', size=80, linewidth=1.5, figsize = (7,8),fontsize=20,legend_loc="right", save= 'KO_myo_top10_dynamical_2.pdf')
```

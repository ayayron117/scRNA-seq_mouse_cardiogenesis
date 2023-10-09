Quality Control & Cell-cycle Scoring
================
Aaron Mohammed

``` r
library(Seurat)
library(DoubletFinder)
library(SingleCellExperiment)
library(scDblFinder)
library(ggplot2)
library(qpdf)

proj_path <- dirname(getwd())
counts_dir <- file.path(proj_path, "filtered_feature_bc_matrices")

plots_dir <- file.path(getwd(), "Plots")
dir.create(plots_dir)
```

``` r
##################################### WT #######################################

# Import count matrix for the wildtype 
WT_mat <- Read10X(file.path(counts_dir, 
                            "WT", 
                            "filtered_feature_bc_matrix"),  
                  gene.column = 2,  
                  cell.column = 1,  
                  unique.features = TRUE,  
                  strip.suffix = FALSE)

# Create seurat object for the wildtype
WT_s <- CreateSeuratObject(WT_mat, 
                           project = "WT")

# Calculate percentages of mitochondrial and ribosomal genes
WT_s <- PercentageFeatureSet(WT_s, 
                             "^mt-", 
                             col.name = "pct_mito")
WT_s <- PercentageFeatureSet(WT_s, 
                             "^Rp[sl]\\d+", 
                             col.name = "pct_ribo")

# Number of cells
nrow(WT_s@meta.data) # 15953
```

``` r
##################################### KO #######################################

# Import count matrix for the knockout  
KO_mat <- Read10X(file.path(counts_dir, 
                            "KO", 
                            "filtered_feature_bc_matrix"),  
                  gene.column = 2,  
                  cell.column = 1,  
                  unique.features = TRUE,  
                  strip.suffix = FALSE)

# Create seurat object for the knockout
KO_s <- CreateSeuratObject(KO_mat, 
                           project = "KO")

# Calculate percentages of mitochondrial and ribosomal genes
KO_s <- PercentageFeatureSet(KO_s, 
                             "^mt-", 
                             col.name = "pct_mito")

KO_s <- PercentageFeatureSet(KO_s, 
                             "^Rp[sl]\\d+", 
                             col.name = "pct_ribo")

# Number of cells
nrow(KO_s@meta.data) # 8666
```

``` r
################################## WT Plots ####################################

# Ranges of counts, expressed genes, and percent mitochondria gene expression
range(WT_s$nCount_RNA) # min = 500; max = 101497
range(WT_s$nFeature_RNA) # min = 36; max = 10055
range(WT_s$pct_mito) # min = 0.04%; max = 96.42%

# Create directory for WT plots
WT_plots <- file.path(plots_dir, "WT")
dir.create(WT_plots)

# Scatter plot - nFeatures vs nCounts
pdf(file.path(WT_plots, "WT_unfiltered_scatter_1.pdf"), height = 8, width = 10)
FeatureScatter(WT_s,
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA", 
               group.by = "orig.ident") +
  ggtitle("WT - Unfiltered") + 
  theme(legend.position="none")
dev.off()

# Scatter plot - pct mito vs nCounts
pdf(file.path(WT_plots, "WT_unfiltered_scatter_2.pdf"), height = 8, width = 10)
FeatureScatter(WT_s,
               feature1 = "nFeature_RNA", 
               feature2 = "nCount_RNA", 
               group.by = "orig.ident") +
  ggtitle("WT - Unfiltered") + 
  theme(legend.position="none")
dev.off()

# Scatter plot - pct mito vs nCounts
pdf(file.path(WT_plots, "WT_unfiltered_scatter_3.pdf"), height = 8, width = 10)
FeatureScatter(WT_s,
               feature1 = "nCount_RNA", 
               feature2 = "pct_mito", 
               group.by = "orig.ident") +
  ggtitle("WT - Unfiltered") + 
  theme(legend.position="none")
dev.off()

# Scatter plot - pct mito vs nFeatures
pdf(file.path(WT_plots, "WT_unfiltered_scatter_4.pdf"), height = 8, width = 10)
FeatureScatter(WT_s,
               feature1 = "nFeature_RNA", 
               feature2 = "pct_mito", 
               group.by = "orig.ident") +
  ggtitle("WT - Unfiltered") + 
  theme(legend.position="none")
dev.off()

# Violin plot - nCounts - nFeature_RNA - pct_mito
pdf(file.path(WT_plots, "WT_unfiltered_vln_plots.pdf"), height = 6, width = 12)
VlnPlot(WT_s, 
        features = c("nCount_RNA","nFeature_RNA","pct_mito")) + 
  theme(legend.position="none")
dev.off()

pdf_combine(input = c(file.path(WT_plots, "WT_unfiltered_scatter_1.pdf"),
                      file.path(WT_plots, "WT_unfiltered_scatter_2.pdf"), 
                      file.path(WT_plots, "WT_unfiltered_scatter_3.pdf"),
                      file.path(WT_plots, "WT_unfiltered_scatter_4.pdf")),
              output = file.path(WT_plots, "WT_unfiltered_scatter_plots.pdf"))

file.remove(c(file.path(WT_plots, "WT_unfiltered_scatter_1.pdf"),
                      file.path(WT_plots, "WT_unfiltered_scatter_2.pdf"), 
                      file.path(WT_plots, "WT_unfiltered_scatter_3.pdf"),
                      file.path(WT_plots, "WT_unfiltered_scatter_4.pdf")))
```
<p align="center">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/c19487bb-4f94-4db2-98a3-98e20b7967c6">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/15a01277-9ea0-411e-ac39-b8fd8ea4cd98">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/71da9814-f350-4645-998d-ed3a3cd3ef73">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/fb68c4af-655f-4c56-b64c-923631d83ea3">
  <img width="800" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/544b8a07-80b4-4dbb-badb-e78b693081c8">
</p>

``` r
################################## KO Plots ####################################

# Ranges of counts, expressed genes, and percent mitochondria gene expression
range(KO_s$nCount_RNA) # min = 500; max = 94972
range(KO_s$nFeature_RNA) # min = 60; max = 9927
range(KO_s$pct_mito) # min = 0%; max = 97.5%

# Create directory for KO plots
KO_plots <- file.path(plots_dir, "KO")
dir.create(KO_plots)

# Scatter plot - nFeatures vs nCounts
pdf(file.path(KO_plots, "KO_unfiltered_scatter_1.pdf"), height = 8, width = 10)
FeatureScatter(KO_s,
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA", 
               group.by = "orig.ident") +
  ggtitle("KO - Unfiltered") + 
  theme(legend.position="none")
dev.off()

# Scatter plot - pct mito vs nCounts
pdf(file.path(KO_plots, "KO_unfiltered_scatter_2.pdf"), height = 8, width = 10)
FeatureScatter(KO_s,
               feature1 = "nFeature_RNA", 
               feature2 = "nCount_RNA", 
               group.by = "orig.ident") +
  ggtitle("KO - Unfiltered") + 
  theme(legend.position="none")
dev.off()

# Scatter plot - pct mito vs nCounts
pdf(file.path(KO_plots, "KO_unfiltered_scatter_3.pdf"), height = 8, width = 10)
FeatureScatter(KO_s,
               feature1 = "nCount_RNA", 
               feature2 = "pct_mito", 
               group.by = "orig.ident") +
  ggtitle("KO - Unfiltered") + 
  theme(legend.position="none")
dev.off()

# Scatter plot - pct mito vs nFeatures
pdf(file.path(KO_plots, "KO_unfiltered_scatter_4.pdf"), height = 8, width = 10)
FeatureScatter(KO_s,
               feature1 = "nFeature_RNA", 
               feature2 = "pct_mito", 
               group.by = "orig.ident") +
  ggtitle("KO - Unfiltered") + 
  theme(legend.position="none")
dev.off()

# Violin plot - nCounts - nFeature_RNA - pct_mito
pdf(file.path(KO_plots, "KO_unfiltered_vln_plots.pdf"), height = 6, width = 12)
VlnPlot(KO_s, 
        features = c("nCount_RNA","nFeature_RNA","pct_mito")) + 
  theme(legend.position="none")
dev.off()

pdf_combine(input = c(file.path(KO_plots, "KO_unfiltered_scatter_1.pdf"),
                      file.path(KO_plots, "KO_unfiltered_scatter_2.pdf"), 
                      file.path(KO_plots, "KO_unfiltered_scatter_3.pdf"),
                      file.path(KO_plots, "KO_unfiltered_scatter_4.pdf")),
              output = file.path(KO_plots, "KO_unfiltered_scatter_plots.pdf"))

file.remove(c(file.path(KO_plots, "KO_unfiltered_scatter_1.pdf"),
                      file.path(KO_plots, "KO_unfiltered_scatter_2.pdf"), 
                      file.path(KO_plots, "KO_unfiltered_scatter_3.pdf"),
                      file.path(KO_plots, "KO_unfiltered_scatter_4.pdf")))
```
<p align="center">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/3a70ed73-a065-48cc-8aa5-8da6a057dce1">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/e730578b-f5a9-47cb-ae42-83f295802627">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/660cd53d-b02e-4a4e-9f01-8d0dc1ba3210">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/cc4a422a-9e2e-45e4-bfdc-4c98c586dd06">
  <img width="800" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/245e9f61-68c5-4d43-ba44-46b56b98cc12">
</p>

``` r
################################ WT Filtering ##################################

WT_f <- subset(WT_s, 
               subset = nFeature_RNA <= 7000 &
                 nFeature_RNA >= 700 &
                 nCount_RNA >= 1000 &
                 nCount_RNA <= 40000 &
                 pct_mito <= 25)

# Number of cells
nrow(WT_f@meta.data) # 7213

# Ranges after filtering
range(WT_f$nCount_RNA) # min = 1000; max = 39639
range(WT_f$nFeature_RNA) # min = 700; max = 6981
range(WT_f$pct_mito) # min = 0.07%; max = 24.8%
```

``` r
################################ KO Filtering ##################################

KO_f <- subset(KO_s, 
               subset = nFeature_RNA <= 8000 &
                 nFeature_RNA >= 700 &
                 nCount_RNA >= 1000 &
                 nCount_RNA <= 50000 &
                 pct_mito <= 25)

# Number of cells
nrow(KO_f@meta.data) # 6493

# Ranges after filtering
range(KO_f$nCount_RNA) # min = 1020; max = 49056
range(KO_f$nFeature_RNA) # min = 700; max = 7966
range(KO_f$pct_mito) # min = 0%; max = 25%

```

``` r
########################### WT Cell Cycle Scoring ##############################

# Dimentionality reduction
WT_f <- NormalizeData(WT_f)
WT_f <- FindVariableFeatures(WT_f, 
                             selection.method = "vst", 
                             nfeatures = 2000)
WT_f <- ScaleData(WT_f)
WT_f <- RunPCA(WT_f)

s.genes <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2",
             "Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Hells","Rfc2","Rpa2","Nasp",
             "Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2",
             "Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm",
             "Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8")

g2m.genes <- c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80",
               "Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Pimreg",
               "Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e",
               "Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Jpt1","Cdc20","Ttk",
               "Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8",
               "Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5",
               "Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")

WT_f <- CellCycleScoring(WT_f, 
                         s.features = s.genes, 
                         g2m.features = g2m.genes, 
                         set.ident = FALSE)
```

``` r
############################## WT Doublet Removal ##############################

# Determine number of dimensions that capture at least 80% of the variation
pct <- WT_f[["pca"]]@stdev / sum(WT_f[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
PC_80_pct <- which(cumu >= 80)[1]
PC_80_pct #33

# UMAP and clustering
WT_f <- FindNeighbors(WT_f, 
                      dims = 1:PC_80_pct)
WT_f <- RunUMAP(WT_f, 
                dims=1:PC_80_pct)

WT_f <- FindClusters(WT_f, 
                     resolution = 1)

seurat_clusters <- WT_f@meta.data$seurat_clusters

# Doublet detection using DoubletFinder

set.seed(123)

# Parameter sweep
param_sweep_results <- paramSweep_v3(WT_f, 
                                     1:PC_80_pct, 
                                     sct = FALSE)
param_sweep_summarized <- summarizeSweep(param_sweep_results, 
                                         GT = FALSE)
bcmnv <- find.pK(param_sweep_summarized)

# Determine which pK value is associated with the highest value of BCmetric
pK <- as.numeric(as.character(bcmnv$pK[which(bcmnv$BCmetric == max(bcmnv$BCmetric))])) # 0.17

# Plot of BCmetric vs pK
WT_BC_pK <- ggplot(bcmnv, aes(pK, BCmetric, group = 1)) +
                   geom_point() +
                   geom_line()

pdf(file.path(WT_plots, "WT_BC_pK.pdf"), height = 8, width = 10)
WT_BC_pK
dev.off()

WT_BC_pK

# Adjust the expected number of doublets using estimated proportion of homotypic 
# doublets
nExp_poi <- round(0.0577*nrow(WT_f@meta.data))
homotypic_prop <- modelHomotypic(seurat_clusters)
nExp_poi_adj <- round(nExp_poi*(1-homotypic_prop))

# Run DoubletFinder
WT_f <- doubletFinder_v3(WT_f,
                         PCs = 1:PC_80_pct,
                         pN = 0.25,
                         pK = pK,
                         nExp = nExp_poi_adj,
                         reuse.pANN = FALSE,
                         sct = FALSE)

################################################################################

# Doublet detection using scDblFinder

WT_sce <- as.SingleCellExperiment(WT_f, 
                                   assay = "RNA")

set.seed(123)

WT_sce <- scDblFinder(WT_sce, 
                      clusters="seurat_clusters")

identical(WT_sce@colData@rownames, row.names(WT_f@meta.data)) # TRUE
identical(WT_sce$seurat_clusters, WT_f@meta.data[,"seurat_clusters"]) # TRUE

WT_f <- AddMetaData(WT_f, WT_sce$scDblFinder.class, col.name = "scDblFinder.class")
WT_f <- AddMetaData(WT_f, WT_sce$scDblFinder.score, col.name = "scDblFinder.score")
WT_f <- AddMetaData(WT_f, WT_sce$scDblFinder.weighted, col.name = "scDblFinder.weighted")
WT_f <- AddMetaData(WT_f, WT_sce$scDblFinder.difficulty, col.name = "scDblFinder.difficulty")
WT_f <- AddMetaData(WT_f, WT_sce$scDblFinder.cxds_score, col.name = "scDblFinder.cxds_score")
WT_f <- AddMetaData(WT_f, WT_sce$scDblFinder.mostLikelyOrigin, col.name = "scDblFinder.mostLikelyOrigin")
WT_f <- AddMetaData(WT_f, WT_sce$scDblFinder.originAmbiguous, col.name = "scDblFinder.originAmbiguous")

################################################################################

# Determine which cells were classified as doublets using both methods
WT_f <- AddMetaData(WT_f, rep(NA, nrow(WT_f@meta.data)), "dbl_both_methods")

TRUE -> WT_f$dbl_both_methods[which(WT_f$DF.classifications_0.25_0.17_392 == "Doublet" &
                                      WT_f$scDblFinder.class == "doublet")]

FALSE -> WT_f$dbl_both_methods[which(is.na(WT_f$dbl_both_methods))]

# Save filtered seurat with doublet classifications
saveRDS(WT_f, file.path(getwd(), "WT_before_doublet_removal.rds"))

# Filter out the cells classified as doublets using both doubletFinder and scDblFinder
Idents(WT_f) <- "dbl_both_methods"
WT_f <- subset(WT_f, idents = FALSE)

# Number of cells
nrow(WT_f@meta.data) # 7084

# Save final QC'd seurat
saveRDS(WT_f, file.path(getwd(), "WT_seurat.rds"))
```

``` r
############################## WT Filtered Plots ###############################

# Scatter plot - nFeatures vs nCounts
pdf(file.path(WT_plots, "WT_filtered_scatter_1.pdf"), height = 8, width = 10)
FeatureScatter(WT_f,
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA", 
               group.by = "orig.ident") +
  ggtitle("WT - Filtered") + 
  theme(legend.position="none")
dev.off()

# Scatter plot - pct mito vs nCounts
pdf(file.path(WT_plots, "WT_filtered_scatter_2.pdf"), height = 8, width = 10)
FeatureScatter(WT_f,
               feature1 = "nFeature_RNA", 
               feature2 = "nCount_RNA", 
               group.by = "orig.ident") +
  ggtitle("WT - Filtered") + 
  theme(legend.position="none")
dev.off()

# Scatter plot - pct mito vs nCounts
pdf(file.path(WT_plots, "WT_filtered_scatter_3.pdf"), height = 8, width = 10)
FeatureScatter(WT_f,
               feature1 = "nCount_RNA", 
               feature2 = "pct_mito", 
               group.by = "orig.ident") +
  ggtitle("WT - Filtered") + 
  theme(legend.position="none")
dev.off()

# Scatter plot - pct mito vs nFeatures
pdf(file.path(WT_plots, "WT_filtered_scatter_4.pdf"), height = 8, width = 10)
FeatureScatter(WT_f,
               feature1 = "nFeature_RNA", 
               feature2 = "pct_mito", 
               group.by = "orig.ident") +
  ggtitle("WT - Filtered") + 
  theme(legend.position="none")
dev.off()

# Violin plot - nCounts - nFeature_RNA - pct_mito
pdf(file.path(WT_plots, "WT_filtered_vln_plots.pdf"), height = 6, width = 12)
VlnPlot(WT_f, 
        features = c("nCount_RNA","nFeature_RNA","pct_mito")) + 
  theme(legend.position="none")
dev.off()

pdf_combine(input = c(file.path(WT_plots, "WT_filtered_scatter_1.pdf"),
                      file.path(WT_plots, "WT_filtered_scatter_2.pdf"), 
                      file.path(WT_plots, "WT_filtered_scatter_3.pdf"),
                      file.path(WT_plots, "WT_filtered_scatter_4.pdf")),
              output = file.path(WT_plots, "WT_filtered_scatter_plots.pdf"))

file.remove(c(file.path(WT_plots, "WT_filtered_scatter_1.pdf"),
                      file.path(WT_plots, "WT_filtered_scatter_2.pdf"), 
                      file.path(WT_plots, "WT_filtered_scatter_3.pdf"),
                      file.path(WT_plots, "WT_filtered_scatter_4.pdf")))
```
<p align="center">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/424dd031-8a6c-4e2e-bc90-7c63eff53ed3">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/1657cd84-f653-428d-bc52-53f2efd51630">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/845beaee-4fc3-401d-bf0e-8e12a603c7df">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/6e10a3fc-4255-480a-b4dd-5e801cfe07b7">
</p>

``` r
########################### KO Cell Cycle Scoring ##############################

# Dimentionality reduction
KO_f <- NormalizeData(KO_f)
KO_f <- FindVariableFeatures(KO_f, 
                             selection.method = "vst", 
                             nfeatures = 2000)
KO_f <- ScaleData(KO_f)
KO_f <- RunPCA(KO_f)

s.genes <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2",
             "Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Hells","Rfc2","Rpa2","Nasp",
             "Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2",
             "Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm",
             "Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8")

g2m.genes <- c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80",
               "Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Pimreg",
               "Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e",
               "Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Jpt1","Cdc20","Ttk",
               "Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8",
               "Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5",
               "Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")

KO_f <- CellCycleScoring(KO_f, 
                         s.features = s.genes, 
                         g2m.features = g2m.genes, 
                         set.ident = FALSE)
```

``` r
############################## KO Doublet Removal ##############################

# Determine number of dimensions that capture at least 80% of the variation
pct <- KO_f[["pca"]]@stdev / sum(KO_f[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
PC_80_pct <- which(cumu >= 80)[1]
PC_80_pct #31

# UMAP and clustering
KO_f <- FindNeighbors(KO_f, 
                      dims = 1:PC_80_pct)
KO_f <- RunUMAP(KO_f, 
                dims=1:PC_80_pct)
KO_f <- FindClusters(KO_f, 
                     resolution = 1)

seurat_clusters <- KO_f@meta.data$seurat_clusters

# Doublet detection using DoubletFinder

set.seed(123)

# Parameter sweep
param_sweep_results <- paramSweep_v3(KO_f, 
                                     1:PC_80_pct, 
                                     sct = FALSE)
param_sweep_summarized <- summarizeSweep(param_sweep_results, 
                                         GT = FALSE)
bcmnv <- find.pK(param_sweep_summarized)

# Determine which pK value is associated with the highest value of BCmetric
pK <- as.numeric(as.character(bcmnv$pK[which(bcmnv$BCmetric == max(bcmnv$BCmetric))])) # 0.18

# Plot of BCmetric vs pK
KO_BC_pK <- ggplot(bcmnv, aes(pK, BCmetric, group = 1)) +
                   geom_point() +
                   geom_line()

pdf(file.path(KO_plots, "KO_BC_pK.pdf"), height = 8, width = 10)
KO_BC_pK
dev.off()

KO_BC_pK

# Adjust the expected number of doublets using estimated proportion of homotypic 
# doublets
nExp_poi <- round(0.051944*nrow(KO_f@meta.data))
homotypic_prop <- modelHomotypic(seurat_clusters)
nExp_poi_adj <- round(nExp_poi*(1-homotypic_prop))

# Run DoubletFinder
KO_f <- doubletFinder_v3(KO_f,
                         PCs = 1:PC_80_pct,
                         pN = 0.25,
                         pK = pK,
                         nExp = nExp_poi_adj,
                         reuse.pANN = FALSE,
                         sct = FALSE)

################################################################################

# Doublet detection using scDblFinder

KO_sce <- as.SingleCellExperiment(KO_f, 
                                  assay = "RNA")

set.seed(123)

KO_sce <- scDblFinder(KO_sce, 
                       clusters="seurat_clusters")

identical(KO_sce@colData@rownames, row.names(KO_f@meta.data)) # TRUE
identical(KO_sce$seurat_clusters, KO_f@meta.data[,"seurat_clusters"]) # TRUE

KO_f <- AddMetaData(KO_f, KO_sce$scDblFinder.class, col.name = "scDblFinder.class")
KO_f <- AddMetaData(KO_f, KO_sce$scDblFinder.score, col.name = "scDblFinder.score")
KO_f <- AddMetaData(KO_f, KO_sce$scDblFinder.weighted, col.name = "scDblFinder.weighted")
KO_f <- AddMetaData(KO_f, KO_sce$scDblFinder.difficulty, col.name = "scDblFinder.difficulty")
KO_f <- AddMetaData(KO_f, KO_sce$scDblFinder.cxds_score, col.name = "scDblFinder.cxds_score")
KO_f <- AddMetaData(KO_f, KO_sce$scDblFinder.mostLikelyOrigin, col.name = "scDblFinder.mostLikelyOrigin")
KO_f <- AddMetaData(KO_f, KO_sce$scDblFinder.originAmbiguous, col.name = "scDblFinder.originAmbiguous")

################################################################################

# Determine which cells were classified as doublets using both methods
KO_f <- AddMetaData(KO_f, rep(NA, nrow(KO_f@meta.data)), "dbl_both_methods")

TRUE -> KO_f$dbl_both_methods[which(KO_f$DF.classifications_0.25_0.18_315 == "Doublet" &
                                      KO_f$scDblFinder.class == "doublet")]

FALSE -> KO_f$dbl_both_methods[which(is.na(KO_f$dbl_both_methods))]

# Save filtered seurat with doublet classifications
saveRDS(KO_f, file.path(getwd(), "KO_before_doublet_removal.rds"))

# Filter out the cells classified as doublets using both doubletFinder and scDblFinder
Idents(KO_f) <- "dbl_both_methods"
KO_f <- subset(KO_f, idents = FALSE)

# Number of cells
nrow(KO_f@meta.data) # 6376

# Save final QC'd seurat
saveRDS(KO_f, file.path(getwd(), "KO_seurat.rds"))
```

``` r
############################## KO Filtered Plots ###############################

# Scatter plot - nFeatures vs nCounts
pdf(file.path(KO_plots, "KO_filtered_scatter_1.pdf"), height = 8, width = 10)
FeatureScatter(KO_f,
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA", 
               group.by = "orig.ident") +
  ggtitle("KO - Filtered") + 
  theme(legend.position="none")
dev.off()

# Scatter plot - pct mito vs nCounts
pdf(file.path(KO_plots, "KO_filtered_scatter_2.pdf"), height = 8, width = 10)
FeatureScatter(KO_f,
               feature1 = "nFeature_RNA", 
               feature2 = "nCount_RNA", 
               group.by = "orig.ident") +
  ggtitle("KO - Filtered") + 
  theme(legend.position="none")
dev.off()

# Scatter plot - pct mito vs nCounts
pdf(file.path(KO_plots, "KO_filtered_scatter_3.pdf"), height = 8, width = 10)
FeatureScatter(KO_f,
               feature1 = "nCount_RNA", 
               feature2 = "pct_mito", 
               group.by = "orig.ident") +
  ggtitle("KO - Filtered") + 
  theme(legend.position="none")
dev.off()

# Scatter plot - pct mito vs nFeatures
pdf(file.path(KO_plots, "KO_filtered_scatter_4.pdf"), height = 8, width = 10)
FeatureScatter(KO_f,
               feature1 = "nFeature_RNA", 
               feature2 = "pct_mito", 
               group.by = "orig.ident") +
  ggtitle("KO - Filtered") + 
  theme(legend.position="none")
dev.off()

# Violin plot - nCounts - nFeature_RNA - pct_mito
pdf(file.path(KO_plots, "KO_filtered_vln_plots.pdf"), height = 6, width = 12)
VlnPlot(KO_f, 
        features = c("nCount_RNA","nFeature_RNA","pct_mito")) + 
  theme(legend.position="none")
dev.off()

pdf_combine(input = c(file.path(KO_plots, "KO_filtered_scatter_1.pdf"),
                      file.path(KO_plots, "KO_filtered_scatter_2.pdf"), 
                      file.path(KO_plots, "KO_filtered_scatter_3.pdf"),
                      file.path(KO_plots, "KO_filtered_scatter_4.pdf")),
              output = file.path(KO_plots, "KO_filtered_scatter_plots.pdf"))

file.remove(c(file.path(KO_plots, "KO_filtered_scatter_1.pdf"),
                      file.path(KO_plots, "KO_filtered_scatter_2.pdf"), 
                      file.path(KO_plots, "KO_filtered_scatter_3.pdf"),
                      file.path(KO_plots, "KO_filtered_scatter_4.pdf")))
```

``` r
writeLines(capture.output(sessionInfo()),  
           file.path(getwd(), "1_sessionInfo.txt"))
```

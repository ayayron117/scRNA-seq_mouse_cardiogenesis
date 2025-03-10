WT Quality Control and Cell Cycle Scoring
================
Aaron Mohammed

``` r
library(Seurat)
library(SingleCellExperiment)
library(DropletQC)
library(scDblFinder)
library(ggplot2)
library(qpdf)
library(openxlsx)

proj_path <- dirname(getwd())

WT_dir <- file.path(getwd(), "WT_QC")
dir.create(WT_dir)

WT_plots <- file.path(WT_dir, "WT_Plots")
dir.create(WT_plots)

WT_rds <- file.path(WT_dir, "WT_RDS_Files")
dir.create(WT_rds)

WT_markers <- file.path(WT_dir, "WT_Markers")
dir.create(WT_markers)
```

## Create WT Seurat Object

``` r
# Import count matrix for the wildtype 
WT_mat <- Read10X(file.path(proj_path, 
                            "WT",
                            "outs",
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

# Normalize 
WT_s <- NormalizeData(WT_s)

# Number of cells
nrow(WT_s@meta.data) # 15953
```

## WT DropletQC

``` r
WT_nf <- nuclear_fraction_tags(outs = file.path(proj_path,
                                                "WT",
                                                "outs"),
                               tiles = 1, cores = 9, verbose = FALSE)

saveRDS(WT_nf, file.path(WT_rds, "0_WT_nuclear_fraction.rds"))

identical(row.names(WT_s@meta.data), row.names(WT_nf)) # True

WT_s <- AddMetaData(WT_s, WT_nf$nuclear_fraction, col.name = "nuclear_fraction")

################################################################################

WT_s.nf.umi <- data.frame(nf=WT_s$nuclear_fraction,
                         umi=WT_s$nCount_RNA)

# Identify empty droplets
WT_s.ed <- identify_empty_drops(nf_umi=WT_s.nf.umi)

################################################################################

WT_s <- FindVariableFeatures(WT_s, 
                             selection.method = "vst", 
                             nfeatures = 2000)
WT_s <- ScaleData(WT_s)
WT_s <- RunPCA(WT_s)

# Determine number of dimensions that capture at least 80% of the variation
pct <- WT_s[["pca"]]@stdev / sum(WT_s[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
PC_80_pct <- which(cumu >= 80)[1]
PC_80_pct #32

# UMAP and clustering
WT_s <- FindNeighbors(WT_s, 
                      dims = 1:PC_80_pct)
WT_s <- RunUMAP(WT_s, 
                dims=1:PC_80_pct)

WT_s <- FindClusters(WT_s, 
                     resolution = 0.6)

WT_s <- AddMetaData(WT_s, WT_s$seurat_clusters, col.name = "WT_SC_1")

################################################################################

WT_s.ed$cell_type <- WT_s$WT_SC_1

# Identify damaged cells
WT_s.ed.dc <- identify_damaged_cells(WT_s.ed, verbose = FALSE, output_plots = FALSE)

WT_s <- AddMetaData(WT_s, WT_s.ed.dc$df$cell_status, col.name = "DropletQC")

# Save the seurat object
saveRDS(WT_s, file.path(WT_rds, "0_WT_seurat.rds"))
```

## Unfiltered WT Plots

``` r
# Ranges of counts, expressed genes, and percent mitochondria gene expression
range(WT_s$nCount_RNA) # min = 500; max = 101497
range(WT_s$nFeature_RNA) # min = 36; max = 10055
range(WT_s$pct_mito) # min = 0.04%; max = 96.42%
range(WT_s$pct_ribo) # min = 0.11%; max = 70.49%


pdf(file.path(WT_plots, "0_WT_unfiltered_scatter_plots.pdf"), height = 8, width = 12)
FeatureScatter(WT_s,
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA", 
               group.by = "DropletQC") +
  ggtitle("WT - Unfiltered")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_s,
               feature1 = "nCount_RNA", 
               feature2 = "pct_mito", 
               group.by = "DropletQC") +
  ggtitle("WT - Unfiltered")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_s,
               feature1 = "nCount_RNA", 
               feature2 = "pct_ribo", 
               group.by = "DropletQC") +
  ggtitle("WT - Unfiltered")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_s,
               feature1 = "nFeature_RNA", 
               feature2 = "nCount_RNA", 
               group.by = "DropletQC") +
  ggtitle("WT - Unfiltered")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_s,
               feature1 = "nFeature_RNA", 
               feature2 = "pct_mito", 
               group.by = "DropletQC") +
  ggtitle("WT - Unfiltered")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_s,
               feature1 = "nFeature_RNA", 
               feature2 = "pct_ribo", 
               group.by = "DropletQC") +
  ggtitle("WT - Unfiltered")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_s,
               feature1 = "nuclear_fraction", 
               feature2 = "nCount_RNA", 
               group.by = "DropletQC") +
  ggtitle("WT - Unfiltered")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_s,
               feature1 = "nuclear_fraction", 
               feature2 = "nFeature_RNA", 
               group.by = "DropletQC") +
  ggtitle("WT - Unfiltered")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_s,
               feature1 = "nuclear_fraction", 
               feature2 = "pct_mito", 
               group.by = "DropletQC") +
  ggtitle("WT - Unfiltered")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_s,
               feature1 = "nuclear_fraction", 
               feature2 = "pct_ribo", 
               group.by = "DropletQC") +
  ggtitle("WT - Unfiltered")
dev.off()
```

<img src="WT_QC/WT_QC_images/0_WT_unfiltered_scatter_plots/0001.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_scatter_plots/0002.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_scatter_plots/0003.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_scatter_plots/0004.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_scatter_plots/0005.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_scatter_plots/0006.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_scatter_plots/0007.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_scatter_plots/0008.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_scatter_plots/0009.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_scatter_plots/0010.png" width="70%" height="70%" />

``` r
# Violin plots
pdf(file.path(WT_plots, "WT_unfiltered_vln_plots_1.pdf"), height = 6, width = 8)
VlnPlot(WT_s, 
        features = "nCount_RNA", 
        group.by = "orig.ident") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # 
VlnPlot(WT_s, 
        features = "nFeature_RNA", 
        group.by = "orig.ident") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # 
VlnPlot(WT_s, 
        features = "pct_mito", 
        group.by = "orig.ident") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # 
VlnPlot(WT_s, 
        features = "pct_ribo", 
        group.by = "orig.ident") + 
  theme(legend.position="none")
dev.off()

# Violin plots grouped by DoubletQC classifications
pdf(file.path(WT_plots, "WT_unfiltered_vln_plots_2.pdf"), height = 6, width = 12)
VlnPlot(WT_s, 
        features = "nCount_RNA", 
        group.by = "DropletQC") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_s, 
        features = "nFeature_RNA", 
        group.by = "DropletQC") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_s, 
        features = "pct_mito", 
        group.by = "DropletQC") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_s, 
        features = "pct_ribo", 
        group.by = "DropletQC") + 
  theme(legend.position="none")
dev.off()


pdf(file.path(WT_plots, "WT_unfiltered_vln_plots_3.pdf"), height = 8, width = 15)
VlnPlot(WT_s, 
        features = "nCount_RNA", 
        group.by = "WT_SC_1") + 
  theme(legend.position="none")
# # # # # # # # # # # # # #
VlnPlot(WT_s, 
        features = "nFeature_RNA", 
        group.by = "WT_SC_1") + 
  theme(legend.position="none")
# # # # # # # # # # # # # #
VlnPlot(WT_s, 
        features = "pct_mito", 
        group.by = "WT_SC_1") + 
  theme(legend.position="none")
# # # # # # # # # # # # # #
VlnPlot(WT_s, 
        features = "pct_ribo", 
        group.by = "WT_SC_1") + 
  theme(legend.position="none")
# # # # # # # # # # # # # #
VlnPlot(WT_s, 
        features = "nuclear_fraction", 
        group.by = "WT_SC_1") + 
  theme(legend.position="none")
dev.off()


# Combine violin plot pdfs
pdf_combine(input = c(file.path(WT_plots, "WT_unfiltered_vln_plots_1.pdf"),
                      file.path(WT_plots, "WT_unfiltered_vln_plots_2.pdf"),
                      file.path(WT_plots, "WT_unfiltered_vln_plots_3.pdf")),
              output = file.path(WT_plots, "0_WT_unfiltered_vln_plots.pdf"))

file.remove(c(file.path(WT_plots, "WT_unfiltered_vln_plots_1.pdf"),
              file.path(WT_plots, "WT_unfiltered_vln_plots_2.pdf"),
              file.path(WT_plots, "WT_unfiltered_vln_plots_3.pdf")))
```

<img src="WT_QC/WT_QC_images/0_WT_unfiltered_vln_plots/0001.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_vln_plots/0002.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_vln_plots/0003.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_vln_plots/0004.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_vln_plots/0005.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_vln_plots/0006.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_vln_plots/0007.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_vln_plots/0008.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_vln_plots/0009.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_vln_plots/0010.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_vln_plots/0011.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_vln_plots/0012.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_vln_plots/0013.png" width="70%" height="70%" />

``` r
# UMAP plots
pdf(file.path(WT_plots, "0_WT_unfiltered_UMAPs.pdf"), height = 7, width = 12)
UMAPPlot(WT_s, group.by = "WT_SC_1", label =TRUE) + ggtitle("WT Unfiltered | Clusters")
UMAPPlot(WT_s, group.by = "DropletQC") + ggtitle("WT Unfiltered | DropletQC Classifications")
FeaturePlot(WT_s, features = "nCount_RNA") + ggtitle("WT Unfiltered | # of Counts")
FeaturePlot(WT_s, features = "nFeature_RNA") + ggtitle("WT Unfiltered | # of Genes")
FeaturePlot(WT_s, features = "pct_mito") + ggtitle("WT Unfiltered | % of Mito Genes")
FeaturePlot(WT_s, features = "pct_ribo") + ggtitle("WT Unfiltered | % of Ribo Genes")
FeaturePlot(WT_s, features = "nuclear_fraction") + ggtitle("WT Unfiltered | Nuclear Fraction")
dev.off()
```

<img src="WT_QC/WT_QC_images/0_WT_unfiltered_UMAPs/0001.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_UMAPs/0002.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_UMAPs/0003.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_UMAPs/0004.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_UMAPs/0005.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_UMAPs/0006.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/0_WT_unfiltered_UMAPs/0007.png" width="70%" height="70%" />

## Cluster Markers

``` r
get.all.markers <- function (seurat, ident, path, fname) {
  
  DefaultAssay(seurat) <- "RNA"
  seurat <- ScaleData(seurat, features = row.names(seurat))
  
  Idents(seurat) <- ident
  
  all_markers <- FindAllMarkers(seurat, logfc.threshold = 0.1)
  
  if (is.factor(seurat@meta.data[,ident])) {
    clusters <- levels(seurat@meta.data[,ident])
  } else if (is.character(seurat@meta.data[,ident])) {
    clusters <- sort(unique(seurat@meta.data[,ident]))
  }
  
  all_markers_list <- vector(mode="list", length = length(clusters))
  names(all_markers_list) <- clusters
  
  for (clust in clusters) {
    all_markers_list[[clust]] <- all_markers[which(all_markers$cluster == clust), 
                                             c("gene", "p_val", "avg_log2FC", 
                                               "pct.1", "pct.2", "p_val_adj")]
  }
  write.xlsx(all_markers_list, 
             file.path(path, fname), 
             rowNames = FALSE)
  
  return(all_markers_list)
  
}
```

``` r
all_markers <- get.all.markers(seurat = WT_s,
                               ident = "seurat_clusters",
                               path = WT_markers,
                               fname = "0_WT_cluster_markers.xlsx")
```

## WT Low-quality Filtering

``` r
WT_f <- subset(WT_s, 
               subset = nFeature_RNA <= 8000 &
                 nFeature_RNA >= 1500 & 
                 nCount_RNA >= 1000 &
                 nCount_RNA <= 40000 &
                 pct_mito <= 25)

Idents(WT_f) <- "DropletQC"
WT_f <- subset(WT_f, idents = "cell")

# Number of cells
nrow(WT_f@meta.data) # 5524
```

## WT Doublet Classifications

``` r
WT_f <- FindVariableFeatures(WT_f, 
                             selection.method = "vst", 
                             nfeatures = 2000)
WT_f <- ScaleData(WT_f)
WT_f <- RunPCA(WT_f)

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
                     resolution = 0.6)

WT_f <- AddMetaData(WT_f, WT_f$seurat_clusters, col.name = "WT_SC_2")

################################################################################

# Doublet detection using scDblFinder
WT_sce <- as.SingleCellExperiment(WT_f, 
                                   assay = "RNA")

set.seed(123)

WT_sce <- scDblFinder(WT_sce, 
                      clusters="WT_SC_2")

identical(WT_sce@colData@rownames, row.names(WT_f@meta.data)) # TRUE
identical(WT_sce$WT_SC_2, WT_f@meta.data[,"WT_SC_2"]) # TRUE

WT_f <- AddMetaData(WT_f, WT_sce$scDblFinder.class, col.name = "scDblFinder.class")
WT_f <- AddMetaData(WT_f, WT_sce$scDblFinder.score, col.name = "scDblFinder.score")

saveRDS(WT_f, file.path(WT_rds, "1_WT_seurat.rds"))
```

## WT Plots with scDblFinder Classifications

``` r
pdf(file.path(WT_plots, "1_WT_doublet_scatter_plots.pdf"), height = 8, width = 12)
FeatureScatter(WT_f,
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA", 
               group.by = "scDblFinder.class")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_f,
               feature1 = "nCount_RNA", 
               feature2 = "pct_mito", 
               group.by = "scDblFinder.class")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_f,
               feature1 = "nCount_RNA", 
               feature2 = "pct_ribo", 
               group.by = "scDblFinder.class")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_f,
               feature1 = "nFeature_RNA", 
               feature2 = "nCount_RNA", 
               group.by = "scDblFinder.class")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_f,
               feature1 = "nFeature_RNA", 
               feature2 = "pct_mito", 
               group.by = "scDblFinder.class")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_f,
               feature1 = "nFeature_RNA", 
               feature2 = "pct_ribo", 
               group.by = "scDblFinder.class")
dev.off()
```

<img src="WT_QC/WT_QC_images/1_WT_doublet_scatter_plots/0001.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_scatter_plots/0002.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_scatter_plots/0003.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_scatter_plots/0004.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_scatter_plots/0005.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_scatter_plots/0006.png" width="70%" height="70%" />

``` r
# Violin plots grouped by scDblFinder classifications
pdf(file.path(WT_plots, "1_WT_doublet_vln_plots.pdf"), height = 6, width = 12)
VlnPlot(WT_f, 
        features = "nCount_RNA", 
        group.by = "scDblFinder.class") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "nFeature_RNA", 
        group.by = "scDblFinder.class") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "pct_mito", 
        group.by = "scDblFinder.class") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "pct_ribo", 
        group.by = "scDblFinder.class") + 
  theme(legend.position="none")
dev.off()
```

<img src="WT_QC/WT_QC_images/1_WT_doublet_vln_plots/0001.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_vln_plots/0002.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_vln_plots/0003.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_vln_plots/0004.png" width="70%" height="70%" />

``` r
# Violin plots grouped by scDblFinder clusters
pdf(file.path(WT_plots, "1_WT_vln_scdblfinder_clusters.pdf"), height = 6, width = 12)
VlnPlot(WT_f, 
        features = "nCount_RNA", 
        group.by = "WT_SC_2") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "nFeature_RNA", 
        group.by = "WT_SC_2") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "pct_mito", 
        group.by = "WT_SC_2") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "pct_ribo", 
        group.by = "WT_SC_2") + 
  theme(legend.position="none")
dev.off()
```

<img src="WT_QC/WT_QC_images/1_WT_vln_scdblfinder_clusters/0001.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_vln_scdblfinder_clusters/0002.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_vln_scdblfinder_clusters/0003.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_vln_scdblfinder_clusters/0004.png" width="70%" height="70%" />

``` r
# Violin plots grouped by original clusters
pdf(file.path(WT_plots, "1_WT_vln_original_clusters.pdf"), height = 6, width = 12)
VlnPlot(WT_f, 
        features = "nCount_RNA", 
        group.by = "WT_SC_1") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "nFeature_RNA", 
        group.by = "WT_SC_1") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "pct_mito", 
        group.by = "WT_SC_1") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "pct_ribo", 
        group.by = "WT_SC_1") + 
  theme(legend.position="none")
dev.off()
```

<img src="WT_QC/WT_QC_images/1_WT_vln_original_clusters/0001.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_vln_original_clusters/0002.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_vln_original_clusters/0003.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_vln_original_clusters/0004.png" width="70%" height="70%" />

``` r
# UMAP plots
pdf(file.path(WT_plots, "1_WT_doublet_UMAPs.pdf"), height = 7, width = 12)
UMAPPlot(WT_f, group.by = "WT_SC_2", label =TRUE) + ggtitle("WT Doublets | Clusters")
UMAPPlot(WT_f, group.by = "scDblFinder.class") + ggtitle("WT Doublets | DropletQC Classifications")
FeaturePlot(WT_f, features = "nCount_RNA") + ggtitle("WT Doublets | # of Counts")
FeaturePlot(WT_f, features = "nFeature_RNA") + ggtitle("WT Doublets | # of Genes")
FeaturePlot(WT_f, features = "pct_mito") + ggtitle("WT Doublets | % of Mito Genes")
FeaturePlot(WT_f, features = "pct_ribo") + ggtitle("WT Doublets | % of Ribo Genes")
FeaturePlot(WT_f, features = "scDblFinder.score") + ggtitle("WT Doublets | scDblFinder Score")
dev.off()
```

<img src="WT_QC/WT_QC_images/1_WT_doublet_UMAPs/0001.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_UMAPs/0002.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_UMAPs/0003.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_UMAPs/0004.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_UMAPs/0005.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_UMAPs/0006.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/1_WT_doublet_UMAPs/0007.png" width="70%" height="70%" />

## WT Doublet Removal

``` r
# Filter out the cells classified as doublets
Idents(WT_f) <- "scDblFinder.class"
WT_f <- subset(WT_f, idents = "singlet")

# Number of cells
nrow(WT_f@meta.data) # 5055
```

## Dim Reduction & Clustering After QC

``` r
WT_f <- FindVariableFeatures(WT_f, 
                             selection.method = "vst", 
                             nfeatures = 2000)
WT_f <- ScaleData(WT_f)
WT_f <- RunPCA(WT_f)

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
                     resolution = 0.6)

WT_f <- AddMetaData(WT_f, WT_f$seurat_clusters, col.name = "WT_SC_3")

saveRDS(WT_f, file.path(WT_rds, "2_WT_seurat.rds"))
```

## WT Cell Cycle Scoring

``` r
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

saveRDS(WT_f, file.path(WT_rds, "3_WT_seurat.rds"))
```

## WT Filtered Plots

``` r
pdf(file.path(WT_plots, "2_WT_scatter_plots.pdf"), height = 8, width = 12)
FeatureScatter(WT_f,
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA", 
               group.by = "orig.ident")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_f,
               feature1 = "nCount_RNA", 
               feature2 = "pct_mito", 
               group.by = "orig.ident")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_f,
               feature1 = "nCount_RNA", 
               feature2 = "pct_ribo", 
               group.by = "orig.ident")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_f,
               feature1 = "nFeature_RNA", 
               feature2 = "nCount_RNA", 
               group.by = "orig.ident")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_f,
               feature1 = "nFeature_RNA", 
               feature2 = "pct_mito", 
               group.by = "orig.ident")
# # # # # # # # # # # # # # # # # # # # # # # #
FeatureScatter(WT_f,
               feature1 = "nFeature_RNA", 
               feature2 = "pct_ribo", 
               group.by = "orig.ident")
dev.off()
```

<img src="WT_QC/WT_QC_images/2_WT_scatter_plots/0001.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_scatter_plots/0002.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_scatter_plots/0003.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_scatter_plots/0004.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_scatter_plots/0005.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_scatter_plots/0006.png" width="70%" height="70%" />

``` r
# Violin plots grouped by scDblFinder classifications
pdf(file.path(WT_plots, "2_WT_vln_plots.pdf"), height = 6, width = 12)
VlnPlot(WT_f, 
        features = "nCount_RNA", 
        group.by = "orig.ident") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "nFeature_RNA", 
        group.by = "orig.ident") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "pct_mito", 
        group.by = "orig.ident") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "pct_ribo", 
        group.by = "orig.ident") + 
  theme(legend.position="none")
dev.off()
```

<img src="WT_QC/WT_QC_images/2_WT_vln_plots/0001.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_vln_plots/0002.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_vln_plots/0003.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_vln_plots/0004.png" width="70%" height="70%" />

``` r
# Violin plots grouped by final clusters
pdf(file.path(WT_plots, "2_WT_vln_final_clusters.pdf"), height = 6, width = 12)
VlnPlot(WT_f, 
        features = "nCount_RNA", 
        group.by = "WT_SC_3") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "nFeature_RNA", 
        group.by = "WT_SC_3") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "pct_mito", 
        group.by = "WT_SC_3") + 
  theme(legend.position="none")
# # # # # # # # # # # # # # # # # # # # # # # #
VlnPlot(WT_f, 
        features = "pct_ribo", 
        group.by = "WT_SC_3") + 
  theme(legend.position="none")
dev.off()
```

<img src="WT_QC/WT_QC_images/2_WT_vln_final_clusters/0001.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_vln_final_clusters/0002.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_vln_final_clusters/0003.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_vln_final_clusters/0004.png" width="70%" height="70%" />

``` r
# UMAP plots
pdf(file.path(WT_plots, "2_WT_UMAPs.pdf"), height = 7, width = 12)
UMAPPlot(WT_f, group.by = "WT_SC_3", label =TRUE) + ggtitle("WT | Clusters")
UMAPPlot(WT_f, group.by = "Phase", label =FALSE) + ggtitle("WT | Cell Cycle")
FeaturePlot(WT_f, features = "nCount_RNA") + ggtitle("WT | # of Counts")
FeaturePlot(WT_f, features = "nFeature_RNA") + ggtitle("WT | # of Genes")
FeaturePlot(WT_f, features = "pct_mito") + ggtitle("WT | % of Mito Genes")
FeaturePlot(WT_f, features = "pct_ribo") + ggtitle("WT | % of Ribo Genes")
dev.off()
```

<img src="WT_QC/WT_QC_images/2_WT_UMAPs/0001.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_UMAPs/0002.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_UMAPs/0003.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_UMAPs/0004.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_UMAPs/0005.png" width="70%" height="70%" /><img src="WT_QC/WT_QC_images/2_WT_UMAPs/0006.png" width="70%" height="70%" />

``` r
# Clean up the seurat (it'll later be merged with the KO seurat, this is just my way of making the 
# final merged seurat organized)
WT_f$RNA_snn_res.0.6 <- NULL 
WT_f$seurat_clusters <- NULL 

saveRDS(WT_f, file.path(getwd(), "WT_seurat.rds"))
```

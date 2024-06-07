Integration & Final Annotations
================
Aaron Mohammed

``` r
library(Seurat)
library(ggplot2)
library(qpdf)
library(openxlsx)
library(plotly)

merged_dir <- file.path(getwd(), "Merged")
dir.create(merged_dir)

plots_dir <- file.path(merged_dir, "Plots")
dir.create(plots_dir)

markers_dir <- file.path(merged_dir, "Markers")
dir.create(markers_dir)

RDS_dir <- file.path(merged_dir, "RDS_files")
dir.create(RDS_dir)
```

``` r
WT_s <- readRDS(file.path(getwd(), "WT", "RDS_files", "2_WT_seurat.rds"))
KO_s <- readRDS(file.path(getwd(), "KO", "RDS_files", "2_KO_seurat.rds"))
```

``` r
# Merge wildtype and knockout seurats
merged_s <- merge(x = WT_s, 
                  y = KO_s,
                  add.cell.ids = c("WT",
                                   "KO"),
                  project = "scRNAseq_Ren")

# Clean up the merged seurat object
strings <- rep("RNA_snn_res.", length = length(res))
columns <- paste(strings, seq(0.1,1.5,0.1), sep = "")

merged_s@meta.data[ , columns] <- NULL
merged_s$seurat_clusters <- NULL

colnames(merged_s@meta.data)

merged_s@meta.data <- merged_s@meta.data[ , c("orig.ident", "nCount_RNA", "nFeature_RNA", 
                                              "pct_mito", "pct_ribo", "nuclear_fraction", 
                                              "WT_SC_1", "KO_SC_1", "DropletQC", 
                                              "WT_SC_2", "KO_SC_2", "scDblFinder.class", "scDblFinder.score", 
                                              "WT_SC_3", "KO_SC_3", "S.Score", "G2M.Score", "Phase", 
                                              "WT_clusters", "merged_SC_1", 
                                              "broad_annot_1", "broad_annot_2", 
                                              "specific_annot_1")]
```

``` r
integrate.and.umap <- function (seurat, nfeat = 2000, cc_regress = FALSE, 
                                pca_npcs = 50, cumm_var = NULL, dims = NULL,
                                verb = FALSE) {
  
  if (is.null(cumm_var) & is.null(dims)) {
    stop("Number of PCs to use has not been set, use either 'cumm_var' or 'dims' to set it")
  } else if (!is.null(cumm_var) & !is.null(dims)) {
    stop("Set the number of PCs to use with either 'cumm_var' or 'dims', not both")
  } 
  
  DefaultAssay(seurat) <- "RNA"

  seurat.list <- SplitObject(seurat, 
                             split.by = "orig.ident")
  
  seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    
    x <- NormalizeData(x, 
                       verbose = verb)
    x <- FindVariableFeatures(x, 
                              selection.method = "vst", 
                              nfeatures = nfeat,
                              verbose = verb)
    
  })
  
  features <- SelectIntegrationFeatures(object.list = seurat.list,
                                        nfeatures = nfeat,
                                        verbose = verb)
  
  anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                    anchor.features = features,
                                    verbose = verb)
  
  seurat <- IntegrateData(anchorset = anchors,
                          verbose = verb)
  
  DefaultAssay(seurat) <- "integrated"
  
  if (cc_regress == TRUE) {
    seurat <- ScaleData(seurat, 
                        vars.to.regress = c("S.Score", "G2M.Score"),
                        verbose = verb)
  } else if (cc_regress == FALSE) {
    seurat <- ScaleData(seurat,
                        verbose = verb)
  }
  
  # Principal component analysis
  seurat <- RunPCA(seurat, 
                   npcs = pca_npcs,
                   verbose = verb)
  
  if (!is.null(cumm_var)) {
    
  # Percent variation associated with each PC
  pct <- seurat[["pca"]]@stdev / sum(seurat[["pca"]]@stdev) * 100
  
  # Cumulative variation with increasing number of PCs
  cumu <- cumsum(pct)
  
  PC <- which(cumu >= cumm_var)[1] 
  
  seurat <- FindNeighbors(seurat, 
                          dims = 1:PC,
                          verbose = verb)
  seurat <- RunUMAP(seurat,
                    dims = 1:PC,
                    verbose = verb)
  
  } else if (!is.null(dims)) {
    
    seurat <- FindNeighbors(seurat, 
                            dims = 1:dims,
                            verbose = verb)
    seurat <- RunUMAP(seurat,
                      dims = 1:dims,
                      verbose = verb)
    
  }
  
  return(seurat)

}
```

``` r
sweep.cluster.res <- function (seurat, res, path, fname) {
  
  i = 1
  file_paths <- vector(length = length(res))
  
  dims <- length(seurat@commands[["RunUMAP.integrated.pca"]]@params[["dims"]])
  
  for (r in res) {
    
    seurat <- FindClusters(seurat, 
                           resolution = r,
                           verbose = FALSE)
    
    title <- paste0("Seurat Clusters", sep= " - ",
                    dims, sep= " ",
                    "Dimensions", sep= " - ",
                    r, sep= " ",
                    "Resolution")
    
    file_name <- paste0("cluster_umap", sep= "_",
                        r, sep= "_",
                        "res.pdf")
    
    file_paths[i] <- file.path(path, file_name)
    
    pdf(file_paths[i], height = 7, width = 12)
    print(UMAPPlot(seurat, 
                   group.by = "seurat_clusters", 
                   label = TRUE) +
            ggtitle(title))
    dev.off()
    
    i = i + 1
    
  }
  
  pdf_combine(input = file_paths, 
              output = file.path(path, fname))
  
  file.remove(file_paths)
  
}
```

``` r
get.cons.markers <- function (seurat, ident, group, path, fname) {
  
  DefaultAssay(seurat) <- "RNA"
  seurat <- ScaleData(seurat, features = row.names(seurat))
  
  Idents(seurat) <- ident
  
  if (is.factor(seurat@meta.data[,ident])) {
    clusters <- levels(seurat@meta.data[,ident])
  } else if (is.character(seurat@meta.data[,ident])) {
    clusters <- sort(unique(seurat@meta.data[,ident]))
  }
  
  conserved_markers <- vector(mode = "list", length = length(clusters))
  names(conserved_markers) <- clusters
  
  for (clust in clusters) {
    conserved_markers[[clust]] <- FindConservedMarkers(seurat, 
                                                       ident.1 = clust,
                                                       grouping.var = group)
  }
  
  write.xlsx(conserved_markers, 
             file.path(path, fname), 
             rowNames = TRUE)
  
  return(conserved_markers)
  
}
```

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
merged_s <- integrate.and.umap(seurat = merged_s,
                               nfeat = 2000,
                               cc_regress = TRUE,
                               pca_npcs = 50,
                               cumm_var = 80,
                               verb = FALSE)
```

``` r
res <- seq(0.1,1.5,0.1)

sweep.cluster.res(seurat = merged_s,
                  res = res,
                  path = plots_dir,
                  fname = "1_cluster_res_sweep.pdf")
```

``` r
merged_s <- FindClusters(merged_s, 
                         resolution = 1,
                         verbose = FALSE)

merged_s$merged_SC_1 <- merged_s$seurat_clusters

saveRDS(merged_s, file.path(RDS_dir, "1_merged_seurat.rds"))

pdf(file.path(plots_dir, "2_clusters.pdf"), height = 7, width = 12)
UMAPPlot(merged_s, group.by = "merged_SC_1")
UMAPPlot(merged_s, group.by = "broad_annot_2")
UMAPPlot(merged_s, group.by = "orig.ident")
UMAPPlot(merged_s, group.by = "specific_annot_1", label = TRUE)
UMAPPlot(merged_s, group.by = "Phase")
dev.off()
```

``` r
merged_s <- integrate.and.umap(seurat = merged_s,
                               nfeat = 2000,
                               cc_regress = FALSE,
                               pca_npcs = 50,
                               cumm_var = 80,
                               verb = FALSE)

pdf(file.path(plots_dir, "2_clusters_no_cc_reg.pdf"), height = 7, width = 12)
UMAPPlot(merged_s, group.by = "merged_SC_1")
UMAPPlot(merged_s, group.by = "broad_annot_2")
UMAPPlot(merged_s, group.by = "orig.ident")
UMAPPlot(merged_s, group.by = "specific_annot_1", label = TRUE)
UMAPPlot(merged_s, group.by = "Phase")
dev.off()

merged_s <- readRDS(file.path(RDS_dir, "1_merged_seurat.rds"))
```

``` r
all_markers <- get.all.markers(seurat = merged_s,
                               ident = "merged_SC_1",
                               path = markers_dir,
                               fname = "1_all_markers.xlsx")

saveRDS(all_markers, file.path(markers_dir, "1_all_markers.rds"))

################################################################################

cons_markers <- get.cons.markers(seurat = merged_s,
                                 ident = "merged_SC_1",
                                 group = "orig.ident",
                                 path = markers_dir,
                                 fname = "1_conserved_markers.xlsx")

saveRDS(cons_markers, file.path(markers_dir, "1_conserved_markers.rds"))
```

``` r
"Myocardial" -> merged_s@meta.data[which(merged_s$merged_SC_1 == 2 |
                                         merged_s$merged_SC_1 == 16), 
                                   "broad_annot_3"]

"Epicardial" -> merged_s@meta.data[which(merged_s$merged_SC_1 == 1 |
                                         merged_s$merged_SC_1 == 4 |
                                         merged_s$merged_SC_1 == 5 |
                                         merged_s$merged_SC_1 == 6 |
                                         merged_s$merged_SC_1 == 9 |
                                         merged_s$merged_SC_1 == 12 |
                                         merged_s$merged_SC_1 == 13 |
                                         merged_s$merged_SC_1 == 23 |
                                         merged_s$merged_SC_1 == 25), 
                                   "broad_annot_3"]


"Endocardial" -> merged_s@meta.data[which(merged_s$merged_SC_1 == 0 |
                                          merged_s$merged_SC_1 == 7 |
                                          merged_s$merged_SC_1 == 18 |
                                          merged_s$merged_SC_1 == 22), 
                                    "broad_annot_3"]

"Mesenchymal" -> merged_s@meta.data[which(merged_s$merged_SC_1 == 3), 
                                    "broad_annot_3"]

"MP" -> merged_s@meta.data[which(merged_s$merged_SC_1 == 10 |
                                 merged_s$merged_SC_1 == 11 |
                                 merged_s$merged_SC_1 == 14 |
                                 merged_s$merged_SC_1 == 20 |
                                 merged_s$merged_SC_1 == 26 |
                                 merged_s$merged_SC_1 == 27), 
                           "broad_annot_3"]

"Blood" -> merged_s@meta.data[which(merged_s$merged_SC_1 == 8 |
                                    merged_s$merged_SC_1 == 10 |
                                    merged_s$merged_SC_1 == 15 |
                                    merged_s$merged_SC_1 == 17 |
                                    merged_s$merged_SC_1 == 19 |
                                    merged_s$merged_SC_1 == 21), 
                              "broad_annot_3"]

"Hepatocyte" -> merged_s@meta.data[which(merged_s$merged_SC_1 == 24), 
                                   "broad_annot_3"]

pdf(file.path(plots_dir, "3_broad_annot.pdf"), height = 7, width = 12)
UMAPPlot(merged_s, group.by = "broad_annot_3")
UMAPPlot(merged_s, group.by = "broad_annot_3", label = TRUE)
dev.off()
```

``` r
merged_umap <- UMAPPlot(merged_s, group.by = "broad_annot_3", label = TRUE)

myo_cells <- CellSelector(merged_umap)

epi_cells <- CellSelector(merged_umap)
epi_cells <- c(epi_cells, CellSelector(merged_umap))

endo_cells <- CellSelector(merged_umap)
endo_cells <- c(endo_cells, CellSelector(merged_umap))

mes_cells <- CellSelector(merged_umap)

MP_cells <- CellSelector(merged_umap)
MP_cells <- c(MP_cells, CellSelector(merged_umap))

blo_cells <- CellSelector(merged_umap)
blo_cells <- c(blo_cells, CellSelector(merged_umap))

hep_cells <- CellSelector(merged_umap)


merged_s <- AddMetaData(merged_s, merged_s$broad_annot_3, col.name = "Broad_Annotation")

"Myocardial" -> merged_s@meta.data[myo_cells, "Broad_Annotation"]
"Epicardial" -> merged_s@meta.data[epi_cells, "Broad_Annotation"]
"Endocardial" -> merged_s@meta.data[endo_cells, "Broad_Annotation"]
"Mesenchymal" -> merged_s@meta.data[mes_cells, "Broad_Annotation"]
"MP" -> merged_s@meta.data[MP_cells, "Broad_Annotation"]
"Blood" -> merged_s@meta.data[blo_cells, "Broad_Annotation"]
"Hepatocyte" -> merged_s@meta.data[hep_cells, "Broad_Annotation"]

pdf(file.path(plots_dir, "3_Broad_Annotation.pdf"), height = 7, width = 12)
UMAPPlot(merged_s, group.by = "Broad_Annotation")
UMAPPlot(merged_s, group.by = "orig.ident")
UMAPPlot(merged_s, group.by = "Phase")
dev.off()

saveRDS(merged_s, file.path(RDS_dir, "2_merged_seurat.rds"))
```

``` r
Idents(merged_s) <- "Broad_Annotation"

merged_s <- subset(merged_s, 
                   idents= c("Myocardial", "Epicardial", "Endocardial", "MP",
                             "Mesenchymal"))

merged_s <- integrate.and.umap(seurat = merged_s,
                               nfeat = 2000,
                               cc_regress = TRUE,
                               pca_npcs = 50,
                               cumm_var = 70,
                               verb = FALSE)

pdf(file.path(plots_dir, "4_cardiac.pdf"), height = 7, width = 12)
UMAPPlot(merged_s, group.by = "Broad_Annotation", label=FALSE)
UMAPPlot(merged_s, group.by = "Broad_Annotation", label=TRUE)
UMAPPlot(merged_s, group.by = "specific_annot_1", label=FALSE)
UMAPPlot(merged_s, group.by = "specific_annot_1", label=TRUE)
UMAPPlot(merged_s, group.by = "orig.ident")
UMAPPlot(merged_s, group.by = "Phase")
dev.off()

saveRDS(merged_s, file.path(RDS_dir, "3_merged_seurat.rds"))
```

``` r
merged_s <- integrate.and.umap(seurat = merged_s,
                               nfeat = 2000,
                               cc_regress = FALSE,
                               pca_npcs = 50,
                               cumm_var = 70,
                               verb = FALSE)

pdf(file.path(plots_dir, "4_cardiac_no_cc_reg.pdf"), height = 7, width = 12)
UMAPPlot(merged_s, group.by = "Broad_Annotation", label=FALSE)
UMAPPlot(merged_s, group.by = "Broad_Annotation", label=TRUE)
UMAPPlot(merged_s, group.by = "specific_annot_1", label=FALSE)
UMAPPlot(merged_s, group.by = "specific_annot_1", label=TRUE)
UMAPPlot(merged_s, group.by = "orig.ident")
UMAPPlot(merged_s, group.by = "Phase")
dev.off()

merged_s <- readRDS(file.path(RDS_dir, "3_merged_seurat.rds"))
```

``` r
DefaultAssay(merged_s) <- "RNA"
pdf(file.path(plots_dir, "5_CPC_markers.pdf"), height = 7, width = 12)
FeaturePlot(merged_s, features = "Gata4")
FeaturePlot(merged_s, features = "Nkx2-5")
FeaturePlot(merged_s, features = "Isl1")
FeaturePlot(merged_s, features = "Tbx5")
FeaturePlot(merged_s, features = "Tbx20")
FeaturePlot(merged_s, features = "Mef2c")
FeaturePlot(merged_s, features = "Hand1")
FeaturePlot(merged_s, features = "Hand2")
dev.off()
```

``` r
Idents(merged_s) <- "Broad_Annotation"
  
myo_s <- subset(merged_s, idents = "Myocardial")

myo_s <- integrate.and.umap(seurat = myo_s,
                             nfeat = 2000,
                             cc_regress = TRUE,
                             pca_npcs = 50,
                             cumm_var = 80,
                             verb = FALSE)

sweep.cluster.res(seurat = myo_s,
                  res = res,
                  path = plots_dir,
                  fname = "6_myo_res_sweep.pdf")
```

``` r
myo_s <- FindClusters(myo_s, 
                      resolution = 0.9,
                      verbose = FALSE)

myo_s$myo_SC_1 <- myo_s$seurat_clusters

# saveRDS(myo_s, file.path(RDS_dir, "4_myo_seurat.rds"))

pdf(file.path(plots_dir, "7_myo_subset.pdf"), height = 7, width = 12)
UMAPPlot(myo_s, group.by = "myo_SC_1")
UMAPPlot(myo_s, group.by = "specific_annot_1")
UMAPPlot(myo_s, group.by = "specific_annot_1", label = TRUE)
UMAPPlot(myo_s, group.by = "orig.ident")
dev.off()
```

``` r
all_markers <- get.all.markers(seurat = myo_s,
                               ident = "myo_SC_1",
                               path = markers_dir,
                               fname = "2_myo_all_markers.xlsx")

saveRDS(all_markers, file.path(markers_dir, "2_myo_all_markers.rds"))

################################################################################

cons_markers <- get.cons.markers(seurat = myo_s,
                                 ident = "myo_SC_1",
                                 group = "orig.ident",
                                 path = markers_dir,
                                 fname = "2_myo_conserved_markers.xlsx")

saveRDS(cons_markers, file.path(markers_dir, "2_myo_conserved_markers.rds"))
```

``` r
"Ven_1" -> myo_s@meta.data[which(myo_s$myo_SC_1 == 0 |
                                        myo_s$myo_SC_1 == 1 |
                                        myo_s$myo_SC_1 == 2), 
                                  "myo_annot"]

"Ven_2" -> myo_s@meta.data[which(myo_s$myo_SC_1 == 6), "myo_annot"]

"SAN" -> myo_s@meta.data[which(myo_s$myo_SC_1 == 7), "myo_annot"]

"OFT" -> myo_s@meta.data[which(myo_s$myo_SC_1 == 5), "myo_annot"]

"AVC" -> myo_s@meta.data[which(myo_s$myo_SC_1 == 4), "myo_annot"]

"Atrial" -> myo_s@meta.data[which(myo_s$myo_SC_1 == 3), "myo_annot"]

pdf(file.path(plots_dir, "8_myo_annotated.pdf"), height = 7, width = 12)
UMAPPlot(myo_s, group.by = "myo_SC_1", label = TRUE)
UMAPPlot(myo_s, group.by = "myo_annot")
UMAPPlot(myo_s, group.by = "myo_annot", label = TRUE)
UMAPPlot(myo_s, group.by = "orig.ident")
UMAPPlot(myo_s, group.by = "Phase")
dev.off()

saveRDS(myo_s, file.path(RDS_dir, "5_myo_seurat.rds"))

# myo_s <- readRDS(file.path(RDS_dir, "5_myo_seurat.rds"))
```

``` r
merged_s <- AddMetaData(merged_s, merged_s$specific_annot_1, col.name = "Specific_Annotation")

for (cell in row.names(myo_s@meta.data)) {
  myo_s@meta.data[cell, "myo_annot"] -> merged_s@meta.data[which(row.names(merged_s@meta.data) == cell), "Specific_Annotation"]
}

umap <- UMAPPlot(merged_s, group.by = "Specific_Annotation", label=TRUE)

cells <- CellSelector(umap)
"Epithelial" -> merged_s@meta.data[cells, "Specific_Annotation"]

cells <- CellSelector(umap)
"EpiMT" -> merged_s@meta.data[cells, "Specific_Annotation"]

cells <- CellSelector(umap)
"Mesenchymal" -> merged_s@meta.data[cells, "Specific_Annotation"]

cells <- CellSelector(umap)
"EndMT" -> merged_s@meta.data[cells, "Specific_Annotation"]

"EndMT" -> merged_s$Specific_Annotation[which(merged_s$Specific_Annotation == "Blood")]
"MP" -> merged_s$Specific_Annotation[which(merged_s$Specific_Annotation == "SAN_2")]
"SAN" -> merged_s$Specific_Annotation[which(merged_s$Specific_Annotation == "SAN_1")]

pdf(file.path(plots_dir, "9_cardiac_annotated.pdf"), height = 7, width = 12)
UMAPPlot(merged_s, group.by = "Specific_Annotation")
UMAPPlot(merged_s, group.by = "Specific_Annotation", label=TRUE)
UMAPPlot(merged_s, group.by = "orig.ident")
UMAPPlot(merged_s, group.by = "Phase")
dev.off()

saveRDS(merged_s, file.path(getwd(), "final_cardiac_seurat.rds"))
```

``` r
full_s <- readRDS(file.path(RDS_dir, "2_merged_seurat.rds"))

full_s <- AddMetaData(full_s, full_s$specific_annot_1, col.name = "Specific_Annotation")

for (cell in row.names(merged_s@meta.data)) {
  merged_s@meta.data[cell, "Specific_Annotation"] -> full_s@meta.data[which(row.names(full_s@meta.data) == cell), "Specific_Annotation"]
}

pdf(file.path(plots_dir, "10_specific_annot.pdf"), height = 7, width = 12)
UMAPPlot(full_s, group.by = "Specific_Annotation")
UMAPPlot(full_s, group.by = "Specific_Annotation", label=TRUE)
UMAPPlot(full_s, group.by = "orig.ident")
UMAPPlot(full_s, group.by = "Phase")
dev.off()

saveRDS(full_s, file.path(getwd(), "final_seurat.rds"))
```

``` r
pdf(file.path(plots_dir, "11_broad_annot.pdf"), height = 7, width = 12)
UMAPPlot(full_s, group.by = "Broad_Annotation")
UMAPPlot(full_s, group.by = "Broad_Annotation", label=TRUE)
UMAPPlot(full_s, group.by = "orig.ident")
UMAPPlot(full_s, group.by = "Phase")
dev.off()
```

``` r
pdf(file.path(plots_dir, "12_All_clusters_cc_regressed.pdf"), height = 7, width = 12)
UMAPPlot(full_s, group.by = "Broad_Annotation")
UMAPPlot(full_s, group.by = "Broad_Annotation", label=TRUE)
UMAPPlot(full_s, group.by = "Specific_Annotation")
UMAPPlot(full_s, group.by = "Specific_Annotation", label=TRUE)
UMAPPlot(full_s, group.by = "orig.ident")
UMAPPlot(full_s, group.by = "Phase") + ggtitle("Phase - Cell Cycle Regressed")
FeaturePlot(full_s, features = "nCount_RNA")
FeaturePlot(full_s, features = "nFeature_RNA")
FeaturePlot(full_s, features = "pct_mito")
FeaturePlot(full_s, features = "pct_ribo")
dev.off()
```

``` r
cardiac_s <- readRDS(file.path(getwd(), "final_cardiac_seurat.rds"))

pdf(file.path(plots_dir, "13_Cardiac_cc_regressed.pdf"), height = 7, width = 12)
UMAPPlot(cardiac_s, group.by = "Broad_Annotation")
UMAPPlot(cardiac_s, group.by = "Broad_Annotation", label=TRUE)
UMAPPlot(cardiac_s, group.by = "Specific_Annotation")
UMAPPlot(cardiac_s, group.by = "Specific_Annotation", label=TRUE)
UMAPPlot(cardiac_s, group.by = "orig.ident")
UMAPPlot(cardiac_s, group.by = "Phase") + ggtitle("Phase - Cell Cycle Regressed")
FeaturePlot(cardiac_s, features = "nCount_RNA")
FeaturePlot(cardiac_s, features = "nFeature_RNA")
FeaturePlot(cardiac_s, features = "pct_mito")
FeaturePlot(cardiac_s, features = "pct_ribo")
dev.off()
```

``` r
myo_s <- readRDS(file.path(getwd(), "final_myo_seurat.rds"))

pdf(file.path(plots_dir, "14_Myocardial_cc_regressed.pdf"), height = 7, width = 12)
UMAPPlot(myo_s, group.by = "myo_annot")
UMAPPlot(myo_s, group.by = "myo_annot", label=TRUE)
UMAPPlot(myo_s, group.by = "orig.ident")
UMAPPlot(myo_s, group.by = "Phase") + ggtitle("Phase - Cell Cycle Regressed")
FeaturePlot(myo_s, features = "nCount_RNA")
FeaturePlot(myo_s, features = "nFeature_RNA")
FeaturePlot(myo_s, features = "pct_mito")
FeaturePlot(myo_s, features = "pct_ribo")
dev.off()
```

``` r
full_s <- readRDS(file.path(getwd(), "final_seurat.rds"))

full_s <- integrate.and.umap(seurat = full_s,
                               nfeat = 2000,
                               cc_regress = FALSE,
                               pca_npcs = 50,
                               cumm_var = 80,
                               verb = FALSE)

pdf(file.path(plots_dir, "15_All_clusters_no_cc_reg.pdf"), height = 7, width = 12)
UMAPPlot(full_s, group.by = "Broad_Annotation")
UMAPPlot(full_s, group.by = "Broad_Annotation", label=TRUE)
UMAPPlot(full_s, group.by = "Specific_Annotation")
UMAPPlot(full_s, group.by = "Specific_Annotation", label=TRUE)
UMAPPlot(full_s, group.by = "orig.ident")
UMAPPlot(full_s, group.by = "Phase") + ggtitle("Phase - No CC Regression")
FeaturePlot(full_s, features = "nCount_RNA")
FeaturePlot(full_s, features = "nFeature_RNA")
FeaturePlot(full_s, features = "pct_mito")
FeaturePlot(full_s, features = "pct_ribo")
dev.off()

cardiac_s <- readRDS(file.path(getwd(), "final_cardiac_seurat.rds"))
```

``` r
cardiac_s <- integrate.and.umap(seurat = cardiac_s,
                               nfeat = 2000,
                               cc_regress = FALSE,
                               pca_npcs = 50,
                               cumm_var = 80,
                               verb = FALSE)

pdf(file.path(plots_dir, "16_Cardiac_no_cc_reg.pdf"), height = 7, width = 12)
UMAPPlot(cardiac_s, group.by = "Broad_Annotation")
UMAPPlot(cardiac_s, group.by = "Broad_Annotation", label=TRUE)
UMAPPlot(cardiac_s, group.by = "Specific_Annotation")
UMAPPlot(cardiac_s, group.by = "Specific_Annotation", label=TRUE)
UMAPPlot(cardiac_s, group.by = "orig.ident")
UMAPPlot(cardiac_s, group.by = "Phase") + ggtitle("Phase - No CC Regression")
FeaturePlot(cardiac_s, features = "nCount_RNA")
FeaturePlot(cardiac_s, features = "nFeature_RNA")
FeaturePlot(cardiac_s, features = "pct_mito")
FeaturePlot(cardiac_s, features = "pct_ribo")
dev.off()
```

``` r
myo_s <- readRDS(file.path(getwd(), "final_myo_seurat.rds"))

myo_s <- integrate.and.umap(seurat = myo_s,
                               nfeat = 2000,
                               cc_regress = FALSE,
                               pca_npcs = 50,
                               cumm_var = 80,
                               verb = FALSE)

pdf(file.path(plots_dir, "17_Myocardial_no_cc_reg.pdf"), height = 7, width = 12)
UMAPPlot(myo_s, group.by = "myo_annot")
UMAPPlot(myo_s, group.by = "myo_annot", label=TRUE)
UMAPPlot(myo_s, group.by = "orig.ident")
UMAPPlot(myo_s, group.by = "Phase") + ggtitle("Phase - No CC Regression")
FeaturePlot(myo_s, features = "nCount_RNA")
FeaturePlot(myo_s, features = "nFeature_RNA")
FeaturePlot(myo_s, features = "pct_mito")
FeaturePlot(myo_s, features = "pct_ribo")
dev.off()
```

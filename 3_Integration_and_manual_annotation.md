Integration & Manual Annotation
================
Aaron Mohammed

``` r
library(Seurat)
library(ggplot2)
library(qpdf)
library(openxlsx)
library(plotly)

seurat_input_dir <- file.path(dirname(getwd()),
                              "2_singleR_annotation",
                              "singleR_annotated_seurats")

plots_dir <- file.path(getwd(), "Plots")
dir.create(plots_dir)

markers_dir <- file.path(getwd(), "Markers")
dir.create(markers_dir)

RDS_dir <- file.path(getwd(), "RDS_files")
dir.create(RDS_dir)
```

``` r
# Read in WT and KO seurats
WT_s <- readRDS(file.path(seurat_input_dir, "WT_singleR_annotated.rds"))
KO_s <- readRDS(file.path(seurat_input_dir, "KO_singleR_annotated.rds"))

WT_s <- AddMetaData(WT_s, rep(NA, nrow(WT_s@meta.data)), col.name = "DF.classifications")
WT_s$DF.classifications_0.25_0.17_392 -> WT_s$DF.classifications

WT_s <- AddMetaData(WT_s, rep(NA, nrow(WT_s@meta.data)), col.name = "DF.pANN")
WT_s$pANN_0.25_0.17_392 -> WT_s$DF.pANN

KO_s <- AddMetaData(KO_s, rep(NA, nrow(KO_s@meta.data)), col.name = "DF.classifications")
KO_s$DF.classifications_0.25_0.18_315 -> KO_s$DF.classifications

KO_s <- AddMetaData(KO_s, rep(NA, nrow(KO_s@meta.data)), col.name = "DF.pANN")
KO_s$pANN_0.25_0.18_315 -> KO_s$DF.pANN

Idents(WT_s) <- "orig.ident"
Idents(KO_s) <- "orig.ident"
```

``` r
# Merge WT and KO seurats
merged_s <- merge(x = WT_s, 
                  y = KO_s,
                  add.cell.ids = c("WT",
                                   "KO"),
                  project = "scRNAseq_Ren")

# Re-organize the merged seurat
merged_s$RNA_snn_res.1 <- NULL
merged_s$seurat_clusters <- NULL

# Re-order columns in meta.data so that the doubletFinder columns are grouped 
# together
merged_s@meta.data <- merged_s@meta.data[,c("orig.ident",
                                            "nCount_RNA", 
                                            "nFeature_RNA", 
                                            "pct_mito", 
                                            "pct_ribo", 
                                            "S.Score", 
                                            "G2M.Score", 
                                            "Phase", 
                                            "pANN_0.25_0.17_392", 
                                            "DF.classifications_0.25_0.17_392",
                                            "pANN_0.25_0.18_315", 
                                            "DF.classifications_0.25_0.18_315",
                                            "DF.pANN",
                                            "DF.classifications",
                                            "scDblFinder.class",
                                            "scDblFinder.score",
                                            "scDblFinder.weighted",
                                            "scDblFinder.difficulty",
                                            "scDblFinder.cxds_score",
                                            "scDblFinder.mostLikelyOrigin",
                                            "scDblFinder.originAmbiguous",
                                            "dbl_both_methods",
                                            "Cao.Spielmann_celltypes", 
                                            "Feng_celltypes",
                                            "deSoysa_celltypes", 
                                            "All_3_refs",
                                            "All_3_refs_renamed")]
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
  
  for (r in res) {
    
    seurat <- FindClusters(seurat, 
                           resolution = r,
                           verbose = FALSE)
    
    title <- paste0("Seurat Clusters", sep= " - ",
                    PC_80_pct, sep= " ",
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
}
```

``` r
get.all.markers <- function (seurat, ident, group, path, fname) {
  
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
}
```

``` r
merged_s <- integrate.and.umap(seurat = merged_s,
                               nfeat = 2000,
                               cc_regress = TRUE,
                               pca_npcs = 50,
                               cumm_var = 80,
                               verb = FALSE)
                        
# saveRDS(merged_s, file.path(RDS_dir, "0_integrated_seurat.rds"))
```

``` r
All_3_umap <- UMAPPlot(merged_s, group.by = "All_3_refs_renamed") + 
                       ggtitle("All 3 Refs Renamed")

CS_umap <- UMAPPlot(merged_s, group.by = "Cao.Spielmann_celltypes") + 
                    ggtitle("Cao - Spielmann") +
                    guides(color = guide_legend(override.aes = list(size=4), ncol=1))

Feng_umap <- UMAPPlot(merged_s, group.by = "Feng_celltypes") + 
                      ggtitle("Feng")

deSoysa_umap <- UMAPPlot(merged_s, group.by = "deSoysa_celltypes") + 
                         ggtitle("de Soysa")

pdf(file.path(plots_dir, "1_UMAPs_singleR_annotated_CTs.pdf"), height = 7, width = 12)
All_3_umap
CS_umap
Feng_umap
deSoysa_umap
dev.off()

All_3_umap
CS_umap
Feng_umap
deSoysa_umap
```
<p align="center">
  <img width="800" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/e152044f-f8c5-492c-bd28-de643d554ee8">
  <img width="800" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/a3c4a614-f17f-46bc-806d-dcb8e739634e">
  <img width="800" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/76e2ac8b-659d-4e91-9461-319f04160965">
  <img width="800" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/5edd749d-dfc8-49a8-934c-00fb1fa7f175">
</p>

``` r
DefaultAssay(merged_s) <- "RNA"
merged_s <- ScaleData(merged_s, features = row.names(merged_s))

# Blood markers
pdf(file.path(plots_dir, "2_blood_markers_fp.pdf"), height = 7, width = 12)
FeaturePlot(merged_s, features = "Hba-a1") # Erythroid
FeaturePlot(merged_s, features = "Hba-a2") # Erythroid
FeaturePlot(merged_s, features = "Hbb-bt") # Erythroid
FeaturePlot(merged_s, features = "Alas2") # Erythroid
FeaturePlot(merged_s, features = "Fcgr1") # Immune
FeaturePlot(merged_s, features = "C1qa") # Immune
FeaturePlot(merged_s, features = "Cd68") # Immune
FeaturePlot(merged_s, features = "Rgs18") # Platelet
FeaturePlot(merged_s, features = "Pf4") # Platelet
FeaturePlot(merged_s, features = "Tubb1") # Platelet
dev.off()
```
<p align="center">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/b860e269-ee55-4c56-99cc-2977af5a85a2">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/a0dbaf63-73d0-4804-beb8-e005eb9820bb">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/a4d8e7d6-b0d5-4702-95ac-c82382de092b">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/c5c14bc4-47bb-4fd4-a866-27c19ab7a546">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/2e224309-274f-4e8b-b6dd-c4682b9e3167">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/393137c5-d58d-4d3c-b0a5-00dcbd8ae8cf">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/a0e31047-40f4-4475-875f-386453cdb300">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/fd35ca63-2d52-4cf7-a882-cf47d9ccf831">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/6a2516ff-c68c-40d4-90f5-74866a029e14">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/11cce9dc-b58f-4570-9b31-ad6ad4541ec4">
</p>

``` r
# Remove blood cells

blood_cells <- CellSelector(All_3_umap)
blood_cells <- c(blood_cells, CellSelector(All_3_umap))

# saveRDS(blood_cells, file.path(RDS_dir, "1_blood_cells.rds"))

merged_s <- AddMetaData(merged_s, rep(NA, nrow(merged_s@meta.data)), col.name = "blood")

TRUE -> merged_s@meta.data[blood_cells, "blood"]
FALSE -> merged_s$blood[which(is.na(merged_s$blood))]

Idents(merged_s) <- "blood"
merged_s <- subset(merged_s, idents = FALSE)
```

``` r
# DoubletFinder classifications
df_umap <- UMAPPlot(merged_s, group.by = "DF.classifications") + 
                    ggtitle("Doublets according to DoubletFinder")
df_umap

# scDblFinder classifications
scdbl_umap_1 <- UMAPPlot(merged_s, group.by = "scDblFinder.class") + 
                ggtitle("Doublets according to scDblFinder")
scdbl_umap_1

# The doublet classifications of DoubletFinder appear mostly in 1 cluster,
# while the classifications from scDblFinder are more spread out among
# the clusters. Removing doublets according to scDblFinder.
Idents(merged_s) <- "scDblFinder.class"
merged_s <- subset(merged_s, idents = "singlet")

g <- ggplot_build(scdbl_umap_1)
cols <- unique(g$data[[1]]["colour"])$colour

scdbl_umap_2 <- UMAPPlot(merged_s, group.by = "scDblFinder.class") +
                         scale_color_manual(values = cols[1]) + 
                         ggtitle("Doublets from scDblFinder removed")

scdbl_umap_1
scdbl_umap_2 

# Removing some outlier cells that were clustered with the doublets that were
# removed
clustered_w_dbl_cells <- CellSelector(scdbl_umap_2)
clustered_w_dbl_cells <- unique(c(clustered_w_dbl_cells, CellSelector(scdbl_umap_2)))

merged_s <- AddMetaData(merged_s, rep(NA, nrow(merged_s@meta.data)), col.name = "outliers_clustered_w_scdbl")
1 -> merged_s@meta.data[clustered_w_dbl_cells, "outliers_clustered_w_scdbl"]
0 -> merged_s$outliers_clustered_w_scdbl[which(is.na(merged_s$outliers_clustered_w_scdbl))]

# saveRDS(clustered_w_dbl_cells, file.path(RDS_dir, "1_cells_clustered_w_scDblFinder_doublets.rds"))

scdbl_umap_1
scdbl_umap_2 
scdbl_umap_3 <- UMAPPlot(merged_s, group.by = "outliers_clustered_w_scdbl") +
                         scale_color_manual(values = cols) + 
                         ggtitle("Outliers were clustered with scDblFinder doublets")

Idents(merged_s) <- "outliers_clustered_w_scdbl"
merged_s <- subset(merged_s, idents = 0)

scdbl_umap_4 <- UMAPPlot(merged_s, group.by = "scDblFinder.class") +
                         scale_color_manual(values = cols[1]) + 
                         ggtitle("Final")

pdf(file.path(plots_dir, "3_doublets_DoubletFinder.pdf"), height = 7, width = 12)
df_umap
dev.off()

pdf(file.path(plots_dir, "4_doublets_scDblFinder.pdf"), height = 7, width = 12)
scdbl_umap_1
scdbl_umap_2
scdbl_umap_3 
scdbl_umap_4 
dev.off()
```
<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/56010bd0-c2c4-46eb-a2fe-e4fc8344e365">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/6ee2f22f-153a-4c73-8730-c50115c7a65c">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/741cb697-d1e2-4dd0-8e70-bb5b25db079f">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/d9d6d12f-04f5-48b8-bbbd-c7c9b1d1c6a9">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/1937acb1-0699-4f44-acea-1ca001fac322">
</p>

``` r
merged_s <- integrate.and.umap(seurat = merged_s,
                               nfeat = 2000,
                               cc_regress = TRUE,
                               pca_npcs = 50,
                               cumm_var = 80,
                               verb = FALSE)
```

``` r
All_3_umap_2 <- UMAPPlot(merged_s, group.by = "All_3_refs_renamed") + 
                       ggtitle("All 3 Refs Renamed")

CS_umap_2 <- UMAPPlot(merged_s, group.by = "Cao.Spielmann_celltypes") + 
                    ggtitle("Cao - Spielmann") +
                    guides(color = guide_legend(override.aes = list(size=4), ncol=1))

Feng_umap_2 <- UMAPPlot(merged_s, group.by = "Feng_celltypes") + 
                      ggtitle("Feng")

deSoysa_umap_2 <- UMAPPlot(merged_s, group.by = "deSoysa_celltypes") + 
                         ggtitle("de Soysa")

pdf(file.path(plots_dir, "5_singleR_CTs_after_removing_scDblfinder_dbls_and_blood.pdf"), height = 7, width = 12)
All_3_umap_2
CS_umap_2
Feng_umap_2
deSoysa_umap_2
dev.off()
```

<p align="center">
  <img width="800" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/d98f3cb3-6540-424b-b49d-2a2e3b17febd">
  <img width="800" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/d169da99-b51e-42ca-8fd7-664fbc937554">
  <img width="800" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/a8ad1ae0-5a2e-4050-b546-c35442ccb782">
  <img width="800" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/0a1c90f7-e5dc-49fa-9ad8-ddf505b18fcc">
</p>

``` r
res <- c(seq(0.1,0.9,by=0.1),seq(1,3,by=0.2))

merged_s <- FindClusters(merged_s, 
                         resolution = res,
                         verbose = FALSE)

sweep.cluster.res(seurat = merged_s,
                  res = res,
                  path = plots_dir,
                  fname = "6_cluster_res_sweep.pdf")

merged_s <- FindClusters(merged_s, 
                         resolution = 2.2,
                         verbose = FALSE)

clust_25_cells_2.2_res <- names(merged_s$seurat_clusters[which(merged_s$seurat_clusters == 25)])

merged_s <- FindClusters(merged_s, 
                         resolution = 0.4,
                         verbose = FALSE)

merged_s <- AddMetaData(merged_s, 
                        merged_s$seurat_clusters, 
                        col.name = "clusters")

levels(merged_s$clusters) <- c(levels(merged_s$clusters), 18)

18 -> merged_s$clusters[clust_25_cells_2.2_res]

cluster_umap <- UMAPPlot(merged_s, group.by = "clusters", label = TRUE)

pdf(file.path(plots_dir, "7_cluster_umap.pdf"), height = 7, width = 12)
cluster_umap
dev.off()

# saveRDS(merged_s, file.path(RDS_dir, "1_seurat_after_removing_scDblfinder_doublets_and_blood.rds"))
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/34a28464-dd2b-4e8d-a68e-683db3eaef0e">
</p>

``` r
get.cons.markers(seurat = merged_s,
                 ident = "clusters",
                 group = "orig.ident",
                 path = markers_dir,
                 fname = "1_conserved_markers_of_clusters.xlsx")

get.all.markers(seurat = merged_s,
                ident = "clusters",
                group = "orig.ident",
                path = markers_dir,
                fname = "2_all_markers_of_clusters.xlsx")
```

``` r
# Violin plots for determining if there are clusters made up of damaged cells
vln_1 <- VlnPlot(merged_s,
                 features = "nCount_RNA",
                 group.by = "clusters")

vln_2 <- VlnPlot(merged_s,
                 features = "nFeature_RNA",
                 group.by = "clusters")

vln_3 <- VlnPlot(merged_s,
                 features = "pct_mito",
                 group.by = "clusters")

pdf(file.path(plots_dir, "8_violin_plots_to_check_cell_quality.pdf"), height = 7, width = 12)
vln_1
vln_2
vln_3
dev.off()
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/789ada07-0d92-4814-86a6-78b64e9f36c0">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/c9dcddb0-08a3-4c67-ad88-2aa6fdb099b0">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/b8351f18-4b3d-4380-ab17-fd80d8146ec9">
</p>

``` r
# Removal of cells clustered in 6, 7, 10, and 18: Those cells appear to be made 
# up of damaged cells based on low UMI counts, low number of expressed genes, and 
# high expression of mitochondrial genes

DefaultAssay(merged_s) <- "integrated"
Idents(merged_s) <- "clusters"

nrow(merged_s@meta.data) # 9766
length(merged_s$orig.ident[which(merged_s$orig.ident == "WT")]) # 5261
length(merged_s$orig.ident[which(merged_s$orig.ident == "KO")]) # 4505

damaged_cells <- CellSelector(cluster_umap)

length(damaged_cells) # 1341

# saveRDS(damaged_cells, file.path(RDS_dir, "2_damaged_cells.rds"))

merged_s <- AddMetaData(merged_s, rep(NA, nrow(merged_s@meta.data)), col.name = "damaged_cells")

TRUE -> merged_s@meta.data[damaged_cells, "damaged_cells"]
FALSE -> merged_s$damaged_cells[which(is.na(merged_s$damaged_cells))]

Idents(merged_s) <- "damaged_cells"
merged_s <- subset(merged_s, idents = FALSE)

nrow(merged_s@meta.data) # 8425
length(merged_s$orig.ident[which(merged_s$orig.ident == "WT")]) # 4578
length(merged_s$orig.ident[which(merged_s$orig.ident == "KO")]) # 3847

removed_dmg_cells_umap <- UMAPPlot(merged_s, group.by = "All_3_refs_renamed")

pdf(file.path(plots_dir, "9_umap_after_removing_damaged_cells.pdf"), height = 7, width = 12)
removed_dmg_cells_umap
dev.off()
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/9b2d25d3-8d2d-402c-8783-0df881e03bcd">
</p>


``` r
# Removal of cluster 15: This cluster appears to be hepatocytes based on
# DGE markers results
DefaultAssay(merged_s) <- "RNA"

pdf(file.path(plots_dir, "10_hepatocyte_markers_fp.pdf"), height = 7, width = 12)
FeaturePlot(merged_s, features = "Alb")
FeaturePlot(merged_s, features = "Afp")
FeaturePlot(merged_s, features = "Apoa1")
FeaturePlot(merged_s, features = "Apom")
FeaturePlot(merged_s, features = "Hnf4a")
dev.off()

hepatocyte_cells <- CellSelector(cluster_umap)

length(hepatocyte_cells) # 54

# saveRDS(hepatocyte_cells, file.path(RDS_dir, "2_hepatocyte_cells.rds"))

merged_s <- AddMetaData(merged_s, rep(NA, nrow(merged_s@meta.data)), col.name = "hepatocyte_cells")

TRUE -> merged_s@meta.data[hepatocyte_cells, "hepatocyte_cells"]
FALSE -> merged_s$hepatocyte_cells[which(is.na(merged_s$hepatocyte_cells))]

Idents(merged_s) <- "hepatocyte_cells"
merged_s <- subset(merged_s, idents = FALSE)

nrow(merged_s@meta.data) # 8371
length(merged_s$orig.ident[which(merged_s$orig.ident == "WT")]) # 4578
length(merged_s$orig.ident[which(merged_s$orig.ident == "KO")]) # 3793
```

<p align="center">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/c8c5a08e-b2ea-4583-aaa8-c680317fcbe8">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/5bf9527c-0527-4aff-8750-2cbe688ef67d">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/59187ea4-93e3-4731-a22d-09880b9a646c">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/a2a005cc-e948-4ea8-ac3c-8ef5aafb87a9">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/bc371cce-3411-492e-9b34-b08604ce14b6">
</p>

``` r
merged_s <- integrate.and.umap(seurat = merged_s,
                               nfeat = 2000,
                               cc_regress = TRUE,
                               pca_npcs = 50,
                               cumm_var = 80,
                               verb = FALSE)

umap_after_reintegration <- UMAPPlot(merged_s, group.by = "All_3_refs_renamed")

pdf(file.path(plots_dir, "11_reintegration_after_removing_hepatocytes_and_damaged_cells.pdf"), height = 7, width = 12)
umap_after_reintegration
dev.off()

# saveRDS(merged_s, file.path(RDS_dir, "2_seurat_after_removing_damaged_cells_and_hepatocytes.rds"))
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/2d12abc1-e003-4240-898f-b502b65fc069">
</p>

``` r
outliers <- CellSelector(umap_after_reintegration)

length(outliers) #14

merged_s <- AddMetaData(merged_s, rep(NA, nrow(merged_s@meta.data)), col.name = "outliers")

1 -> merged_s@meta.data[outliers, "outliers"]
0 -> merged_s$outliers[which(is.na(merged_s$outliers))]

outliers_highlighted <- UMAPPlot(merged_s, group.by = "outliers")

g <- ggplot_build(outliers_highlighted)
cols <- unique(g$data[[1]]["colour"])$colour

outliers_highlighted <- outliers_highlighted + scale_color_manual(values = cols[c(2,1)])

Idents(merged_s) <- "outliers"
outlier_markers <- FindMarkers(merged_s, 
                               ident.1 = 1,
                               ident.2 = 0,
                               logfc.threshold = 0.1)

write.xlsx(outlier_markers, 
           file.path(markers_dir, "3_small_outlier_cluster_markers.xlsx"), 
           rowNames = TRUE)

pdf(file.path(plots_dir, "12_small_outlier_cluster.pdf"), height = 7, width = 12)
outliers_highlighted
dev.off()
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/38ab621b-c68b-48ac-b109-5fe0265ccaea">
</p>

``` r
# Removal of small outlier cluster: Based on the DGE markers of that cluster,
# those cells appear to be smooth muscle cells

# saveRDS(outliers, file.path(RDS_dir, "3_SMC_cells.rds"))

# Smooth muscle cell markers
DefaultAssay(merged_s) <- "RNA"
merged_s <- ScaleData(merged_s, features = row.names(merged_s))

pdf(file.path(plots_dir, "13_SMC_markers.pdf"), height = 7, width = 12)
FeaturePlot(merged_s, features = "Hhip")
FeaturePlot(merged_s, features = "Myh11")
FeaturePlot(merged_s, features = "Hpse2")
FeaturePlot(merged_s, features = "Mylk")
dev.off()

merged_s <- subset(merged_s, idents = 0)

Idents(merged_s) <- "outliers"

nrow(merged_s@meta.data) # 8357
length(merged_s$orig.ident[which(merged_s$orig.ident == "WT")]) # 4564
length(merged_s$orig.ident[which(merged_s$orig.ident == "KO")]) # 3793
```

<p align="center">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/d1c3d38f-a772-49b3-a25f-f62ac7ced0ef">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/4b2561b3-9a1f-4825-9760-50355a593302">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/b452691c-ee9a-4133-a575-db0e6e6e8e33">
  <img width="500" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/ce98b1ec-0542-484b-85cf-000e24b30d0b">
</p>


``` r
merged_s <- integrate.and.umap(seurat = merged_s,
                               nfeat = 2000,
                               cc_regress = TRUE,
                               pca_npcs = 50,
                               cumm_var = 80,
                               verb = FALSE)

umap_after_reintegration_2 <- UMAPPlot(merged_s, group.by = "All_3_refs_renamed")

pdf(file.path(plots_dir, "14_reintegration_after_removing_SMC_cells.pdf"), height = 7, width = 12)
umap_after_reintegration_2
dev.off()

# saveRDS(merged_s, file.path(RDS_dir, "3_seurat_after_removing_SMC_cells.rds"))
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/cb5a5594-0595-4c43-857e-7cd6e599733f">
</p>


``` r
res <- c(seq(0.1,0.9,by=0.1),seq(1,3,by=0.2))

merged_s <- FindClusters(merged_s, 
                         resolution = res,
                         verbose = FALSE)

sweep.cluster.res(seurat = merged_s,
                  res = res,
                  path = plots_dir,
                  fname = "15_cluster_res_sweep.pdf")

merged_s <- FindClusters(merged_s,
                         resolution = 0.3,
                         verbose = FALSE)

# UMAPPlot(merged_s, group.by = "seurat_clusters", label = TRUE)

merged_s <- AddMetaData(merged_s, 
                        merged_s$seurat_clusters, 
                        col.name = "clusters_2")
```

``` r
get.cons.markers(seurat = merged_s,
                 ident = "clusters_2",
                 group = "orig.ident",
                 path = markers_dir,
                 fname = "4_conserved_markers.xlsx")

get.all.markers(seurat = merged_s,
                ident = "clusters_2",
                group = "orig.ident",
                path = markers_dir,
                fname = "5_all_markers_of_clusters.xlsx")
```

``` r
# Myocardium markers
pdf(file.path(plots_dir, "16_myocardial_markers.pdf"), height = 7, width = 12)
FeaturePlot(merged_s, features = "Actc1")
FeaturePlot(merged_s, features = "Myl4")
FeaturePlot(merged_s, features = "Tnnt2")
FeaturePlot(merged_s, features = "Nebl")
FeaturePlot(merged_s, features = "Ttn")
FeaturePlot(merged_s, features = "Tnni3k")
FeaturePlot(merged_s, features = "Tnnc1")
FeaturePlot(merged_s, features = "Acta2")
FeaturePlot(merged_s, features = "Cryab")
dev.off()

# Endocardium markers
pdf(file.path(plots_dir, "17_endocardial_markers.pdf"), height = 7, width = 12)
FeaturePlot(merged_s, features = "Pecam1")
FeaturePlot(merged_s, features = "Egfl7")
FeaturePlot(merged_s, features = "Emcn")
FeaturePlot(merged_s, features = "Plxnd1")
FeaturePlot(merged_s, features = "Plvap")
FeaturePlot(merged_s, features = "Cdh5")
FeaturePlot(merged_s, features = "Icam2")
FeaturePlot(merged_s, features = "Ecscr")
FeaturePlot(merged_s, features = "Npr3")
FeaturePlot(merged_s, features = "Klf2")
dev.off()

# Epicardium markers
pdf(file.path(plots_dir, "18_epicardial_markers.pdf"), height = 7, width = 12)
FeaturePlot(merged_s, features = "Wt1")
FeaturePlot(merged_s, features = "Upk3b")
FeaturePlot(merged_s, features = "Aldh1a2")
FeaturePlot(merged_s, features = "Tbx18")
FeaturePlot(merged_s, features = "Sparc")
FeaturePlot(merged_s, features = "Upk1b")
FeaturePlot(merged_s, features = "Tmem255a")
FeaturePlot(merged_s, features = "Kcne1l")
dev.off()

# Mesenchymal markers
pdf(file.path(plots_dir, "19_mesenchymal_markers.pdf"), height = 7, width = 12)
FeaturePlot(merged_s, features = "Postn")
FeaturePlot(merged_s, features = "Cthrc1")
FeaturePlot(merged_s, features = "Sox9")
FeaturePlot(merged_s, features = "Pdgfra")
FeaturePlot(merged_s, features = "Papss2")
dev.off()

# Multipotent Progenitor markers
pdf(file.path(plots_dir, "20_multipotent_progenitor_markers.pdf"), height = 7, width = 12)
FeaturePlot(merged_s, features = "Osr1")
FeaturePlot(merged_s, features = "Foxf1")
FeaturePlot(merged_s, features = "Rgs5")
FeaturePlot(merged_s, features = "Isl1")
FeaturePlot(merged_s, features = "Tbx1")
FeaturePlot(merged_s, features = "Fgf10")
dev.off()
```

``` r
# Annotation 
myo_cells <- CellSelector(umap_after_reintegration_2)
endo_cells <- CellSelector(umap_after_reintegration_2)
epi_cells <- CellSelector(umap_after_reintegration_2)
mes_cells <- CellSelector(umap_after_reintegration_2)
mp_cells <- CellSelector(umap_after_reintegration_2)

merged_s <- AddMetaData(merged_s, rep(NA, nrow(merged_s@meta.data)), col.name = "broad_celltypes")

"Myocardial" -> merged_s@meta.data[myo_cells, "broad_celltypes"]
"Endocardial" -> merged_s@meta.data[endo_cells, "broad_celltypes"]
"Epicardial" -> merged_s@meta.data[epi_cells, "broad_celltypes"]
"Mesenchymal" -> merged_s@meta.data[mes_cells, "broad_celltypes"]
"MP" -> merged_s@meta.data[mp_cells, "broad_celltypes"]

broad_annotation_umap <- UMAPPlot(merged_s, 
                                  group.by = "broad_celltypes", 
                                  label = TRUE)

broad_annotation_umap <- broad_annotation_umap + ggtitle("Broad Cell Types")

condition_umap <- UMAPPlot(merged_s, group.by = "orig.ident") + ggtitle("Conditions")

phase_umap_cc_regressed <- UMAPPlot(merged_s, group.by = "Phase") + ggtitle("Phases - CC Regressed")

pdf(file.path(plots_dir, "21_broad_celltypes.pdf"), height = 7, width = 12)
broad_annotation_umap
broad_annotation_umap <- broad_annotation_umap + ggtitle("Broad Cell Types")
dev.off()

pdf(file.path(plots_dir, "22_conditions.pdf"), height = 7, width = 12)
condition_umap
dev.off()

pdf(file.path(plots_dir, "23_phase_cc_regressed.pdf"), height = 7, width = 12)
phase_umap_cc_regressed
dev.off()
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/f4c7a2a7-f1cf-4d28-8ceb-76af9b20e6ab">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/73135d40-cf8f-4fba-9188-7b0544d2c8c0">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/a94b4ee7-4ff7-4ab1-be95-2e44602255d1">
</p>

``` r
merged_s <- integrate.and.umap(seurat = merged_s,
                               nfeat = 2000,
                               cc_regress = FALSE,
                               pca_npcs = 50,
                               cumm_var = 80,
                               verb = FALSE)

phase_umap <- UMAPPlot(merged_s, group.by = "Phase") + ggtitle("Phases")
condition_umap_2 <- UMAPPlot(merged_s, group.by = "orig.ident") + ggtitle("Conditions")
broad_annotation_umap_2 <- broad_annotation_umap <- UMAPPlot(merged_s, group.by = "broad_celltypes") + ggtitle("Broad Cell Types")

pdf(file.path(plots_dir, "24_phase_CTs_condition.pdf"), height = 7, width = 12)
phase_umap
broad_annotation_umap_2
condition_umap_2
dev.off()
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/5e69aa89-3e34-46df-acf3-1ffcf9ba77a9">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/c7006a6d-33f5-491e-ae0d-36835f1e5f2f">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/51ad9dbe-25f0-4eba-a03d-90645025975b">
</p>

``` r
merged_s <- ScaleData(merged_s, 
                      vars.to.regress = c("S.Score", "G2M.Score"),
                      verbose = FALSE)

# Principal component analysis
merged_s <- RunPCA(merged_s, 
                   npcs = 50,
                   verbose = FALSE)

# Percent variation associated with each PC
pct <- merged_s[["pca"]]@stdev / sum(merged_s[["pca"]]@stdev) * 100

# Cumulative variation with increasing number of PCs
cumu <- cumsum(pct)

PC_80_pct <- which(cumu >= 80)[1] # 33 dims 80.69% cumulative variation

merged_s <- FindNeighbors(merged_s, 
                               dims = 1:PC_80_pct,
                               verbose = FALSE)

merged_s <- RunUMAP(merged_s,
                    dims = 1:PC_80_pct,
                    verbose = FALSE)

UMAPPlot(merged_s, group.by = "broad_celltypes")

# saveRDS(merged_s, file.path(RDS_dir, "4_seurat_broad_celltypes_annotated.rds"))
```

``` r
Idents(merged_s) <- "broad_celltypes"
myo_s <- subset(merged_s, idents = "Myocardial")
```

``` r
myo_s <- integrate.and.umap(seurat = myo_s,
                            nfeat = 2000,
                            cc_regress = TRUE,
                            pca_npcs = 50,
                            cumm_var = 80,
                            verb = FALSE)

# saveRDS(myo_s, file.path(RDS_dir, "5_myo_sub_clusters_seurat.rds"))
```

``` r
res <- c(seq(0.1,0.9,by=0.1),seq(1,3,by=0.2))

myo_s <- FindClusters(myo_s, 
                         resolution = res,
                         verbose = FALSE)

sweep.cluster.res(seurat = myo_s,
                  res = res,
                  path = plots_dir,
                  fname = "25_myo_cluster_res_sweep.pdf")

myo_s <- FindClusters(myo_s,
                         resolution = 0.9,
                         verbose = FALSE)

myo_s <- AddMetaData(myo_s, 
                        myo_s$seurat_clusters, 
                        col.name = "myo_clusters")


myo_clusters_umap <- UMAPPlot(myo_s, group.by = "myo_clusters", label = TRUE)

pdf(file.path(plots_dir, "26_myo_clusters.pdf"), height = 7, width = 12)
myo_clusters_umap
dev.off()
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/2b22a9fb-1971-4d57-9b1a-473882e978fd">
</p>

``` r
get.cons.markers(seurat = myo_s,
                 ident = "myo_clusters",
                 group = "orig.ident",
                 path = markers_dir,
                 fname = "6_myo_conserved_markers.xlsx")

get.all.markers(seurat = myo_s,
                ident = "myo_clusters",
                group = "orig.ident",
                path = markers_dir,
                fname = "7_all_markers_of_myo_clusters.xlsx")
```

``` r
# Load the seurat before removing scDblFinder doublets, subset myocardial 
# cluster, integrate, and then compare to myo_s

merged__s_0 <- readRDS(file.path(RDS_dir, "0_integrated_seurat.rds"))

merged__s_0 <- AddMetaData(merged__s_0, rep(NA, nrow(merged__s_0@meta.data)), col.name = "Myocardial")

TRUE -> merged__s_0@meta.data[myo_cells, "Myocardial"]
FALSE -> merged__s_0$Myocardial[is.na(which(merged__s_0$Myocardial))]

myocells_0 <- CellSelector(original_umap_myo_highlighted)

# Only keep the cells that were not annotated as MP cells
myocells_0 <- myocells_0[! myocells_0 %in% mp_cells]

merged__s_0 <- AddMetaData(merged__s_0, rep(NA, nrow(merged__s_0@meta.data)), col.name = "Myocardial")

TRUE -> merged__s_0@meta.data[myocells_0, "Myocardial"]
FALSE -> merged__s_0$Myocardial[is.na(which(merged__s_0$Myocardial))]

original_umap_myo_highlighted <- UMAPPlot(merged__s_0, group.by = "Myocardial") + 
                                 ggtitle("Original dataset with myocardial cells highlighted")
original_umap_scdblfinder_dbls <- UMAPPlot(merged__s_0, group.by = "scDblFinder.class") + 
                                  ggtitle("Original dataset with scDblFinder doublets")

Idents(merged__s_0) <- "Myocardial"
myo__s_0 <- subset(merged__s_0, idents=TRUE)
```

``` r
myo__s_0 <- integrate.and.umap(seurat = myo__s_0,
                               nfeat = 2000,
                               cc_regress = TRUE,
                               pca_npcs = 50,
                               cumm_var = 80,
                               verb = FALSE)

# UMAPPlot(myo__s_0, group.by = "Myocardial")
# UMAPPlot(myo__s_0, group.by = "scDblFinder.class")
```

``` r
myo__s_0 <- AddMetaData(myo__s_0, rep(NA, nrow(myo__s_0@meta.data)), col.name = "myo_clusters")

myo_s$myo_clusters -> myo__s_0$myo_clusters

orig_myo_clusters_w_scDblFinder_dbls <- UMAPPlot(myo__s_0, 
                                                 group.by = "myo_clusters", 
                                                 label = TRUE) + 
                                        ggtitle("Myo clusters w/ doublets that were removed (NA cells are dbls)")

orig_myo_scDblFinder_dbls <- UMAPPlot(myo__s_0, 
                                      group.by = "scDblFinder.class") + 
                             ggtitle("Myo scDblFinder doublets")

orig_myo_w_scDblFinder_dbls_scores <- FeaturePlot(myo__s_0, 
                                                  features = "scDblFinder.score") + 
                                      ggtitle("Myo scDblFinder scores")

orig_myo_conditions_w_scDblFinder_dbls <- UMAPPlot(myo__s_0, 
                                                   group.by = "orig.ident") + 
                                          ggtitle("Myo clusters w/ doublets that were removed - Conditions")

myo_cells_clust_w_scDblF_dbls <- CellSelector(orig_myo_clusters_w_scDblFinder_dbls)

myo_scdblfinder_dbls <- names(myo__s_0$myo_clusters[which(is.na(myo__s_0$myo_clusters))])

myo_cells_clust_w_scDblF_dbls_2 <- myo_cells_clust_w_scDblF_dbls[! myo_cells_clust_w_scDblF_dbls %in% myo_scdblfinder_dbls]

# saveRDS(myo_cells_clust_w_scDblF_dbls_2, file.path(RDS_dir, "7_myo_cells_that_clustered_w_scDblF_dbls.rds"))

myo__s_0 <- AddMetaData(myo__s_0, rep(NA, nrow(myo__s_0@meta.data)), col.name = "Cells_to_be_removed")

TRUE -> myo__s_0@meta.data[myo_cells_clust_w_scDblF_dbls_2, "Cells_to_be_removed"]
FALSE -> myo__s_0$Cells_to_be_removed[is.na(which(myo__s_0$Cells_to_be_removed))]

# saveRDS(myo__s_0, file.path(RDS_dir, "6_myo_clusters_w_scDblFinder_doublets_that_were_removed.rds"))

myo_cells_that_clustered_w_dbls_umap <- UMAPPlot(myo__s_0, group.by="Cells_to_be_removed")  + 
                                        ggtitle("Cells that clustered with scDblFinder doublets")
```

``` r
# Remove myo cells that clustered with scDblFinder doublets

myo_s <- AddMetaData(myo_s, rep(NA, nrow(myo_s@meta.data)), col.name = "myo_cells_clust_w_scDblF_dbls")

TRUE -> myo_s@meta.data[myo_cells_clust_w_scDblF_dbls_2, "myo_cells_clust_w_scDblF_dbls"]
FALSE -> myo_s$myo_cells_clust_w_scDblF_dbls[which(is.na(myo_s$myo_cells_clust_w_scDblF_dbls))]

myo_cells_that_will_be_removed_umap <- UMAPPlot(myo_s, group.by="myo_cells_clust_w_scDblF_dbls") + 
                                       ggtitle("Cells to be removed from myo seurat")

g <- ggplot_build(myo_cells_that_will_be_removed_umap)
cols <- unique(g$data[[1]]["colour"])$colour

myo_cells_that_will_be_removed_umap <- myo_cells_that_will_be_removed_umap + scale_color_manual(values = cols[c(2,1)])

Idents(myo_s) <- "myo_cells_clust_w_scDblF_dbls"
myo_s <- subset(myo_s, idents = FALSE)

myo_clusters_after_removal <- UMAPPlot(myo_s, group.by="myo_clusters") + 
                              ggtitle("Myo clusters after removing cells that clustered w/ scDblFinder dbls")

# saveRDS(myo_s, file.path(RDS_dir, "7_myo_seurat_after_removing_cells_that_clustered_w_dbls.rds"))

pdf(file.path(plots_dir, "27_removal_of_myo_cells_that_clustered_w_scDblF_dbls.pdf"), height = 7, width = 12)
original_umap_myo_highlighted
original_umap_scdblfinder_dbls
orig_myo_clusters_w_scDblFinder_dbls
orig_myo_scDblFinder_dbls
orig_myo_w_scDblFinder_dbls_scores
orig_myo_conditions_w_scDblFinder_dbls
myo_cells_that_clustered_w_dbls_umap
myo_cells_that_will_be_removed_umap
myo_clusters_after_removal
dev.off()
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/b4883a02-81ad-45aa-b968-1b766b9da32f">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/49ec154a-93fb-4dcb-9d1b-9173e1a5414f">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/6e2f2c4f-1d34-46a5-b330-b9b7f45e0c8c">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/4b221241-76a3-4690-860a-c8adc99daa43">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/29858c0d-a773-445c-ab0f-f95fbab6ab28">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/64662d37-2df2-435c-986d-cf0cf05e9977">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/a3b24385-a15f-46ed-bb34-c551b9b0cf5c">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/94bc58da-1a1c-4421-9e8f-5d0c3f668b4e">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/8cd6a873-f520-4eda-81af-0668859124cf">
</p>

``` r
myo_s <- integrate.and.umap(seurat = myo_s,
                            nfeat = 2000,
                            cc_regress = TRUE,
                            pca_npcs = 50,
                            cumm_var = 80,
                            verb = FALSE)
```

``` r
res <- c(seq(0.1,0.9,by=0.1),seq(1,3,by=0.2))

myo_s <- FindClusters(myo_s, 
                      resolution = res,
                      verbose = FALSE)

sweep.cluster.res(seurat = myo_s,
                  res = res,
                  path = plots_dir,
                  fname = "28_myo_cluster_res_sweep.pdf")

myo_s <- FindClusters(myo_s,
                         resolution = 0.7,
                         verbose = FALSE)

myo_s <- AddMetaData(myo_s, 
                     myo_s$seurat_clusters, 
                     col.name = "myo_clusters_2")

cells <- CellSelector(UMAPPlot(myo_s, group.by = "myo_clusters_2", label = TRUE))

7 -> myo_s@meta.data[cells, "myo_clusters_2"]

cells <- CellSelector(UMAPPlot(myo_s, group.by = "myo_clusters_2", label = TRUE))

5 -> myo_s@meta.data[cells, "myo_clusters_2"]

cells <- CellSelector(UMAPPlot(myo_s, group.by = "myo_clusters_2", label = TRUE))

3 -> myo_s@meta.data[cells, "myo_clusters_2"]



myo_clusters_umap_2 <- UMAPPlot(myo_s, group.by = "myo_clusters_2", label = TRUE) + 
                       ggtitle("Myocardial Clusters")

pdf(file.path(plots_dir, "29_myo_clusters.pdf"), height = 7, width = 12)
myo_clusters_umap_2
dev.off()

# saveRDS(myo_s, file.path(RDS_dir, "8_myo_seurat_after_reintegration.rds"))
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/6e612929-b337-48be-9f42-bcbdc2a2fb2f">
</p>

``` r
get.cons.markers(seurat = myo_s,
                 ident = "myo_clusters_2",
                 group = "orig.ident",
                 path = markers_dir,
                 fname = "8_myo_conserved_markers_after_removing_cells_that_clustered_w_dbls.xlsx")

get.all.markers(seurat = myo_s,
                ident = "myo_clusters_2",
                group = "orig.ident",
                path = markers_dir,
                fname = "9_all_markers_of_myo_after_removing_cells_that_clustered_w_dbls.xlsx")
```

``` r
# Atrial
pdf(file.path(plots_dir, "30_atrial_markers.pdf"), height = 7, width = 12)
FeaturePlot(myo_s, features = "Nppa")
FeaturePlot(myo_s, features = "Nr2f2")
FeaturePlot(myo_s, features = "Nr2f1")
FeaturePlot(myo_s, features = "Cav1")
FeaturePlot(myo_s, features = "Stard10")
FeaturePlot(myo_s, features = "Angpt1")
dev.off()

# Ventricular
pdf(file.path(plots_dir, "31_ventricular_markers.pdf"), height = 7, width = 12)
FeaturePlot(myo_s, features = "Hey2")
FeaturePlot(myo_s, features = "Myh7")
FeaturePlot(myo_s, features = "Myl2")
FeaturePlot(myo_s, features = "Mpped2")
FeaturePlot(myo_s, features = "Pln")
dev.off()

# OFT markers
pdf(file.path(plots_dir, "32_OFT_markers.pdf"), height = 7, width = 12)
FeaturePlot(myo_s, features = "Dlk1")
FeaturePlot(myo_s, features = "Meg3")
FeaturePlot(myo_s, features = "Bmp4")
FeaturePlot(myo_s, features = "Itm2a")
dev.off()

# SAN
pdf(file.path(plots_dir, "33_SAN_markers.pdf"), height = 7, width = 12)
FeaturePlot(myo_s, features = "Shox2")
FeaturePlot(myo_s, features = "Vsnl1")
FeaturePlot(myo_s, features = "Smoc2")
FeaturePlot(myo_s, features = "Hcn1")
FeaturePlot(myo_s, features = "Igfbp5")
dev.off()

# AVC
pdf(file.path(plots_dir, "34_AVC_markers.pdf"), height = 7, width = 12)
FeaturePlot(myo_s, features = "Bmp2")
FeaturePlot(myo_s, features = "Tbx3")
FeaturePlot(myo_s, features = "Tbx5")
FeaturePlot(myo_s, features = "Rspo3")
FeaturePlot(myo_s, features = "Tbx2")
dev.off()

# Early CM marker
pdf(file.path(plots_dir, "35_early_CM_marker.pdf"), height = 7, width = 12)
FeaturePlot(myo_s, features = "Gata4")
dev.off()

# Mature CM marker
pdf(file.path(plots_dir, "36_mature_CM_marker.pdf"), height = 7, width = 12)
FeaturePlot(myo_s, features = "Tnni3")
dev.off()
```

``` r
# Annotate myocardial subclusters
myo_s <- AddMetaData(myo_s, rep(NA, nrow(myo_s@meta.data)), col.name = "myocardial_celltypes")

"Ventricular" -> myo_s@meta.data[which(myo_s$myo_clusters_2 == 0), "myocardial_celltypes"]
"Ventricular" -> myo_s@meta.data[which(myo_s$myo_clusters_2 == 1), "myocardial_celltypes"]
"Ventricular" -> myo_s@meta.data[which(myo_s$myo_clusters_2 == 2), "myocardial_celltypes"]
"Ventricular" -> myo_s@meta.data[which(myo_s$myo_clusters_2 == 7), "myocardial_celltypes"]
"OFT" -> myo_s@meta.data[which(myo_s$myo_clusters_2 == 3), "myocardial_celltypes"]
"Atrial" -> myo_s@meta.data[which(myo_s$myo_clusters_2 == 4), "myocardial_celltypes"]
"SAN" -> myo_s@meta.data[which(myo_s$myo_clusters_2 == 6), "myocardial_celltypes"]
"AVC" -> myo_s@meta.data[which(myo_s$myo_clusters_2 == 5), "myocardial_celltypes"]

cells <- CellSelector(FeaturePlot(myo_s, features = "Shox2"))

"SAN" -> myo_s@meta.data[cells, "myocardial_celltypes"]

# saveRDS(myo_s, file.path(RDS_dir, "9_final_myocardial_seurat.rds"))
```

``` r
myo_CT_umap <- UMAPPlot(myo_s, group.by = "myocardial_celltypes", label = TRUE) +
               ggtitle("Myocardial Cell Types")

myo_condition_umap <- UMAPPlot(myo_s, group.by = "orig.ident") +
               ggtitle("Myocardial Conditions")

myo_phase_umap <- UMAPPlot(myo_s, group.by = "Phase", label = TRUE) +
               ggtitle("Myocardial Phases - CC regressed")

pdf(file.path(plots_dir, "37_myo_umaps.pdf"), height = 7, width = 12)
myo_CT_umap
myo_condition_umap
myo_phase_umap
dev.off()
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/4177a1d8-529e-4707-a26e-bf62f41e3e07">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/2e9d3754-15ad-41f0-8145-b9e41d6a53d6">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/f547a6fe-be45-4f04-af4f-3e143b6e99de">
</p>

``` r
get.cons.markers(seurat = myo_s,
                 ident = "myocardial_celltypes",
                 group = "orig.ident",
                 path = markers_dir,
                 fname = "10_Myocardial_conserved_markers.xlsx")

get.all.markers(seurat = myo_s,
                ident = "myocardial_celltypes",
                group = "orig.ident",
                path = markers_dir,
                fname = "11_Myocardial_all_markers.xlsx")
```

``` r
merged_s <- readRDS(file.path(RDS_dir, "4_seurat_broad_celltypes_annotated.rds"))

# Remove the myo cells that were clustered with scDblFinder doublets from the 
# main seurat
merged_s <- AddMetaData(merged_s, rep(NA, nrow(merged_s@meta.data)), col.name = "myo_cells_clust_w_scDblF_dbls")

TRUE -> merged_s@meta.data[myo_cells_clust_w_scDblF_dbls_2, "myo_cells_clust_w_scDblF_dbls"]
FALSE -> merged_s@meta.data[which(is.na(merged_s$myo_cells_clust_w_scDblF_dbls)),
                            "myo_cells_clust_w_scDblF_dbls"]

Idents(merged_s) <- "myo_cells_clust_w_scDblF_dbls"
merged_s <- subset(merged_s, idents = FALSE)

merged_s <- AddMetaData(merged_s, rep(NA, nrow(merged_s@meta.data)), col.name = "myocardial_celltypes")

myo_s$myocardial_celltypes -> merged_s$myocardial_celltypes
```

``` r
merged_s <- integrate.and.umap(seurat = merged_s,
                               nfeat = 2000,
                               cc_regress = TRUE,
                               pca_npcs = 50,
                               cumm_var = 80,
                               verb = FALSE)

broad_cts_umap_2 <- UMAPPlot(merged_s, group.by = "broad_celltypes") + ggtitle("Broad Cell Types")

cells <- CellSelector(broad_cts_umap_2)
"Epicardial" -> merged_s@meta.data[cells, "broad_celltypes"]

cells <- CellSelector(broad_cts_umap_2)
"Mesenchymal" -> merged_s@meta.data[cells, "broad_celltypes"]

cells <- CellSelector(broad_cts_umap_2)
"MP" -> merged_s@meta.data[cells, "broad_celltypes"]

cells <- CellSelector(broad_cts_umap_2)
"Myocardial" -> merged_s@meta.data[cells, "broad_celltypes"]

cells <- CellSelector(broad_cts_umap_2)
"Endocardial" -> merged_s@meta.data[cells, "broad_celltypes"]

broad_cts_umap_2 <- UMAPPlot(merged_s, group.by = "broad_celltypes") + ggtitle("Broad Cell Types")

conditions_umap_2 <- UMAPPlot(merged_s, group.by = "orig.ident") + ggtitle("Conditions")

pdf(file.path(plots_dir, "38_broad_celltypes_after_removing_myo_cells_clust_w_scdbls.pdf"), height = 7, width = 12)
broad_cts_umap_2
dev.off()

pdf(file.path(plots_dir, "39_conditions_after_removing_myo_cells_clust_w_scdbls.pdf"), height = 7, width = 12)
conditions_umap_2
dev.off()

# saveRDS(merged_s, file.path(RDS_dir, "10_final_broad_celltypes_seurat.rds"))
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/a2931b43-2dc0-4b38-ac4d-c6916fa3ad6b">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/ddb17d7b-0701-4e7f-a2dd-76650779f120">
</p>

``` r
get.cons.markers(seurat = merged_s,
                 ident = "broad_celltypes",
                 group = "orig.ident",
                 path = markers_dir,
                 fname = "12_broad_conserved_markers.xlsx")

get.all.markers(seurat = merged_s,
                ident = "broad_celltypes",
                group = "orig.ident",
                path = markers_dir,
                fname = "13_broad_all_markers.xlsx")
```

``` r
DefaultAssay(merged_s) <- "RNA"
merged_s <- ScaleData(merged_s, row.names(merged_s))

Idents(merged_s) <- "broad_celltypes"
as.factor(merged_s$broad_celltypes) -> merged_s$broad_celltypes 
c("Endocardial", "Epicardial", "Mesenchymal", "MP", "Myocardial")-> levels(merged_s$broad_celltypes)

gene_markers <- c("Actc1", "Myl4", "Nebl", "Ttn",
                  "Osr1", "Foxf1", "Rgs5", "Isl1",
                  "Postn", "Cthrc1", "Sox9", "Papss2",
                  "Wt1", "Upk3b", "Aldh1a2", "Tbx18",
                  "Pecam1", "Egfl7", "Cdh5", "Ecscr")

broad_dot_plot <- DotPlot(merged_s, features = gene_markers, group.by = "broad_celltypes", 
                          scale.by = "radius", cols = c("lightgrey", "red")) + 
                          theme(axis.text.x = element_text(angle = 45,
                                                           margin = margin(t = 12))) +
                          xlab("Markers") +
                          ylab("Broad Cell Types") 

gene_markers <- c("Pecam1", "Egfl7", "Cdh5", "Ecscr",
                  "Wt1", "Upk3b", "Aldh1a2", "Tbx18",
                  "Postn", "Cthrc1", "Sox9", "Papss2",
                  "Osr1", "Foxf1", "Rgs5", "Isl1",
                  "Actc1", "Myl4", "Nebl", "Ttn")


broad_heatmap <- DoHeatmap(merged_s, features = gene_markers, label = FALSE) + 
                           theme(axis.text = element_text(face="bold"))

pdf(file.path(plots_dir, "40_broad_CT_markers_dotplot.pdf"), height = 5, width = 10)
broad_dot_plot
dev.off()

pdf(file.path(plots_dir, "41_broad_CT_markers_heatmap.pdf"), height = 5, width = 10)
broad_heatmap
dev.off()
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/83bea2f4-b378-407f-ad49-a5710c9ed0db">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/b0bd80e0-36e2-40a0-9d63-9dc868fd91e1">
</p>

``` r
DefaultAssay(myo_s) <- "RNA"
myo_s <- ScaleData(myo_s, row.names(myo_s))

gene_markers <- c("Hey2", "Myh7", "Myl2",
                  "Vsnl1", "Shox2", "Hcn1",
                  "Dlk1", "Meg3", "Itm2a",
                  "Bmp2", "Tbx3", "Rspo3",
                  "Nppa", "Nr2f1", "Angpt1")

myo_broad_dot_plot <- DotPlot(myo_s, features = gene_markers, group.by = "myocardial_celltypes", 
                              scale.by = "radius", cols = c("lightgrey", "red")) + 
                              theme(axis.text.x = element_text(angle = 45,
                                                           margin = margin(t = 12))) +
                              xlab("Markers") +
                              ylab("Myocardial Cell Types") 

Idents(myo_s) <- "myocardial_celltypes"
as.factor(myo_s$myocardial_celltypes) -> myo_s$myocardial_celltypes 
c("Atrial", "AVC", "OFT", "SAN", "Ventricular")-> levels(myo_s$myocardial_celltypes)

gene_markers <- c("Nppa", "Nr2f1", "Angpt1",
                  "Bmp2", "Tbx3", "Rspo3",
                  "Dlk1", "Meg3", "Itm2a",
                  "Vsnl1", "Shox2", "Hcn1",
                  "Hey2", "Myh7", "Myl2")


myo_broad_heatmap <- DoHeatmap(myo_s, features = gene_markers, label = FALSE) + 
                               theme(axis.text = element_text(face="bold"))

pdf(file.path(plots_dir, "42_myo_CT_markers_dotplot.pdf"), height = 5, width = 10)
myo_broad_dot_plot
dev.off()

pdf(file.path(plots_dir, "43_myo_CT_markers_heatmap.pdf"), height = 5, width = 10)
myo_broad_heatmap
dev.off() 
```

<p align="center">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/4f5b58d3-eb99-4c6e-a1d6-a03fe89a78ca">
  <img width="1000" alt="image" src="https://github.com/ayayron117/scRNA-seq_mouse_cardiogenesis/assets/135864654/1805b36b-d2e3-457a-9d92-b3e5a9efe6f6">
</p>

``` r
writeLines(capture.output(sessionInfo()),  
           file.path(getwd(), "3_sessionInfo.txt"))
```

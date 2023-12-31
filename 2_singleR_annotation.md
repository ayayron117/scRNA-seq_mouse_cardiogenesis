Automatic Annotation Using singleR
================
Aaron Mohammed

``` r
library(Seurat)
library(SingleCellExperiment)
library(SingleR)
library(BiocParallel)

ref_dir <- file.path(getwd(), "References")

seurat_input_dir <- file.path(dirname(getwd()), 
                              "1_QC_and_CC_scoring")

singleR_results_dir <- file.path(getwd(), "SingleR_results")
dir.create(singleR_results_dir)

seurat_output_dir <- file.path(getwd(), "singleR_annotated_seurats")
dir.create(seurat_output_dir)
```

``` r
# Load references
Cao.Spielmann_sce <- readRDS(file.path(ref_dir, "Cao.Spielmann_sce_E11.rds")) 
Feng_sce <- readRDS(file.path(ref_dir, "Feng_sce_E10_E11_E12.rds")) 
deSoysa_sce <- readRDS(file.path(ref_dir, "deSoysa_sce_E9.rds"))

################################# WT singleR ###################################

WT_counts <- readRDS(file.path(getwd(), "WT_counts.rds"))

WT_pred_data <- SingleR(test = WT_counts, 
                        ref = list(Cao.Spielmann = Cao.Spielmann_sce,
                                   Feng = Feng_sce,
                                     deSoysa = deSoysa_sce),
                          labels = list(Cao.Spielmann_sce$celltype,
                                        Feng_sce$cell_type_broad,
                                        deSoysa_sce$celltype),
                          assay.type.test = "counts",
                          assay.type.ref = "logcounts",
                        BPPARAM=MulticoreParam(8))

saveRDS(WT_pred_data, file.path(singleR_results_dir, "WT_pred_data.rds"))

gc()

################################# KO singleR ###################################

KO_counts <- readRDS(file.path(getwd(), "KO_counts.rds"))

KO_pred_data <- SingleR(test = KO_counts, 
                          ref = list(Cao.Spielmann = Cao.Spielmann_sce,
                                     Feng = Feng_sce,
                                     deSoysa = deSoysa_sce),
                          labels = list(Cao.Spielmann_sce$celltype,
                                        Feng_sce$cell_type_broad,
                                        deSoysa_sce$celltype),
                          assay.type.test = "counts",
                          assay.type.ref = "logcounts",
                        BPPARAM=MulticoreParam(8))

saveRDS(KO_pred_data, file.path(singleR_results_dir, "KO_pred_data.rds")) 

gc()
```

``` r
############################### Add results to WT ##############################

WT_s <- readRDS(file.path(seurat_input_dir, 
                          "WT_seurat.rds"))

WT_s <- AddMetaData(WT_s, 
                    WT_pred_data@listData$orig.results$Cao.Spielmann$labels, 
                    col.name = "Cao.Spielmann_celltypes")

WT_s <- AddMetaData(WT_s, 
                    WT_pred_data@listData$orig.results$Feng$labels, 
                    col.name = "Feng_celltypes")

WT_s <- AddMetaData(WT_s, 
                    WT_pred_data@listData$orig.results$deSoysa$labels, 
                    col.name = "deSoysa_celltypes")

WT_s <- AddMetaData(WT_s, 
                    WT_pred_data$labels, 
                    col.name = "All_3_refs")

WT_s <- AddMetaData(WT_s, rep(NA, nrow(WT_s@meta.data)), col.name = "All_3_refs_renamed")

# unique(WT_s$All_3_refs)

"Fibroblast" -> WT_s@meta.data[which(WT_s$All_3_refs == "fibroblast-like"),"All_3_refs_renamed"]
"Endocardium" -> WT_s@meta.data[which(WT_s$All_3_refs == "endothelial"),"All_3_refs_renamed"]
"Endocardium" -> WT_s@meta.data[which(WT_s$All_3_refs == "Endothelial cells"),"All_3_refs_renamed"]
"Myocardium" -> WT_s@meta.data[which(WT_s$All_3_refs == "cardiomyocyte"),"All_3_refs_renamed"]
"Epicardium" -> WT_s@meta.data[which(WT_s$All_3_refs == "epicardial"),"All_3_refs_renamed"]
"Mesoderm" -> WT_s@meta.data[which(WT_s$All_3_refs == "Intermediate Mesoderm"),"All_3_refs_renamed"]
"Progenitor" -> WT_s@meta.data[which(WT_s$All_3_refs == "Schwann cell precursor"),"All_3_refs_renamed"]
"Erythroid" -> WT_s@meta.data[which(WT_s$All_3_refs == "blood"),"All_3_refs_renamed"]
"Immune" -> WT_s@meta.data[which(WT_s$All_3_refs == "immune"),"All_3_refs_renamed"]
"Endocardium" -> WT_s@meta.data[which(WT_s$All_3_refs == "Endocardial"),"All_3_refs_renamed"]
"Mesoderm" -> WT_s@meta.data[which(WT_s$All_3_refs == "Paraxial_mesoderm"),"All_3_refs_renamed"]
"Progenitor" -> WT_s@meta.data[which(WT_s$All_3_refs == "Multipot_progenitors"),"All_3_refs_renamed"]
"Progenitor" -> WT_s@meta.data[which(WT_s$All_3_refs == "Jaw and tooth progenitors"),"All_3_refs_renamed"]
"Progenitor" -> WT_s@meta.data[which(WT_s$All_3_refs == "Chondroctye progenitors"),"All_3_refs_renamed"]
"Erythroid" -> WT_s@meta.data[which(WT_s$All_3_refs == "Definitive erythroid lineage"),"All_3_refs_renamed"]
"Epicardium" -> WT_s@meta.data[which(WT_s$All_3_refs == "Epicardium"),"All_3_refs_renamed"]
"Myocardium" -> WT_s@meta.data[which(WT_s$All_3_refs == "Cardiac muscle lineages"),"All_3_refs_renamed"]
"Neural_Crest" -> WT_s@meta.data[which(WT_s$All_3_refs == "Neural_Crest"),"All_3_refs_renamed"]
"Megakaryocyte" -> WT_s@meta.data[which(WT_s$All_3_refs == "Megakaryocytes"),"All_3_refs_renamed"]
"Immune" -> WT_s@meta.data[which(WT_s$All_3_refs == "White blood cells"),"All_3_refs_renamed"]
"Nerve" -> WT_s@meta.data[which(WT_s$All_3_refs == "Sensory neurons"),"All_3_refs_renamed"]
"Progenitor" -> WT_s@meta.data[which(WT_s$All_3_refs == "Connective tissue progenitors"),"All_3_refs_renamed"]
"Nerve" -> WT_s@meta.data[which(WT_s$All_3_refs == "Excitatory neurons"),"All_3_refs_renamed"]
"Epicardium" -> WT_s@meta.data[which(WT_s$All_3_refs == "Epithelial cells"),"All_3_refs_renamed"]
"Nerve" -> WT_s@meta.data[which(WT_s$All_3_refs == "Cholinergic neurons"),"All_3_refs_renamed"]

saveRDS(WT_s, file.path(seurat_output_dir, "WT_singleR_annotated.rds"))
```

``` r
############################### Add results to KO ##############################

KO_s <- readRDS(file.path(seurat_input_dir, 
                          "KO_seurat.rds"))

KO_s <- AddMetaData(KO_s, 
                    KO_pred_data@listData$orig.results$Cao.Spielmann$labels, 
                    col.name = "Cao.Spielmann_celltypes")

KO_s <- AddMetaData(KO_s, 
                    KO_pred_data@listData$orig.results$Feng$labels, 
                    col.name = "Feng_celltypes")

KO_s <- AddMetaData(KO_s, 
                    KO_pred_data@listData$orig.results$deSoysa$labels, 
                    col.name = "deSoysa_celltypes")

KO_s <- AddMetaData(KO_s, 
                    KO_pred_data$labels, 
                    col.name = "All_3_refs")

KO_s <- AddMetaData(KO_s, rep(NA, nrow(KO_s@meta.data)), col.name = "All_3_refs_renamed")

# unique(KO_s$All_3_refs)

"Epicardium" -> KO_s@meta.data[which(KO_s$All_3_refs == "epicardial"),"All_3_refs_renamed"]
"Myocardium" -> KO_s@meta.data[which(KO_s$All_3_refs == "cardiomyocyte"),"All_3_refs_renamed"]
"Immune" -> KO_s@meta.data[which(KO_s$All_3_refs == "immune"),"All_3_refs_renamed"]
"Endocardium" -> KO_s@meta.data[which(KO_s$All_3_refs == "endothelial"),"All_3_refs_renamed"]
"Fibroblast" -> KO_s@meta.data[which(KO_s$All_3_refs == "fibroblast-like"),"All_3_refs_renamed"]
"Erythroid" -> KO_s@meta.data[which(KO_s$All_3_refs == "blood"),"All_3_refs_renamed"]
"Neural_Crest" -> KO_s@meta.data[which(KO_s$All_3_refs == "Neural_Crest"),"All_3_refs_renamed"]
"Endocardium" -> KO_s@meta.data[which(KO_s$All_3_refs == "Endocardial"),"All_3_refs_renamed"]
"Mesoderm" -> KO_s@meta.data[which(KO_s$All_3_refs == "Paraxial_mesoderm"),"All_3_refs_renamed"]
"Mesoderm" -> KO_s@meta.data[which(KO_s$All_3_refs == "LPM"),"All_3_refs_renamed"]
"Endocardium" -> KO_s@meta.data[which(KO_s$All_3_refs == "Endothelial cells"),"All_3_refs_renamed"]
"Epicardium" -> KO_s@meta.data[which(KO_s$All_3_refs == "Epicardium"),"All_3_refs_renamed"]
"Myocardium" -> KO_s@meta.data[which(KO_s$All_3_refs == "Cardiac muscle lineages"),"All_3_refs_renamed"]
"Megakaryocyte" -> KO_s@meta.data[which(KO_s$All_3_refs == "Megakaryocytes"),"All_3_refs_renamed"]
"Mesoderm" -> KO_s@meta.data[which(KO_s$All_3_refs == "Intermediate Mesoderm"),"All_3_refs_renamed"]
"Skeletal" -> KO_s@meta.data[which(KO_s$All_3_refs == "Chondrocytes & osteoblasts"),"All_3_refs_renamed"]
"Progenitor" -> KO_s@meta.data[which(KO_s$All_3_refs == "Multipot_progenitors"),"All_3_refs_renamed"]
"Progenitor" -> KO_s@meta.data[which(KO_s$All_3_refs == "Chondroctye progenitors"),"All_3_refs_renamed"]
"Erythroid" -> KO_s@meta.data[which(KO_s$All_3_refs == "Primitive erythroid lineage"),"All_3_refs_renamed"]
"Myocardium" -> KO_s@meta.data[which(KO_s$All_3_refs == "Myocardium"),"All_3_refs_renamed"]
"Erythroid" -> KO_s@meta.data[which(KO_s$All_3_refs == "Definitive erythroid lineage"),"All_3_refs_renamed"]
"Immune" -> KO_s@meta.data[which(KO_s$All_3_refs == "White blood cells"),"All_3_refs_renamed"]
"Epicardium" -> KO_s@meta.data[which(KO_s$All_3_refs == "Epithelial cells"),"All_3_refs_renamed"]
"Myocardium" -> KO_s@meta.data[which(KO_s$All_3_refs == "Myocytes"),"All_3_refs_renamed"]
"Mesoderm" -> KO_s@meta.data[which(KO_s$All_3_refs == "Limb mesenchyme"),"All_3_refs_renamed"]
"Hepatocyte" -> KO_s@meta.data[which(KO_s$All_3_refs == "Hepatocytes"),"All_3_refs_renamed"]

saveRDS(KO_s, file.path(seurat_output_dir, "KO_singleR_annotated.rds"))
```

``` r
writeLines(capture.output(sessionInfo()), 
           file.path(getwd(), "2_sessionInfo.txt"))
```

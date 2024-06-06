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
Feng_sce <- readRDS(file.path(ref_dir, "Feng_sce_E10_E11_E12.rds")) 
deSoysa_sce <- readRDS(file.path(ref_dir, "deSoysa_sce_E9.rds"))

################################# WT singleR ###################################

WT_counts <- readRDS(file.path(getwd(), "WT_counts.rds"))

WT_pred_data <- SingleR(test = WT_counts, 
                        ref = list(Feng = Feng_sce,
                                   deSoysa = deSoysa_sce),
                          labels = list(Feng_sce$cell_type_broad,
                                        deSoysa_sce$celltype),
                          assay.type.test = "counts",
                          assay.type.ref = "logcounts",
                        BPPARAM=MulticoreParam(8))

saveRDS(WT_pred_data, file.path(singleR_results_dir, "WT_pred_data.rds"))

gc()

################################# KO singleR ###################################

KO_counts <- readRDS(file.path(getwd(), "KO_counts.rds"))

KO_pred_data <- SingleR(test = KO_counts, 
                          ref = list(Feng = Feng_sce,
                                     deSoysa = deSoysa_sce),
                          labels = list(Feng_sce$cell_type_broad,
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
                    WT_pred_data@listData$orig.results$Feng$labels, 
                    col.name = "Feng_celltypes")

WT_s <- AddMetaData(WT_s, 
                    WT_pred_data@listData$orig.results$deSoysa$labels, 
                    col.name = "deSoysa_celltypes")

WT_s <- AddMetaData(WT_s, 
                    WT_pred_data$labels, 
                    col.name = "Both_refs")

WT_s <- AddMetaData(WT_s, rep(NA, nrow(WT_s@meta.data)), col.name = "SingleR_annots")

unique(WT_s$Both_refs)

"Erythroid" -> WT_s@meta.data[which(WT_s$Both_refs == "blood"),"SingleR_annots"]
"Fibroblast" -> WT_s@meta.data[which(WT_s$Both_refs == "fibroblast-like"),"SingleR_annots"]
"Endocardium" -> WT_s@meta.data[which(WT_s$Both_refs == "endothelial"),"SingleR_annots"]
"Endocardium" -> WT_s@meta.data[which(WT_s$Both_refs == "Endothelial cells"),"SingleR_annots"]
"Myocardium" -> WT_s@meta.data[which(WT_s$Both_refs == "cardiomyocyte"),"SingleR_annots"]
"Epicardium" -> WT_s@meta.data[which(WT_s$Both_refs == "epicardial"),"SingleR_annots"]
"Mesoderm" -> WT_s@meta.data[which(WT_s$Both_refs == "Intermediate Mesoderm"),"SingleR_annots"]
"Progenitor" -> WT_s@meta.data[which(WT_s$Both_refs == "Schwann cell precursor"),"SingleR_annots"]
"Immune" -> WT_s@meta.data[which(WT_s$Both_refs == "immune"),"SingleR_annots"]
"Endocardium" -> WT_s@meta.data[which(WT_s$Both_refs == "Endocardial"),"SingleR_annots"]
"Erythroid" -> WT_s@meta.data[which(WT_s$Both_refs == "Primitive erythroid lineage"),"SingleR_annots"]
"Mesoderm" -> WT_s@meta.data[which(WT_s$Both_refs == "Paraxial_mesoderm"),"SingleR_annots"]
"Myocardium" -> WT_s@meta.data[which(WT_s$Both_refs == "Cardiac muscle lineages"),"SingleR_annots"]
"Progenitor" -> WT_s@meta.data[which(WT_s$Both_refs == "Jaw and tooth progenitors"),"SingleR_annots"]
"Progenitor" -> WT_s@meta.data[which(WT_s$Both_refs == "Chondroctye progenitors"),"SingleR_annots"]
"Erythroid" -> WT_s@meta.data[which(WT_s$Both_refs == "Definitive erythroid lineage"),"SingleR_annots"]
"Epicardium" -> WT_s@meta.data[which(WT_s$Both_refs == "Epicardium"),"SingleR_annots"]
"Megakaryocyte" -> WT_s@meta.data[which(WT_s$Both_refs == "Megakaryocytes"),"SingleR_annots"]
"Progenitor" -> WT_s@meta.data[which(WT_s$Both_refs == "Multipot_progenitors"),"SingleR_annots"]
"Immune" -> WT_s@meta.data[which(WT_s$Both_refs == "White blood cells"),"SingleR_annots"]
"Nerve" -> WT_s@meta.data[which(WT_s$Both_refs == "Sensory neurons"),"SingleR_annots"]
"Progenitor" -> WT_s@meta.data[which(WT_s$Both_refs == "Connective tissue progenitors"),"SingleR_annots"]
"Neural_Crest" -> WT_s@meta.data[which(WT_s$Both_refs == "Neural_Crest"),"SingleR_annots"]
"Nerve" -> WT_s@meta.data[which(WT_s$Both_refs == "Excitatory neurons"),"SingleR_annots"]
"Nerve" -> WT_s@meta.data[which(WT_s$Both_refs == "Cholinergic neurons"),"SingleR_annots"]
"Epicardium" -> WT_s@meta.data[which(WT_s$Both_refs == "Epithelial cells"),"SingleR_annots"]

saveRDS(WT_s, file.path(seurat_output_dir, "WT_singleR_annotated.rds"))
```

``` r
############################### Add results to KO ##############################

KO_s <- readRDS(file.path(seurat_input_dir, 
                          "KO_seurat.rds"))

KO_s <- AddMetaData(KO_s, 
                    KO_pred_data@listData$orig.results$Feng$labels, 
                    col.name = "Feng_celltypes")

KO_s <- AddMetaData(KO_s, 
                    KO_pred_data@listData$orig.results$deSoysa$labels, 
                    col.name = "deSoysa_celltypes")

KO_s <- AddMetaData(KO_s, 
                    KO_pred_data$labels, 
                    col.name = "Both_refs")

KO_s <- AddMetaData(KO_s, rep(NA, nrow(KO_s@meta.data)), col.name = "SingleR_annots")

unique(KO_s$Both_refs)

"Epicardium" -> KO_s@meta.data[which(KO_s$Both_refs == "epicardial"),"SingleR_annots"]
"Myocardium" -> KO_s@meta.data[which(KO_s$Both_refs == "cardiomyocyte"),"SingleR_annots"]
"Immune" -> KO_s@meta.data[which(KO_s$Both_refs == "immune"),"SingleR_annots"]
"Endocardium" -> KO_s@meta.data[which(KO_s$Both_refs == "endothelial"),"SingleR_annots"]
"Fibroblast" -> KO_s@meta.data[which(KO_s$Both_refs == "fibroblast-like"),"SingleR_annots"]
"Erythroid" -> KO_s@meta.data[which(KO_s$Both_refs == "blood"),"SingleR_annots"]
"Neural_Crest" -> KO_s@meta.data[which(KO_s$Both_refs == "Neural_Crest"),"SingleR_annots"]
"Endocardium" -> KO_s@meta.data[which(KO_s$Both_refs == "Endocardial"),"SingleR_annots"]
"Mesoderm" -> KO_s@meta.data[which(KO_s$Both_refs == "Paraxial_mesoderm"),"SingleR_annots"]
"Mesoderm" -> KO_s@meta.data[which(KO_s$Both_refs == "LPM"),"SingleR_annots"]
"Endocardium" -> KO_s@meta.data[which(KO_s$Both_refs == "Endothelial cells"),"SingleR_annots"]
"Epicardium" -> KO_s@meta.data[which(KO_s$Both_refs == "Epicardium"),"SingleR_annots"]
"Myocardium" -> KO_s@meta.data[which(KO_s$Both_refs == "Cardiac muscle lineages"),"SingleR_annots"]
"Megakaryocyte" -> KO_s@meta.data[which(KO_s$Both_refs == "Megakaryocytes"),"SingleR_annots"]
"Progenitor" -> KO_s@meta.data[which(KO_s$Both_refs == "Multipot_progenitors"),"SingleR_annots"]
"Mesoderm" -> KO_s@meta.data[which(KO_s$Both_refs == "Intermediate Mesoderm"),"SingleR_annots"]
"Skeletal" -> KO_s@meta.data[which(KO_s$Both_refs == "Chondrocytes & osteoblasts"),"SingleR_annots"]
"Progenitor" -> KO_s@meta.data[which(KO_s$Both_refs == "Chondroctye progenitors"),"SingleR_annots"]
"Erythroid" -> KO_s@meta.data[which(KO_s$Both_refs == "Primitive erythroid lineage"),"SingleR_annots"]
"Myocardium" -> KO_s@meta.data[which(KO_s$Both_refs == "Myocardium"),"SingleR_annots"]
"Erythroid" -> KO_s@meta.data[which(KO_s$Both_refs == "Definitive erythroid lineage"),"SingleR_annots"]
"Immune" -> KO_s@meta.data[which(KO_s$Both_refs == "White blood cells"),"SingleR_annots"]
"Myocardium" -> KO_s@meta.data[which(KO_s$Both_refs == "Myocytes"),"SingleR_annots"]
"Mesoderm" -> KO_s@meta.data[which(KO_s$Both_refs == "Limb mesenchyme"),"SingleR_annots"]
"Hepatocyte" -> KO_s@meta.data[which(KO_s$Both_refs == "Hepatocytes"),"SingleR_annots"]

saveRDS(KO_s, file.path(seurat_output_dir, "KO_singleR_annotated.rds"))
```

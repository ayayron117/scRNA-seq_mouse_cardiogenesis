Preparation for SingleR
================
Aaron Mohammed

``` r
library(Seurat)
library(SingleCellExperiment)
library(Matrix)
library(data.table)
library(scuttle)

ref_dir <- file.path(getwd(), "References")

seurat_dir <- file.path(dirname(getwd()), 
                        "1_QC_and_CC_scoring")
```

``` r
# The file below was downloaded from here: 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193346

Feng_seurat <- readRDS(file.path(ref_dir,"GSE193346_Feng", "GSE193346_C57BL6_seurat_object.rds"))
Feng_seurat@meta.data$cell_type_broad <- as.character(Feng_seurat@meta.data$cell_type_broad)
Feng_seurat@meta.data$cell_type_spec <- as.character(Feng_seurat@meta.data$cell_type_spec)
Feng_seurat@meta.data$stage <- as.character(Feng_seurat@meta.data$stage)
Feng_seurat@meta.data$zone <- as.character(Feng_seurat@meta.data$zone)

Idents(Feng_seurat) <- "stage"
Feng_seurat <- subset(Feng_seurat, idents = c("E10.5", "E11.5", "E12.5"))
Feng_seurat[["orig.ident"]] <- "Feng"

Feng_E10_E11_E12_sce <- as.SingleCellExperiment(Feng_seurat)
Feng_E10_E11_E12_sce <- logNormCounts(Feng_E10_E11_E12_sce)

saveRDS(Feng_E10_E11_E12_sce, file.path(ref_dir, "Feng_sce_E10_E11_E12.rds"))
```

``` r
# The links below were gathered from here:
# https://cells.ucsc.edu/?ds=mouse-cardiac+wtAll

mat <- fread("https://cells.ucsc.edu/mouse-cardiac/wtAll/exprMatrix.tsv.gz")
meta <- data.frame(fread("https://cells.ucsc.edu/mouse-cardiac/wtAll/meta.tsv"), row.names=1)
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
deSoysa_seurat <- CreateSeuratObject(counts = mat, project = "deSoysa", meta.data=meta)
deSoysa_seurat <- AddMetaData(deSoysa_seurat, meta$Cluster, col.name = "celltype")
deSoysa_seurat$Time.Point <- meta$Time.Point

Idents(deSoysa_seurat) <- "Time.Point"
deSoysa_seurat <- subset(deSoysa_seurat, idents="E9.25")

deSoysa_E9_sce <- as.SingleCellExperiment(deSoysa_seurat)
deSoysa_E9_sce <- logNormCounts(deSoysa_E9_sce)

saveRDS(deSoysa_E9_sce, file.path(ref_dir, "deSoysa_sce_E9.rds"))
```

``` r
wild_s <- readRDS(file.path(seurat_dir, 
                            "WT_seurat.rds"))

wild_counts <- GetAssayData(wild_s, slot = "counts")
saveRDS(wild_counts, file.path(getwd(), "WT_counts.rds"))

knock_s <- readRDS(file.path(seurat_dir, 
                            "KO_seurat.rds"))

knock_counts <- GetAssayData(knock_s, slot = "counts")
saveRDS(knock_counts, file.path(getwd(), "KO_counts.rds"))
```

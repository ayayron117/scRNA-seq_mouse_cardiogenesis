DGE Analysis
================
Aaron Mohammed

``` r
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(KEGG.db)
library(enrichplot)
library(EnhancedVolcano)
library(openxlsx)
library(simplifyEnrichment)
library(scales)

seurat_input_dir <- file.path(dirname(getwd()),
                              "2_Annotation_and_integration")

dge_dir <- file.path(getwd(), "DGE_results")
dir.create(dge_dir)
```

``` r
myo <- readRDS(file.path(seurat_input_dir, "final_myo_seurat.rds"))
```

``` r
DefaultAssay(myo) <- "RNA"
Idents(myo) <- "myo_annot"

atrial <- subset(myo, idents = "Atrial")
atrial <- ScaleData(atrial, row.names(atrial))

ventriclular <- subset(myo, idents = c("Ven_1", "Ven_2"))
ventriclular <- ScaleData(ventriclular, row.names(ventriclular))

SAN <- subset(myo, idents = "SAN")
SAN <- ScaleData(SAN, row.names(SAN))

OFT <- subset(myo, idents = "OFT")
OFT <- ScaleData(OFT, row.names(OFT))

AVC <- subset(myo, idents = "AVC")
AVC <- ScaleData(AVC, row.names(AVC))
```

``` r
go.gsea <- function(df, up_dir, down_dir, fileName, db, pCutoff, dge_method, color_by) {
  
  set.seed(123)

  # Get fold change values with names 
  log2FC_list <- df$avg_log2FC
  names(log2FC_list) <- rownames(df)
  log2FC_list <- na.omit(log2FC_list)
  log2FC_list = sort(log2FC_list, decreasing = TRUE)
  
  # GSEA GO
  gsea_GO_BP <- gseGO(geneList=log2FC_list,
                      ont ="BP", 
                      keyType = "SYMBOL",
                      OrgDb = db,
                      pAdjustMethod = "BH",
                      pvalueCutoff = pCutoff,
                      minGSSize = 3,
                      eps = 0)
  
  gsea_GO_CC <- gseGO(geneList=log2FC_list,
                      ont ="CC", 
                      keyType = "SYMBOL",
                      OrgDb = db,
                      pAdjustMethod = "BH",
                      pvalueCutoff = pCutoff,
                      minGSSize = 3,
                      eps = 0)
  
  gsea_GO_MF <- gseGO(geneList=log2FC_list,
                      ont ="MF", 
                      keyType = "SYMBOL",
                      OrgDb = db,
                      pAdjustMethod = "BH",
                      pvalueCutoff = pCutoff,
                      minGSSize = 3,
                      eps = 0)
  
  onts <- c("BP", "CC", "MF")
  
  ############################################################################
  
  up_gsea <- file.path(up_dir, "GSEA")
  dir.create(up_gsea, 
             showWarnings = FALSE)
  
  up_gsea_go <- file.path(up_gsea, "GO")
  dir.create(up_gsea_go)
  up_gsea_go_m <- file.path(up_gsea_go, dge_method)
  dir.create(up_gsea_go_m)
  
  up_gsea_GO_list <- list(subset_enrichResult(gsea_GO_BP, 
                                              which(gsea_GO_BP@result$NES >= 0)), 
                          subset_enrichResult(gsea_GO_CC,
                                              which(gsea_GO_CC@result$NES >= 0)), 
                          subset_enrichResult(gsea_GO_MF,
                                              which(gsea_GO_MF@result$NES >= 0)))
  names(up_gsea_GO_list) <- onts
  
  up_gsea_GO_df_l <- list(up_gsea_GO_list[[1]]@result,
                          up_gsea_GO_list[[2]]@result,
                          up_gsea_GO_list[[3]]@result)
  
  names(up_gsea_GO_df_l) <- onts
  
  write.xlsx(up_gsea_GO_df_l, 
             file.path(up_gsea_go_m, 
                       paste0("UP", sep = "_",
                              fileName, sep="_", 
                              "GSEA_GO.xlsx")))

  # # # # # # # # # # # # # # # # # # # # # #
  
  for (o in onts) {
    
    if (o == "BP") {
      title <- "Biological Processes"} 
    else if (o == "MF") {
      title <- "Molecular Functions"} 
    else if (o == "CC") {
      title <- "Cellular Components"}
    
    h <- 5
    
    if (nrow(up_gsea_GO_list[[o]]@result) >= 10) {
      if (max(nchar(up_gsea_GO_list[[o]]@result$Description[1:10])) >= 120) {
        h <- h + floor(max(nchar(up_gsea_GO_list[[o]]@result$Description[1:20]))/60)
      }
    } else if (nrow(up_gsea_GO_list[[o]]@result) < 10) { 
        if (max(nchar(up_gsea_GO_list[[o]]@result$Description)) >= 120) {
          h <- h + floor(max(nchar(up_gsea_GO_list[[o]]@result$Description[1:20]))/60)
        }
    }

    
    dp <- dotplot(up_gsea_GO_list[[o]], 
                  showCategory=10, 
                  label_format = 60,
                  color = color_by) + 
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    pdf(file.path(up_gsea_go_m, 
                  paste0("UP", sep="_",
                         o, sep="_",
                         fileName, sep="_", 
                         "GSEA_GO_dotplot.pdf")), 
        height = h, width=12)
    print(dp)
    dev.off()
    
    set.seed(123)
    
    pt <- pairwise_termsim(up_gsea_GO_list[[o]])
    emap <- emapplot(pt, 
                     color = color_by) + 
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    pdf(file.path(up_gsea_go_m, 
                  paste0("UP", sep="_",
                         o, sep="_",
                         fileName, sep="_", 
                         "GSEA_GO_emap.pdf")), 
        height = 15, 
        width=15)
    print(emap)
    dev.off()
  
  }
  
  down_gsea <- file.path(down_dir, "GSEA")
  dir.create(down_gsea, 
             showWarnings = FALSE)
  
  down_gsea_go <- file.path(down_gsea, "GO")
  dir.create(down_gsea_go)
  down_gsea_go_m <- file.path(down_gsea_go, dge_method)
  dir.create(down_gsea_go_m)
  
  down_gsea_GO_list <- list(subset_enrichResult(gsea_GO_BP, 
                                              which(gsea_GO_BP@result$NES <= 0)), 
                          subset_enrichResult(gsea_GO_CC,
                                              which(gsea_GO_CC@result$NES <= 0)), 
                          subset_enrichResult(gsea_GO_MF,
                                              which(gsea_GO_MF@result$NES <= 0)))
  names(down_gsea_GO_list) <- onts
  
  down_gsea_GO_df_l <- list(down_gsea_GO_list[[1]]@result,
                            down_gsea_GO_list[[2]]@result,
                            down_gsea_GO_list[[3]]@result)
  
  names(down_gsea_GO_df_l) <- onts
  
  write.xlsx(down_gsea_GO_df_l, 
             file.path(down_gsea_go_m, 
                       paste0("DOWN", sep = "_",
                              fileName, sep="_", 
                              "GSEA_GO.xlsx")))
  
  # # # # # # # # # # # # # # # # # # # # # #
  
  for (o in onts) {
    
    if (o == "BP") {
      title <- "Biological Processes"} 
    else if (o == "MF") {
      title <- "Molecular Functions"} 
    else if (o == "CC") {
      title <- "Cellular Components"}
    
    h <- 5
    
    if (nrow(down_gsea_GO_list[[o]]@result) >= 10) {
      if (max(nchar(down_gsea_GO_list[[o]]@result$Description[1:10])) >= 120) {
        h <- h + floor(max(nchar(down_gsea_GO_list[[o]]@result$Description[1:20]))/60)
      }
    } else if (nrow(down_gsea_GO_list[[o]]@result) < 10) { 
        if (max(nchar(down_gsea_GO_list[[o]]@result$Description)) >= 120) {
          h <- h + floor(max(nchar(down_gsea_GO_list[[o]]@result$Description[1:20]))/60)
        }
    }
    
    dp <- dotplot(down_gsea_GO_list[[o]], 
                  showCategory=10, 
                  label_format = 60,
                  color = color_by) + 
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    pdf(file.path(down_gsea_go_m, 
                  paste0("DOWN", sep="_",
                         o, sep="_",
                         fileName, sep="_", 
                         "GSEA_GO_dotplot.pdf")), 
        height = h, width=12)
    print(dp)
    dev.off()
    
    set.seed(123)
    
    pt <- pairwise_termsim(down_gsea_GO_list[[o]])
    emap <- emapplot(pt, 
                     color = color_by) + 
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    pdf(file.path(down_gsea_go_m, 
                  paste0("DOWN", sep="_",
                         o, sep="_",
                         fileName, sep="_", 
                         "GSEA_GO_emap.pdf")), 
        height = 15, 
        width=15)
    print(emap)
    dev.off()
  
  }

}
```

``` r
# GSEA KEGG

kegg.gsea <- function(df, up_dir, down_dir, fileName, pCutoff, dge_method, color_by) {
  
  set.seed(123)
  
  df <- df[which(df$p_val <= pCutoff), ]
  
  # Get names of the genes
  gene_list <- row.names(df)
  
  # Get fold change values with ENTREZID as names 
  entrz_id <- bitr(gene_list, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
  log2FC_list <- df$avg_log2FC[which(gene_list %in% entrz_id$SYMBOL)]
  names(log2FC_list) <- entrz_id$ENTREZID
  log2FC_list <- na.omit(log2FC_list)
  log2FC_list = sort(log2FC_list, decreasing = TRUE)
  
  gsea_KEGG <- gseKEGG(geneList=log2FC_list,
                       organism = 'mmu',
                       use_internal_data = T,
                       pAdjustMethod = "BH",
                       pvalueCutoff = pCutoff,
                       minGSSize = 3)
  
  ############################################################################
  
  up_gsea <- file.path(up_dir, "GSEA")
  dir.create(up_gsea, 
             showWarnings = FALSE)
  
  up_gsea_kegg <- file.path(up_gsea, "KEGG")
  dir.create(up_gsea_kegg)
  up_gsea_kegg_m <- file.path(up_gsea_kegg, dge_method)
  dir.create(up_gsea_kegg_m)
  
  # # # # # # # # # # # # # # # # # # # # # #
  
  up_gsea_KEGG <- subset_enrichResult(gsea_KEGG, 
                                      which(gsea_KEGG@result$NES >= 0))
  
  write.xlsx(up_gsea_KEGG@result, 
             file.path(up_gsea_kegg_m, 
                       paste0("UP", sep = "_",
                              fileName, sep="_", 
                              "GSEA_KEGG.xlsx")))
  
  # # # # # # # # # # # # # # # # # # # # # #
  
  h <- 5
  
  if (nrow(up_gsea_KEGG@result) >= 10) {
    if (max(nchar(up_gsea_KEGG@result$Description[1:10])) >= 120) {
      h <- h + floor(max(nchar(up_gsea_KEGG@result$Description[1:20]))/60)
    }
  } else if (nrow(up_gsea_KEGG@result) < 10) { 
      if (max(nchar(up_gsea_KEGG@result$Description)) >= 120) {
        h <- h + floor(max(nchar(up_gsea_KEGG@result$Description[1:20]))/60)
      }
  }
  
  dp <- dotplot(up_gsea_KEGG, 
                showCategory=10, 
                label_format = 60,
                color = color_by)
  
  pdf(file.path(up_gsea_kegg_m, 
                paste0("UP", sep="_",
                       fileName, sep="_", 
                       "GSEA_KEGG_dotplot.pdf")), 
      height = h, width=12)
  print(dp)
  dev.off()
  
  set.seed(123)
  
  pt <- pairwise_termsim(up_gsea_KEGG)
  emap <- emapplot(pt, 
                   color = color_by)
  
  pdf(file.path(up_gsea_kegg_m, 
                paste0("UP", sep="_",
                       fileName, sep="_", 
                       "GSEA_KEGG_emap.pdf")), 
      height = 15, 
      width=15)
  print(emap)
  dev.off()
  
  ############################################################################
  
  down_gsea <- file.path(down_dir, "GSEA")
  dir.create(down_gsea, 
             showWarnings = FALSE)
  
  down_gsea_kegg <- file.path(down_gsea, "KEGG")
  dir.create(down_gsea_kegg)
  down_gsea_kegg_m <- file.path(down_gsea_kegg, dge_method)
  dir.create(down_gsea_kegg_m)
  
  # # # # # # # # # # # # # # # # # # # # # #
  
  down_gsea_KEGG <- subset_enrichResult(gsea_KEGG, 
                                      which(gsea_KEGG@result$NES >= 0))
  
  write.xlsx(down_gsea_KEGG@result, 
             file.path(down_gsea_kegg_m, 
                       paste0("DOWN", sep = "_",
                              fileName, sep="_", 
                              "GSEA_KEGG.xlsx")))
  
  # # # # # # # # # # # # # # # # # # # # # #
  
  h <- 5
  
  if (nrow(down_gsea_KEGG@result) >= 10) {
    if (max(nchar(down_gsea_KEGG@result$Description[1:10])) >= 120) {
      h <- h + floor(max(nchar(down_gsea_KEGG@result$Description[1:20]))/60)
    }
  } else if (nrow(down_gsea_KEGG@result) < 10) { 
      if (max(nchar(down_gsea_KEGG@result$Description)) >= 120) {
        h <- h + floor(max(nchar(down_gsea_KEGG@result$Description[1:20]))/60)
      }
  }
  
  
  dp <- dotplot(down_gsea_KEGG, 
                showCategory=10, 
                label_format = 60,
                color = color_by)
  
  pdf(file.path(down_gsea_kegg_m, 
                paste0("DOWN", sep="_",
                       fileName, sep="_", 
                       "GSEA_KEGG_dotplot.pdf")), 
      height = h, width=12)
  print(dp)
  dev.off()
  
  set.seed(123)
  
  pt <- pairwise_termsim(down_gsea_KEGG)
  emap <- emapplot(pt, 
                   color = color_by)
  
  pdf(file.path(down_gsea_kegg_m, 
                paste0("DOWN", sep="_",
                       fileName, sep="_", 
                       "GSEA_KEGG_emap.pdf")), 
      height = 15, 
      width=15)
  print(emap)
  dev.off()

}
```

``` r
go.ora.plots <- function (enrichResult, ont, dir, up_down, fileName, color_by) {
  
  ont_sub <- subset_enrichResult(enrichResult, 
                                 which(enrichResult@result$ONTOLOGY==ont & 
                                         enrichResult@result$pvalue <= 0.05))
  
  if (ont == "BP") {
    title <- "Biological Processes"} 
  else if (ont == "MF") {
    title <- "Molecular Functions"} 
  else if (ont == "CC") {
    title <- "Cellular Components"}
  
  h <- 10
  
  if (max(nchar(ont_sub@result$Description[1:20])) >= 120) {
    h <- h + floor(max(nchar(ont_sub@result$Description[1:20]))/60)
  }

  bp <- barplot(ont_sub,
        showCategory = 20,
        color = color_by,
        label_format = 60) +
  ggtitle(title) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  pdf(height = h, 
      width=12,
      file.path(dir, paste0(up_down, sep="_", 
                            fileName, sep="_",
                            "ORA", sep= "_",
                            "GO", sep= "_",
                            ont, sep= "_",
                            "barplot.pdf")))
  print(bp)
  dev.off()
  
  set.seed(123)
  
  pt <- pairwise_termsim(ont_sub)
  
  emap <- emapplot(pt, color = color_by) + 
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, 
                                    face = "bold"))
  
  pdf(height = 20, 
      width=20,
      file.path(dir, paste0(up_down, sep="_", 
                            fileName, sep="_",
                            "ORA", sep= "_",
                            "GO", sep= "_",
                            ont, sep= "_",
                            "emap.pdf")))
  print(emap)
  dev.off()
   
}
```

``` r
go.ora <- function(df, up_down, dir, fileName, db, color_by, pCutoff) {
  
  set.seed(123)

  # Get names of the genes
  gene_list <- row.names(df)
  
  # ORA GO
  ora_GO <- enrichGO(gene = gene_list,
                 OrgDb= db,
                 keyType = "SYMBOL",
                 ont = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff = pCutoff) # cutoff for adjusted pvalue
  
  onts <- unique(ora_GO@result$ONTOLOGY)
  
  ora_GO_list <- vector(mode = "list", length = length(onts))
  names(ora_GO_list) <- onts
  
  for (o in onts) {
    ora_GO_list[[o]] <- ora_GO@result[which(ora_GO@result$ONTOLOGY == o), ]
  }
  
  write.xlsx(ora_GO_list, 
             file.path(dir, paste0(up_down, sep="_", 
                                   fileName, sep="_", 
                                   "ORA_GO.xlsx")))
  
  for (o in onts) {
    go.ora.plots(ora_GO, 
              o, 
              dir, 
              up_down, 
              fileName,
              color_by)
  }
  
}
```

``` r
kegg.ora <- function(df, up_down, dir, fileName, db, color_by, pCutoff) {
  
  set.seed(123)

  # Get names of the genes
  gene_list <- row.names(df)
  
  # Convert the names of the genes to ENTREZID for KEGG
  entrz_id <- bitr(gene_list, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
  
  # ORA KEGG
  ora_KEGG <- enrichKEGG(gene = entrz_id$ENTREZID,
                 organism = 'mmu',
                 use_internal_data = T,
                 pAdjustMethod = "BH",
                 pvalueCutoff = pCutoff)
  
  write.xlsx(ora_KEGG@result, file.path(dir, paste0(up_down, sep="_", 
                                                   fileName, sep="_", 
                                                   "ORA_KEGG.xlsx")))

  bp <- barplot(ora_KEGG, 
              showCategory = ora_KEGG@result$Description[1:10],
              color = color_by)
  
  pdf(file.path(dir, 
                paste0(up_down, sep="_", fileName, sep="_", "ORA_KEGG_barplot.pdf")), 
      height = 8, width=10)
  print(bp)
  dev.off()
  
  pt <- pairwise_termsim(ora_KEGG)
  
  emap <- emapplot(pt, color = color_by)
  
  pdf(file.path(dir, 
                paste0(up_down, sep="_", fileName, sep="_", "ORA_KEGG_emap.pdf")), 
      height = 15, 
      width=15)
  print(emap)
  dev.off()

}
```

``` r
ora.enrich <- function (df, dir,up_down, label, pCutoff, dge_method = NULL, uni_int = FALSE, color_by) {
  
  set.seed(123)
  
  ora <- file.path(dir, "ORA")
  dir.create(ora)
  
  ora_go <- file.path(ora, "GO")
  dir.create(ora_go)
  if (uni_int == FALSE) {
    ora_go_m <- file.path(ora_go, dge_method)
    dir.create(ora_go_m)}
  
  ora_kegg <- file.path(ora, "KEGG")
  dir.create(ora_kegg)
  if (uni_int == FALSE) {
    ora_kegg_m <- file.path(ora_kegg, dge_method)
    dir.create(ora_kegg_m)}
  
  ############################################################################
  
  # GO ORA
  if (uni_int == FALSE) {
    try(go.ora(df, up_down= up_down, dir = ora_go_m, fileName= label, 
               db= "org.Mm.eg.db", color_by = color_by, pCutoff = pCutoff))
  } else if (uni_int == TRUE) {
    try(go.ora(df, up_down= up_down, dir = ora_go, fileName= label, 
               db= "org.Mm.eg.db", color_by = color_by, pCutoff = pCutoff))
  }
  
  # KEGG ORA
  if (uni_int == FALSE) {
    try(kegg.ora(df, up_down= up_down, dir = ora_kegg_m, fileName= label, 
                 color_by = color_by, pCutoff = pCutoff))
  } else if (uni_int == TRUE) {
    try(kegg.ora(df, up_down= up_down, dir = ora_kegg, fileName= label, 
                 color_by = color_by, pCutoff = pCutoff))
  }
  
}
```

``` r
dge.analysis <- function (seurat, cluster, methods, ident.1, ident.2,
                          FCcutoff = 0, pCutoff = 0.05, path, color_by) {
  
  set.seed(123)
  
  DefaultAssay(seurat) <- "RNA"
  Idents(seurat) <- "orig.ident"
  
  parent_dir <- file.path(path, cluster)
  dir.create(parent_dir)
  
  results_list <- vector(mode = "list", length = length(methods))
  names(results_list) <- methods
  
  up_list <- vector(mode = "list", length = length(methods))
  names(up_list) <- methods
  
  down_list <- vector(mode = "list", length = length(methods))
  names(down_list) <- methods
  
  ############################################################################
  
  up <- file.path(parent_dir, "1_UP")
  dir.create(up)
  
  down <- file.path(parent_dir, "2_DOWN")
  dir.create(down)
  
  ############################################################################
  
  for (i in 1:length(methods)) {
    
    label <- paste0(methods[i], sep="_",
                    cluster)
    
    # # # # # # # # # # # # # # # # # # # # # #
  
    I1_vs_I2 <- FindMarkers(seurat, 
                            ident.1 = ident.1, 
                            ident.2 = ident.2, 
                            logfc.threshold = FCcutoff, 
                            test.use = methods[i])
  
    I1_vs_I2 <- I1_vs_I2[order(-I1_vs_I2$avg_log2FC), ]
    
    csv_name <- paste0(methods[i], sep="_",
                       ident.1, sep="_vs_", 
                       ident.2, sep="_", 
                       cluster, sep="", 
                       ".csv")
    
    write.csv(I1_vs_I2, 
              file.path(parent_dir, 
                        csv_name))
    
    results_list[[i]] <- I1_vs_I2
    
    #######################################################################
    
    plot_title <- paste0(ident.1, sep=" vs ", 
                         ident.2, sep=" - ", 
                         cluster, sep=" - ", 
                         methods[i])
    
    sub_title <- paste0("FC cutoff =", sep = " " , 
                        FCcutoff, sep= "     ", 
                        "p-value cutoff =", sep=" ", 
                        pCutoff)
    
    upper <- round(max(I1_vs_I2$avg_log2FC[which(I1_vs_I2$p_val <= pCutoff)])) + 1
    lower <- round(min(I1_vs_I2$avg_log2FC[which(I1_vs_I2$p_val <= pCutoff)])) + 1
    lim <- max(upper, abs(lower))
    
    I1_vs_I2_plot <- EnhancedVolcano(I1_vs_I2, 
                                        x="avg_log2FC", y="p_val", 
                                        lab= row.names(I1_vs_I2),
                                        pCutoff = pCutoff, FCcutoff = FCcutoff,
                                        title = plot_title,
                                        xlim = c(-lim,lim),
                                        subtitle = sub_title,
                                        caption = NULL)
    
    I1_vs_I2_plot <- I1_vs_I2_plot + theme(plot.subtitle=element_text(face="italic"))
    
    pdf_name <- paste0(methods[i],sep="_",
                       ident.1, sep="_vs_", 
                       ident.2, sep="_", 
                       cluster, sep="", 
                       ".pdf")
    
    pdf(file.path(parent_dir, pdf_name), height=6, width=8)
    print(I1_vs_I2_plot)
    dev.off()
    
    #######################################################################
    
    try(go.gsea(df = I1_vs_I2, 
                up_dir = up, 
                down_dir = down, 
                fileName = label, 
                pCutoff = pCutoff,
                db= "org.Mm.eg.db",
                dge_method = methods[i],
                color_by = color_by))
    
    try(kegg.gsea(df = I1_vs_I2, 
                  up_dir = up, 
                  down_dir = down, 
                  fileName = label, 
                  pCutoff = pCutoff,
                  dge_method = methods[i],
                  color_by = color_by))
    
    #######################################################################
    
    I1_vs_I2_up <- I1_vs_I2[which(I1_vs_I2$avg_log2FC >= FCcutoff & I1_vs_I2$p_val <= pCutoff), ]
    csv_name <- paste0(methods[i], sep="_",
                       cluster, sep="_",
                       "UP", sep="",
                       ".csv")
    write.csv(I1_vs_I2_up, file.path(up, csv_name))
    
    up_list[[i]] <- I1_vs_I2_up
    
    ora.enrich(df = I1_vs_I2_up, 
              dir = up,
              up_down = "UP", 
              label = label, 
              dge_method = methods[i], 
              uni_int = FALSE,
              color_by = color_by,
              pCutoff = pCutoff)
    
    # # # # # # # # # # # # # # # # # # # # # #
    
    I1_vs_I2_down <- I1_vs_I2[which(I1_vs_I2$avg_log2FC <= FCcutoff & I1_vs_I2$p_val <= pCutoff), ]
    I1_vs_I2_down <- I1_vs_I2_down[order(I1_vs_I2_down$avg_log2FC), ]
    
    csv_name <- paste0(methods[i], sep="_",
                       cluster, sep="_",
                       "DOWN", sep="",
                       ".csv")
    
    write.csv(I1_vs_I2_down, file.path(down, csv_name))
    
    down_list[[i]] <- I1_vs_I2_down
    
    ora.enrich(df = I1_vs_I2_down, 
               dir = down,
               up_down = "DOWN", 
               label = label, 
               dge_method = methods[i], 
               uni_int = FALSE,
               color_by = color_by,
               pCutoff = pCutoff)

  }
  
  ############################################################################
  
  union <- file.path(parent_dir, "3_Union")
  dir.create(union) 

  union_up <- file.path(union, "UP")
  dir.create(union_up)

  union_down <- file.path(union, "DOWN")
  dir.create(union_down)

  # # # # # # # # # # # # # # # # # # # # # #

  inter <- file.path(parent_dir, "4_Intersection")
  dir.create(inter)

  inter_up <- file.path(inter, "UP")
  dir.create(inter_up)

  inter_down <- file.path(inter, "DOWN")
  dir.create(inter_down)

  ##########################################################################

  up_union <- sort(unique(unlist(lapply(up_list, row.names))))
  up_union_ <- as.data.frame(up_union)
  colnames(up_union_) <- "GENES"

  write.xlsx(up_union_,
            file.path(union_up,
                      paste0("UP_union_gene_list", sep="_",
                             cluster, sep="",
                             ".xlsx")),
            rowNames = FALSE)
  rm(up_union_)

  # # # # # # # # # # # # # # # # # # # # # #

  down_union <- sort(unique(unlist(lapply(down_list, row.names))))
  down_union_ <- as.data.frame(down_union)
  colnames(down_union_) <- "GENES"

  write.xlsx(down_union_,
            file.path(union_down,
                      paste0("DOWN_union_gene_list", sep="_",
                             cluster, sep="",
                             ".xlsx")),
            rowNames = FALSE)
  rm(down_union_)

  # # # # # # # # # # # # # # # # # # # # # #

  up_intersection <- sort(Reduce(intersect, lapply(up_list, row.names)))
  up_intersection_ <- as.data.frame(up_intersection)
  colnames(up_intersection_) <- "GENES"

  write.xlsx(up_intersection_,
            file.path(inter_up,
                      paste0("UP_intersected_gene_list", sep="_",
                             cluster, sep="",
                             ".xlsx")),
            rowNames = FALSE)
  rm(up_intersection_)

  # # # # # # # # # # # # # # # # # # # # # #

  down_intersection <- sort(Reduce(intersect, lapply(down_list, row.names)))
  down_intersection_ <- as.data.frame(down_intersection)
  colnames(down_intersection_) <- "GENES"

  write.xlsx(down_intersection_,
            file.path(inter_down,
                      paste0("DOWN_intersected_gene_list", sep="_",
                             cluster, sep="",
                             ".xlsx")),
            rowNames = FALSE)
  rm(down_intersection_)

  ##########################################################################

  up_union_list <- vector(mode = "list", length = length(methods))
  names(up_union_list) <- methods

  down_union_list <- vector(mode = "list", length = length(methods))
  names(down_union_list) <- methods

  # # # # # # # # # # # # # # # # # # # # # #

  up_intersection_list <- vector(mode = "list", length = length(methods))
  names(up_intersection_list) <- methods

  down_intersection_list <- vector(mode = "list", length = length(methods))
  names(down_intersection_list) <- methods

  for (m in methods) {

    all_up <- results_list[[m]][which(results_list[[m]]$avg_log2FC > FCcutoff), ]
    in_df <- which(up_union %in% row.names(all_up))
    df <- all_up[up_union[in_df], ]
    df <- df[order(-df$avg_log2FC), ]
    up_union_list[[m]] <- df

    all_down <- results_list[[m]][which(results_list[[m]]$avg_log2FC < FCcutoff), ]
    in_df <- which(down_union %in% row.names(all_down))
    df <- all_down[down_union[in_df], ]
    df <- df[order(df$avg_log2FC), ]
    down_union_list[[m]] <- df

    # # # # # # # # # # # # # # # # # # # # # #

    df <- up_list[[m]][up_intersection,]
    df <- df[order(-df$avg_log2FC),]
    up_intersection_list[[m]] <- df

    df <- down_list[[m]][down_intersection,]
    df <- df[order(df$avg_log2FC),]
    down_intersection_list[[m]] <- df
    
  }
  
    ora.enrich(df = up_union_list[[1]], 
               dir = union_up,
               up_down = "UP", 
               label = paste0("union", sep="_",
                               cluster), 
               uni_int = TRUE,
               color_by = color_by,
               pCutoff = pCutoff)

    ora.enrich(df = down_union_list[[1]], 
               dir = union_down,
               up_down = "UP", 
               label = paste0("union", sep="_",
                               cluster), 
               uni_int = TRUE,
               color_by = color_by,
               pCutoff = pCutoff)
    
    ora.enrich(df = up_intersection_list[[1]], 
               dir = inter_up,
               up_down = "UP", 
               label = paste0("inter", sep="_",
                              cluster), 
               uni_int = TRUE,
               color_by = color_by,
               pCutoff = pCutoff)
    
    ora.enrich(df = down_intersection_list[[1]], 
               dir = inter_down,
               up_down = "DOWN", 
               label = paste0("inter", sep="_",
                              cluster), 
               uni_int = TRUE,
               color_by = color_by,
               pCutoff = pCutoff)

  ##########################################################################

  write.xlsx(up_union_list,
             file.path(union_up,
                       paste0("UP_union", sep="_",
                              cluster, sep="",
                              ".xlsx")), rowNames = TRUE)

  write.xlsx(down_union_list,
             file.path(union_down,
                       paste0("DOWN_union", sep="_",
                              cluster, sep="",
                              ".xlsx")), rowNames = TRUE)

  # # # # # # # # # # # # # # # # # # # # # #

  write.xlsx(up_intersection_list,
             file.path(inter_up,
                       paste0("UP_intersection", sep="_",
                              cluster, sep="",
                              ".xlsx")), rowNames = TRUE)

  write.xlsx(down_intersection_list,
             file.path(inter_down,
                       paste0("DOWN_intersection", sep="_",
                              cluster, sep="",
                              ".xlsx")), rowNames = TRUE)
  
}
```

``` r
myo_list <- list(atrial, ventriclular, SAN, OFT, AVC)
names(myo_list) <- c("Atrial","Ventriclular", "SAN", "OFT", "AVC")

for (i in 1:length(myo_list)) {
  dge.analysis(seurat = myo_list[[i]], 
               cluster = names(myo_list)[i], 
               methods = c("DESeq2","MAST","wilcox"),
               ident.1 = "KO", ident.2 = "WT",
               FCcutoff = 0.1, pCutoff = 0.05,
               path= dge_dir,
               color_by = "p.adjust")
}
```

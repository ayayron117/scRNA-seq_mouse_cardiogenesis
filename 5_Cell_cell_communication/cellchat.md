Cell-cell communication
================
Aaron Mohammed

``` r
library(Seurat)
library(CellChat)
library(patchwork)
library(ggalluvial) 
library(NMF)
library(reticulate)
library(ComplexHeatmap)
library(wordcloud)
use_condaenv("cellchat")


seurat_input_dir <- file.path(dirname(getwd()),
                              "2_Annotation_and_integration")

CellChatDB <- CellChatDB.mouse

future::plan("multisession", workers = 9) 

WT_dir <- file.path(getwd(), "WT")
dir.create(WT_dir)

KO_dir <- file.path(getwd(), "KO")
dir.create(KO_dir)

com_dir <- file.path(getwd(), "Comparisons")
dir.create(com_dir)
```

``` r
cardiac_s <- readRDS(file.path(seurat_input_dir, "final_cardiac_seurat.rds"))

cnames <- colnames(cardiac_s@meta.data)

`%nin%` <- Negate(`%in%`)

rem <- cnames[which(cnames %nin% c("orig.ident", "nCount_RNA", "nFeature_RNA", "pct_mito", "pct_ribo", "nuclear_fraction", "S.Score", "G2M.Score", "Phase", "Broad_Annotation", "Specific_Annotation"))]

cardiac_s@meta.data[ ,rem] <- NULL


Idents(cardiac_s) <- "orig.ident"

WT_s <- subset(cardiac_s, idents = "WT")
Idents(WT_s) <- "Broad_Annotation"
WT_s <- ScaleData(WT_s, row.names(WT_s))

KO_s <- subset(cardiac_s, idents = "KO")
Idents(KO_s) <- "Broad_Annotation"
KO_s <- ScaleData(KO_s, row.names(KO_s))

cellchat_WT <- createCellChat(object = WT_s, group.by = "Broad_Annotation", assay = "RNA")

cellchat_WT@DB <- CellChatDB

cellchat_WT <- subsetData(cellchat_WT) 

cellchat_WT <- identifyOverExpressedGenes(cellchat_WT)
cellchat_WT <- identifyOverExpressedInteractions(cellchat_WT)
cellchat_WT <- computeCommunProb(cellchat_WT, type = "triMean")
cellchat_WT <- filterCommunication(cellchat_WT, min.cells = 10)
df.net <- subsetCommunication(cellchat_WT, slot.name = "netP")
cellchat_WT <- computeCommunProbPathway(cellchat_WT)
cellchat_WT <- aggregateNet(cellchat_WT)




pdf(file.path(WT_dir, "1_WT_aggregated_network_plots.pdf"), height = 9, width = 11)

groupSize <- as.numeric(table(cellchat_WT@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_WT@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_WT@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat_WT@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pdf(file.path(WT_dir, "2_WT_TGFb_plots.pdf"), height = 5, width = 7)
netVisual_aggregate(cellchat_WT, signaling = "TGFb", layout = "chord")
netVisual_heatmap(cellchat_WT, signaling = "TGFb", color.heatmap = "Reds")
netVisual_aggregate(cellchat_WT, signaling = "TGFb", layout = "circle")
dev.off()

pdf(file.path(WT_dir, "3_WT_BMP_plots.pdf"), height = 5, width = 7)
netVisual_aggregate(cellchat_WT, signaling = "BMP", layout = "chord")
netVisual_heatmap(cellchat_WT, signaling = "BMP", color.heatmap = "Reds")
netVisual_aggregate(cellchat_WT, signaling = "BMP", layout = "circle")
dev.off()

pdf(file.path(WT_dir, "4_WT_NOTCH_plots.pdf"), height = 5, width = 7)
netVisual_aggregate(cellchat_WT, signaling = "NOTCH", layout = "chord")
netVisual_heatmap(cellchat_WT, signaling = "NOTCH", color.heatmap = "Reds")
netVisual_aggregate(cellchat_WT, signaling = "NOTCH", layout = "circle")
dev.off()

pdf(file.path(WT_dir, "5_WT_bubble_plot.pdf"), height = 30, width = 8)
netVisual_bubble(cellchat_WT, sources.use = 1, targets.use = c(2:5), remove.isolate = FALSE)
netVisual_bubble(cellchat_WT, sources.use = 2, targets.use = c(1,3:5), remove.isolate = FALSE)
netVisual_bubble(cellchat_WT, sources.use = 3, targets.use = c(1:2,4:5), remove.isolate = FALSE)
netVisual_bubble(cellchat_WT, sources.use = 4, targets.use = c(1:3,5), remove.isolate = FALSE)
netVisual_bubble(cellchat_WT, sources.use = 5, targets.use = c(1:5), remove.isolate = FALSE)
dev.off()

pdf(file.path(WT_dir, "6_WT_chord_outgoing.pdf"), height = 15, width = 15)
netVisual_chord_gene(cellchat_WT, sources.use = 1, targets.use = c(2:5), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_WT, sources.use = 2, targets.use = c(1,3:5), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_WT, sources.use = 3, targets.use = c(1:2,4:5), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_WT, sources.use = 4, targets.use = c(1:3,5), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_WT, sources.use = 5, targets.use = c(1:4), lab.cex = 0.5,legend.pos.y = 30)
dev.off()

pdf(file.path(WT_dir, "7_WT_network_centrality.pdf"), height = 5, width = 7)
cellchat_WT <- netAnalysis_computeCentrality(cellchat_WT, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat_WT, signaling = c("TGFb", "BMP", "NOTCH"), width = 8, height = 2.5, font.size = 10)
dev.off()

pdf(file.path(WT_dir, "8_WT_dominant_senders_receivers.pdf"), height = 5, width = 8)
netAnalysis_signalingRole_scatter(cellchat_WT)
dev.off()

pdf(file.path(WT_dir, "9_WT_signaling_roles.pdf"), height = 20, width = 15)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_WT, pattern = "outgoing", height = 20, width = 10)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_WT, pattern = "incoming", height = 20, width = 10)
ht1 + ht2
dev.off()


out_pat <- selectK(cellchat_WT, pattern = "outgoing")

pdf(file.path(WT_dir, "10_WT_out_sig_npatterns.pdf"), height = 5, width = 10)
out_pat
dev.off()

pdf(file.path(WT_dir, "11_WT_out_patterns.pdf"), height = 10, width = 15)
cellchat_WT <- identifyCommunicationPatterns(cellchat_WT, pattern = "outgoing", k = 4)
netAnalysis_river(cellchat_WT, pattern = "outgoing")
netAnalysis_dot(cellchat_WT, pattern = "outgoing")
dev.off()

inc_pat <- selectK(cellchat_WT, pattern = "incoming")

pdf(file.path(WT_dir, "12_WT_inc_sig_npatterns.pdf"), height = 5, width = 10)
inc_pat
dev.off()


pdf(file.path(WT_dir, "13_WT_inc_patterns.pdf"), height = 10, width = 15)
cellchat_WT <- identifyCommunicationPatterns(cellchat_WT, pattern = "incoming", k = 4)
netAnalysis_river(cellchat_WT, pattern = "incoming")
netAnalysis_dot(cellchat_WT, pattern = "incoming")
dev.off()


cellchat_WT <- computeNetSimilarity(cellchat_WT, type = "functional")
cellchat_WT <- netEmbedding(cellchat_WT, type = "functional")
cellchat_WT <- netClustering(cellchat_WT, type = "functional")

pdf(file.path(WT_dir, "14_WT_functional_similarity.pdf"), height = 10, width = 15)
netVisual_embedding(cellchat_WT, type = "functional", label.size = 3.5)
dev.off()

cellchat_WT <- computeNetSimilarity(cellchat_WT, type = "structural")
cellchat_WT <- netEmbedding(cellchat_WT, type = "structural")
cellchat_WT <- netClustering(cellchat_WT, type = "structural")

pdf(file.path(WT_dir, "15_WT_structural_similarity.pdf"), height = 10, width = 15)
netVisual_embedding(cellchat_WT, type = "structural", label.size = 3.5)
dev.off()

# saveRDS(cellchat_WT, file.path(getwd(),"WT_cellchat.rds"))
```

``` r
cellchat_KO <- createCellChat(object = KO_s, group.by = "Broad_Annotation", assay = "RNA")

cellchat_KO@DB <- CellChatDB

cellchat_KO <- subsetData(cellchat_KO) # This step is necessary even if using the whole database
cellchat_KO <- identifyOverExpressedGenes(cellchat_KO)
cellchat_KO <- identifyOverExpressedInteractions(cellchat_KO)

cellchat_KO <- computeCommunProb(cellchat_KO, type = "triMean")
cellchat_KO <- filterCommunication(cellchat_KO, min.cells = 10)
df.net <- subsetCommunication(cellchat_KO, slot.name = "netP")
cellchat_KO <- computeCommunProbPathway(cellchat_KO)
cellchat_KO <- aggregateNet(cellchat_KO)

pdf(file.path(KO_dir, "1_KO_aggregated_network_plots.pdf"), height = 9, width = 11)
groupSize <- as.numeric(table(cellchat_KO@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_KO@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_KO@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat_KO@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pdf(file.path(KO_dir, "2_KO_TGFb_plots.pdf"), height = 5, width = 7)
netVisual_aggregate(cellchat_KO, signaling = "TGFb", layout = "chord")
netVisual_heatmap(cellchat_KO, signaling = "TGFb", color.heatmap = "Reds")
netVisual_aggregate(cellchat_KO, signaling = "TGFb", layout = "circle")
dev.off()

pdf(file.path(KO_dir, "3_KO_BMP_plots.pdf"), height = 5, width = 7)
netVisual_aggregate(cellchat_KO, signaling = "BMP", layout = "chord")
netVisual_heatmap(cellchat_KO, signaling = "BMP", color.heatmap = "Reds")
netVisual_aggregate(cellchat_KO, signaling = "BMP", layout = "circle")
dev.off()

pdf(file.path(KO_dir, "4_KO_NOTCH_plots.pdf"), height = 5, width = 7)
netVisual_aggregate(cellchat_KO, signaling = "NOTCH", layout = "chord")
netVisual_heatmap(cellchat_KO, signaling = "NOTCH", color.heatmap = "Reds")
netVisual_aggregate(cellchat_KO, signaling = "NOTCH", layout = "circle")
dev.off()

pdf(file.path(KO_dir, "5_KO_bubble_plot.pdf"), height = 30, width = 8)
netVisual_bubble(cellchat_KO, sources.use = 1, targets.use = c(2:5), remove.isolate = FALSE)
netVisual_bubble(cellchat_KO, sources.use = 2, targets.use = c(1,3:5), remove.isolate = FALSE)
netVisual_bubble(cellchat_KO, sources.use = 3, targets.use = c(1:2,4:5), remove.isolate = FALSE)
netVisual_bubble(cellchat_KO, sources.use = 4, targets.use = c(1:3,5), remove.isolate = FALSE)
netVisual_bubble(cellchat_KO, sources.use = 5, targets.use = c(1:5), remove.isolate = FALSE)
dev.off()

pdf(file.path(KO_dir, "6_KO_chord_outgoing.pdf"), height = 15, width = 15)
netVisual_chord_gene(cellchat_KO, sources.use = 1, targets.use = c(2:5), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_KO, sources.use = 2, targets.use = c(1,3:5), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_KO, sources.use = 3, targets.use = c(1:2,4:5), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_KO, sources.use = 4, targets.use = c(1:3,5), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_KO, sources.use = 5, targets.use = c(1:4), lab.cex = 0.5,legend.pos.y = 30)
dev.off()

pdf(file.path(KO_dir, "7_KO_network_centrality.pdf"), height = 5, width = 7)
cellchat_KO <- netAnalysis_computeCentrality(cellchat_KO, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat_KO, signaling = c("TGFb", "BMP", "NOTCH"), width = 8, height = 2.5, font.size = 10)
dev.off()

pdf(file.path(KO_dir, "8_KO_dominant_senders_receivers.pdf"), height = 5, width = 8)
netAnalysis_signalingRole_scatter(cellchat_KO)
dev.off()

pdf(file.path(KO_dir, "9_KO_signaling_roles.pdf"), height = 20, width = 15)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_KO, pattern = "outgoing", height = 20, width = 10)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_KO, pattern = "incoming", height = 20, width = 10)
ht1 + ht2
dev.off()


out_pat <- selectK(cellchat_KO, pattern = "outgoing")

pdf(file.path(KO_dir, "10_KO_out_sig_npatterns.pdf"), height = 5, width = 10)
out_pat
dev.off()

pdf(file.path(KO_dir, "11_KO_out_patterns.pdf"), height = 10, width = 15)
cellchat_KO <- identifyCommunicationPatterns(cellchat_KO, pattern = "outgoing", k = 4)
netAnalysis_river(cellchat_KO, pattern = "outgoing")
netAnalysis_dot(cellchat_KO, pattern = "outgoing")
dev.off()

inc_pat <- selectK(cellchat_KO, pattern = "incoming")

pdf(file.path(KO_dir, "12_KO_inc_sig_npatterns.pdf"), height = 5, width = 10)
inc_pat
dev.off()


pdf(file.path(KO_dir, "13_KO_inc_patterns.pdf"), height = 10, width = 15)
cellchat_KO <- identifyCommunicationPatterns(cellchat_KO, pattern = "incoming", k = 4)
netAnalysis_river(cellchat_KO, pattern = "incoming")
netAnalysis_dot(cellchat_KO, pattern = "incoming")
dev.off()


cellchat_KO <- computeNetSimilarity(cellchat_KO, type = "functional")
cellchat_KO <- netEmbedding(cellchat_KO, type = "functional")
cellchat_KO <- netClustering(cellchat_KO, type = "functional")

pdf(file.path(KO_dir, "14_KO_functional_similarity.pdf"), height = 10, width = 15)
netVisual_embedding(cellchat_KO, type = "functional", label.size = 3.5)
dev.off()

cellchat_KO <- computeNetSimilarity(cellchat_KO, type = "structural")
cellchat_KO <- netEmbedding(cellchat_KO, type = "structural")
cellchat_KO <- netClustering(cellchat_KO, type = "structural")

pdf(file.path(KO_dir, "15_KO_structural_similarity.pdf"), height = 10, width = 15)
netVisual_embedding(cellchat_KO, type = "structural", label.size = 3.5)
dev.off()

# saveRDS(cellchat_KO, file.path(getwd(),"KO_cellchat.rds"))
```

``` r
object.list <- list(WT = cellchat_WT, KO = cellchat_KO)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB

pdf(file.path(com_dir, "1_interaction_comparison.pdf"), height = 10, width = 15)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

pdf(file.path(com_dir, "2_differential_interactions.pdf"), height = 10, width = 15)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

pdf(file.path(com_dir, "3_differential_heatmap.pdf"), height = 10, width = 15)
gg1 <- netVisual_heatmap(cellchat, color.heatmap = c("blue","darkred"))
gg2 <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = c("blue","darkred"))
gg1 + gg2
dev.off()



group.cellType <- c("Myocardial", "Epicardial", "Endocardial", "MP", "Mesenchymal")
group.cellType <- factor(group.cellType, levels = c("Myocardial", "Epicardial", "Endocardial", "MP", "Mesenchymal"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

pdf(file.path(com_dir, "4_num_interactions.pdf"), height = 10, width = 15)
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()



pdf(file.path(com_dir, "5_diff_num_interac.pdf"), height = 10, width = 15)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)
dev.off()


pdf(file.path(com_dir, "6_major_sources.pdf"), height = 10, width = 15)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
dev.off()


pdf(file.path(com_dir, "7_signaling_changes.pdf"), height = 10, width = 15)
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Myocardial")
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Epicardial")
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Endocardial")
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "MP")
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mesenchymal")
dev.off()

pdf(file.path(com_dir, "8_functional_similarity.pdf"), height = 10, width = 15)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()

pdf(file.path(com_dir, "9_structural_similarity.pdf"), height = 10, width = 15)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()

pdf(file.path(com_dir, "10_functional_sim_ranked.pdf"), height = 10, width = 15)
rankSimilarity(cellchat, type = "functional")
dev.off()

pdf(file.path(com_dir, "11_structural_sim_ranked.pdf"), height = 10, width = 15)
rankSimilarity(cellchat, type = "structural")
dev.off()

pdf(file.path(com_dir, "12_information_flow.pdf"), height = 20, width = 10)
rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
dev.off()


pdf(file.path(com_dir, "13_outgoing_patterns.pdf"), height = 20, width = 15)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], height = 30, width = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], height = 30, width = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf(file.path(com_dir, "14_incoming_patterns.pdf"), height = 20, width = 15)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], height = 30, width = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], height = 30, width = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



pdf(file.path(com_dir, "15_communication_probabilities.pdf"), height = 20, width = 10)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:5),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,3:5),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2,4:5),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:3,5),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:4),  comparison = c(1, 2), angle.x = 45)
dev.off()


gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:5),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:5),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T)

gg3 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,3:5),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T)
gg4 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,3:5),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T)

gg5 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2,4:5),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T)
gg6 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2,4:5),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T)

gg7 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:3,5),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T)
gg8 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:3,5),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T)

gg9 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:4),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T)
gg10 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:4),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T)

pdf(file.path(com_dir, "16_KO_differential_signaling.pdf"), height = 20, width = 10)
gg1 + gg2
gg3 + gg4
gg5 + gg6
gg7 + gg8
gg9 + gg10
dev.off()


pos.dataset = "KO"
features.name = paste0(pos.dataset, ".merged")
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "KO",ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "WT",ligand.logFC = -0.05, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)



pairLR.use.up = net.up[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 1, targets.use = c(2:5), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 1, targets.use = c(2:5), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


gg3 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 2, targets.use = c(1,3:5), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg4 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 2, targets.use = c(1,3:5), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


gg5 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 3, targets.use = c(1:2,4:5), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg6 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 3, targets.use = c(1:2,4:5), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


gg7 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(1:3,5), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg8 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(1:3,5), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

gg9 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 5, targets.use = c(1:4), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg10 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 5, targets.use = c(1:4), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


pdf(file.path(com_dir, "17_KO_differential_L-R.pdf"), height = 20, width = 10)
gg1 + gg2
gg3 + gg4
gg5 + gg6
gg7 + gg8
gg9 + gg10
dev.off()




pdf(file.path(com_dir, "18_KO_differential_L-R_chord.pdf"), height = 13, width = 13)
netVisual_chord_gene(object.list[[2]], sources.use = 1, targets.use = c(2:5), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 1, targets.use = c(2:5), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 2, targets.use = c(1,3:5), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 2, targets.use = c(1,3:5), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 3, targets.use = c(1:2,4:5), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 3, targets.use = c(1:2,4:5), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(1:3,5), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(1:3,5), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = 5, targets.use = c(1:4), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 5, targets.use = c(1:4), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

pdf(file.path(com_dir, "19_KO_downreg_L-R.pdf"), height = 6, width = 6)
computeEnrichmentScore(net.down, species = 'mouse', variable.both = TRUE)
dev.off()

pdf(file.path(com_dir, "20_KO_upreg_L-R.pdf"), height = 6, width = 6)
computeEnrichmentScore(net.up, species = 'mouse', variable.both = TRUE)
dev.off()

# saveRDS(cellchat, file.path(getwd(),"merged_cellchat.rds"))
```

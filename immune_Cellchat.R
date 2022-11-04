# CellChat
library(dplyr)
library(Seurat)
library(future)
library(sctransform)
library(patchwork)
library(ggplot2)
library(dittoSeq)
library(tidyverse)
library(SingleR)
library(SingleCellExperiment)
library(celldex)
library(patchwork)
library(fgsea)
library(monocle3)
library(scRNAseq)
library(fgsea)
library(infercnv)
library(viridis)
library(slingshot)
library("scProportionTest")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
setwd("D:/TAK981_KPC/control_vs_TAK981")
load("immune_Cellchat.RData")

DefaultAssay(immune) <- "RNA"
immune <- NormalizeData(immune)
all.genes <- rownames(immune)
immune <- ScaleData(immune, features = all.genes,verbose = FALSE)
immune <- FindVariableFeatures(immune, selection.method = "vst", nfeatures = 2500,verbose = FALSE)
immune <- RunPCA(immune, verbose = FALSE)
immune <- RunUMAP(immune, dims = 1:40,verbose = FALSE)
immune <- FindNeighbors(immune, dims = 1:2, reduction = "umap", verbose = FALSE)
immune <- FindClusters(immune, resolution = 1)
immune[["SingleR.label"]]<- replace(immune[["SingleR.label"]], immune[["SingleR.label"]]=="cDC1-Ccl22", "DC")
immune[["SingleR.label"]]<- replace(immune[["SingleR.label"]], immune[["SingleR.label"]]=="cDC1-Clec9a", "DC")
immune[["SingleR.label"]]<- replace(immune[["SingleR.label"]], immune[["SingleR.label"]]=="cDC2-Itgax", "DC")
immune[["SingleR.label"]]<- replace(immune[["SingleR.label"]], immune[["SingleR.label"]]=="pDC", "DC")
immune[["SingleR.label"]]<- replace(immune[["SingleR.label"]], immune[["SingleR.label"]]=="PMN-MDSC", "Granulocytes")
DimPlot(immune, reduction="umap",label=TRUE,pt.size=1, group.by = "SingleR.label", repel=T)

# Split for comparison analysis
treatment_split <- SplitObject(immune, split.by = "treatment")
control <- treatment_split$control
TAK981 <- treatment_split$TAK981

# control
control_cellchat <- GetAssayData(control, assay = "RNA", slot = "data") # normalized data matrix
Idents(control) <- "SingleR.label"
labels <- Idents(control)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
control_cellchat <- createCellChat(object = control_cellchat, meta = meta, group.by = "group") 
control_cellchat <- addMeta(control_cellchat, meta = meta, meta.name = "labels")
control_cellchat <- setIdent(control_cellchat, ident.use = "labels")
levels(control_cellchat@idents) 
groupSize <- as.numeric(table(control_cellchat@idents))
control_cellchatDB <- CellChatDB.mouse 
control_cellchat@DB <- control_cellchatDB
control_cellchat <- subsetData(control_cellchat)
control_cellchat <- identifyOverExpressedGenes(control_cellchat)
control_cellchat <- identifyOverExpressedInteractions(control_cellchat)
control_cellchat <- projectData(control_cellchat, PPI.mouse)
control_cellchat <- computeCommunProb(control_cellchat)
# control_cellchat <- filterCommunication(control_cellchat, min.cells = 20)
control_cellchat <- computeCommunProbPathway(control_cellchat)
control_cellchat <- aggregateNet(control_cellchat)
groupSize <- as.numeric(table(control_cellchat@idents))
par(mfrow = c(2,2), xpd=TRUE)
netVisual_circle(control_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(control_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# TAK981
TAK981_cellchat <- GetAssayData(TAK981, assay = "RNA", slot = "data") # normalized data matrix
Idents(TAK981) <- "SingleR.label"
labels <- Idents(TAK981)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
TAK981_cellchat <- createCellChat(object = TAK981_cellchat, meta = meta, group.by = "group")
TAK981_cellchat <- addMeta(TAK981_cellchat, meta = meta, meta.name = "labels")
TAK981_cellchat <- setIdent(TAK981_cellchat, ident.use = "labels")
levels(TAK981_cellchat@idents) 
groupSize <- as.numeric(table(TAK981_cellchat@idents)) 
TAK981_cellchatDB <- CellChatDB.mouse 
TAK981_cellchat@DB <- TAK981_cellchatDB
TAK981_cellchat <- subsetData(TAK981_cellchat)
TAK981_cellchat <- identifyOverExpressedGenes(TAK981_cellchat)
TAK981_cellchat <- identifyOverExpressedInteractions(TAK981_cellchat)
TAK981_cellchat <- projectData(TAK981_cellchat, PPI.mouse)
TAK981_cellchat <- computeCommunProb(TAK981_cellchat)
# TAK981_cellchat <- filterCommunication(TAK981_cellchat, min.cells = 5)
TAK981_cellchat <- computeCommunProbPathway(TAK981_cellchat)
TAK981_cellchat <- aggregateNet(TAK981_cellchat) 
groupSize <- as.numeric(table(TAK981_cellchat@idents))
par(mfrow = c(2,2), xpd=TRUE)
netVisual_circle(TAK981_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(TAK981_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


# Comparison Analysis
group.new = levels(control_cellchat@idents)
TAK981_cellchat <- liftCellChat(TAK981_cellchat, group.new)
object.list <- list(control = control_cellchat, TAK981 = TAK981_cellchat)
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}
combined_cellchat <- mergeCellChat(object.list, add.names = names(object.list))
combined_cellchat
# Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(combined_cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(combined_cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(2,2), xpd=TRUE)
netVisual_diffInteraction(combined_cellchat, weight.scale = T, measure = "weight", top=1)

netVisual_diffInteraction(combined_cellchat, weight.scale = T, measure = "weight", top=1,
                          sources.use = c("T/NK cells"), targets.use = c("T/NK cells"))

tiff("interaction_number_web.jpeg", unit="in", width=6, height=6, res=500)
netVisual_diffInteraction(combined_cellchat, weight.scale = T, measure = "count")
dev.off()

tiff("interaction_strength_web.jpeg", unit="in", width=6, height=6, res=500)
netVisual_diffInteraction(combined_cellchat, weight.scale = T, measure = "weight")
dev.off()


gg1 <- netVisual_heatmap(combined_cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(combined_cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

netVisual_heatmap(combined_cellchat, measure = "weight", 
                  row.show = c("iCAF", "myCAF", "apCAF", "Ductal 1", "Ductal 2", "Endothelial"),
                  col.show = c("iCAF", "myCAF", "apCAF", "Ductal 1", "Ductal 2", "Endothelial"))

# Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

tiff("interaction_weight_comparison.jpeg", unit="in", width=8, height=4, res=500)
patchwork::wrap_plots(plots = gg)
dev.off()

control_cellchat@netP$pathways
TAK981_cellchat@netP$pathways

# Hierarchy plot
# Hierarchy plot
pathways.show <- c("CD80")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
vertex.receiver = seq(1,10) # Left portion of hierarchy plot the shows signaling to dermal cells and right portion shows signaling to epidermal cells
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver,layout = "hierarchy",
                      edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("CXCL")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
vertex.receiver = seq(1,10) # Left portion of hierarchy plot the shows signaling to dermal cells and right portion shows signaling to epidermal cells
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, 
                      edge.weight.max = weight.max[1], edge.width.max = 10, layout = "hierarchy",
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

gg1 <- netAnalysis_signalingChanges_scatter(combined_cellchat, idents.use = "T/NK cells")
gg2 <- netAnalysis_signalingChanges_scatter(combined_cellchat, idents.use = "Mono/MDSC")
patchwork::wrap_plots(plots = list(gg1,gg2))

#Part II: Identify the conserved and context-specific signaling pathways
#functional
combined_cellchat <- computeNetSimilarityPairwise(combined_cellchat, type = "functional")
combined_cellchat <- netEmbedding(combined_cellchat, type = "functional", umap.method = "uwot")
combined_cellchat <- netClustering(combined_cellchat, type = "functional")
netVisual_embeddingPairwise(combined_cellchat, type = "functional", label.size = 3.5)

#structural
combined_cellchat <- computeNetSimilarityPairwise(combined_cellchat, type = "structural")
combined_cellchat <- netEmbedding(combined_cellchat, type = "structural",umap.method = "uwot")
combined_cellchat <- netClustering(combined_cellchat, type = "structural")
netVisual_embeddingPairwise(combined_cellchat, type = "structural", label.size = 2)

# Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(combined_cellchat, type = "functional")

#Identify and visualize the conserved and context-specific signaling pathways
gg1 <- rankNet(combined_cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(combined_cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

unique(combined_cellchat@idents$joint)

#outgoing signal
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#incoming signal
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#Overal
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
#Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(combined_cellchat, thresh=0.05, sources.use = 1, targets.use = c(15,16,17,18,19,20,21,22,23),  comparison = c(1,2), angle.x = 45)

gg1 <- netVisual_bubble(combined_cellchat, sources.use = 5, targets.use = c(2,15,17,18,19,20,21,22,23),  comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(combined_cellchat, sources.use = 5, targets.use = c(2,15,17,18,19,20,21,22,23),  comparison = c(1,2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2



#Identify dysfunctional signaling by using differential expression analysis
pos.dataset = "TAK981"
features.name = pos.dataset
combined_cellchat <- identifyOverExpressedGenes(combined_cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(combined_cellchat, features.name = features.name)
net.up <- subsetCommunication(combined_cellchat, net = net, datasets = "TAK981",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(combined_cellchat, net = net, datasets = "control",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, combined_cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, combined_cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(combined_cellchat, pairLR.use = pairLR.use.up, sources.use = 1, targets.use = c(15,16,17,18,19,20,21,22,23), comparison = c(1, 2), 
                        angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(combined_cellchat, pairLR.use = pairLR.use.down, sources.use = 1, targets.use = c(15,16,17,18,19,20,21,22,23), comparison = c(1, 2),
                        angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 1, targets.use = c(15:23), slot.name = 'net', net = net.up, 
                     lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 1, targets.use = c(15:23), slot.name = 'net', net = net.down, 
                     lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))



# Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("CXCL")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))


# Chord diagram
pathways.show <- c("CD80") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

# NET Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(object.list[[1]]@idents)
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}


unique(combined_cellchat@idents$joint)


#Net chord diagram - gene
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(18), targets.use = c(2,5), lab.cex = 0.5, title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
}

# compare all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(18), targets.use = c(2,5),  title.name = paste0("Signaling received by Inflam.DC and .TC - ", names(object.list)[i]), legend.pos.x = 10)
}

# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(18), targets.use = c(2,5),slot.name = "netP", title.name = paste0("Signaling pathways sending from fibroblast - ", names(object.list)[i]), legend.pos.x = 10)
}

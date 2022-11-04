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
setwd("D:/control_vs_TAK981")
load("cellchat_new.RData")
setwd("D:/control_vs_TAK981/cellchat 07192022")

DefaultAssay(immune) <- "RNA"
DefaultAssay(nonimmune) <- "RNA"
scKPC <- merge(immune, y=c(nonimmune))
scKPC <- DietSeurat(scKPC, assay="RNA")

DefaultAssay(scKPC) <- "RNA"
scKPC <- NormalizeData(scKPC)
all.genes <- rownames(scKPC)
scKPC <- ScaleData(scKPC, features = all.genes,verbose = FALSE)
scKPC <- FindVariableFeatures(scKPC, selection.method = "vst", nfeatures = 2500,verbose = FALSE)
scKPC <- RunPCA(scKPC, verbose = FALSE)
scKPC <- RunUMAP(scKPC, dims = 1:40,verbose = FALSE)
scKPC <- FindNeighbors(scKPC, dims = 1:2, reduction = "umap", verbose = FALSE)
scKPC <- FindClusters(scKPC, resolution = 1)
DimPlot(scKPC, reduction="umap",label=TRUE,pt.size=1, group.by = "SingleR.label", repel=T)

# scKPC.markers <- FindAllMarkers(scKPC, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# scKPC.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC) -> scKPC.top10
# write.csv(scKPC.markers, "scKPC_markers.csv", row.names = T)
# jpeg("scKPC_cluster_heatmap.jpeg", width=3000, height=2000)
# DoHeatmap(scKPC, features = scKPC.top10$gene) + NoLegend()
# dev.off()
scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="cDC1-Ccl22", "DC")
scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="cDC1-Clec9a", "DC")
scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="cDC2-Itgax", "DC")
scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="pDC", "DC")
scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="PMN-MDSC", "Granulocytes")
scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="Macro-Isg15", "Mono-Ly6c")
scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="Macro-C1q", "Macro-MHCII")

scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="Macro-MHCII", "Macrophage")
scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="Macro-Proliferating", "Macrophage")
scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="Macro-Spp1", "Macrophage")

scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="myCAF", "CAF")
scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="iCAF", "CAF")
scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="apCAF", "CAF")
scKPC[["SingleR.label"]]<- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="PSC", "CAF")

scKPC <- subset(scKPC, subset=(SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2" |
                                 SingleR.label=="Endothelial" | SingleR.label=="CAF" |
                                 SingleR.label=="Mono/MDSC" | SingleR.label=="Mono-Ly6c" | SingleR.label=="Macrophage" | SingleR.label=="DC"|
                                 SingleR.label=="T/NK cells" | SingleR.label=="Plasma cells" | SingleR.label=="Granulocytes" |
                                 SingleR.label=="B cells"))

DimPlot(scKPC, reduction="umap",label=TRUE,pt.size=1, group.by = "SingleR.label", repel = T)


####################################################################################################################

# Part I: Data input & processing and initialization of CellChat object
scKPC_cellchat <- GetAssayData(scKPC, assay = "RNA", slot = "data") # normalized data matrix
Idents(scKPC) <- "SingleR.label"
labels <- Idents(scKPC)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat <- createCellChat(object = scKPC_cellchat, meta = meta, group.by = "group")
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

#Set ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

#Preprocessing the expression data for cell-cell communication analysis
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
# future::plan("multiprocess", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)


# Part II: Inference of cell-cell communication network
# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 35)

# Extract the inferred cellular communication network as a data frame
# Option 1
# returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
# # Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# df.net <- subsetCommunication(cellchat) 
# 
# # Option 2
# # gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
# 
# # Option 3
# # gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) 
# 
# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
# USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(2,3), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# tiff("interaction weights strength.tiff", unit="in", width=7, height=6,res=500)
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# dev.off()

# examine the signaling sent from each cell group
mat <- cellchat@net$weight
par(mfrow = c(5,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Part III: Visualization of cell-cell communication network
#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
pathways.list <- as.data.frame(unique(CellChatDB[["interaction"]][["pathway_name"]]))
unique(CellChatDB[["interaction"]][["pathway_name"]])
pathways.show <- c("TGFb") #VCAM, TBFb, VEGF, PD-L1, MHCI, 

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver =seq(1,7)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "hierarchy")

# Circle plot
pathways.show <- c("SPP1") #VCAM, TGFb, VEGF, PD-L1, MHCI, CD80, CD86, SELL
par(mfrow=c(2,2))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

tiff("SPP1.tiff", unit="in", width=7, height=6,res=500)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")


# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
# names(group.cellType) <- levels(cellchat@idents)
# netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
# #Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
# netAnalysis_contribution(cellchat, signaling = pathways.show)
# pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
# LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# # Hierarchy plot
# vertex.receiver = seq(1,4) # a numeric vector
# netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout="hierarchy")

pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = c(1,2,3,4,5,6,7,8,9,10)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

######################################################################################################
# Split for comparison analysis
treatment_split <- SplitObject(scKPC, split.by = "treatment")
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
control_cellchat <- filterCommunication(control_cellchat, min.cells = 20)
control_cellchat <- computeCommunProbPathway(control_cellchat)
control_cellchat <- aggregateNet(control_cellchat)
groupSize <- as.numeric(table(control_cellchat@idents))
par(mfrow = c(2,2), xpd=TRUE)
netVisual_circle(control_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(control_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

setwd("D:/control_vs_TAK981/cellchat 07192022/control")
pathways.show.all <- control_cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(control_cellchat@idents)
vertex.receiver = c(7,8,9,10)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(control_cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(control_cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
# tiff("control_interaction_weights_strength.tiff", unit="in", width=7, height=6,res=500)
# netVisual_circle(control_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# dev.off()


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
TAK981_cellchat <- filterCommunication(TAK981_cellchat, min.cells = 5)
TAK981_cellchat <- computeCommunProbPathway(TAK981_cellchat)
TAK981_cellchat <- aggregateNet(TAK981_cellchat) 
groupSize <- as.numeric(table(TAK981_cellchat@idents))
par(mfrow = c(2,2), xpd=TRUE)
netVisual_circle(TAK981_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(TAK981_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# tiff("TAK981_interaction_weights_strength.tiff", unit="in", width=7, height=6,res=500)
# netVisual_circle(TAK981_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# dev.off()
setwd("D:/control_vs_TAK981/cellchat 07192022/TAK981")
pathways.show.all <- TAK981_cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(TAK981_cellchat@idents)
vertex.receiver = c(7,8,9,10)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(TAK981_cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(TAK981_cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}


##########################################

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
# netVisual_diffInteraction(combined_cellchat, weight.scale = T)
netVisual_diffInteraction(combined_cellchat, weight.scale = T, measure = "weight")

tiff("diffInteraction.tiff", unit="in", width=7, height=6,res=500)
netVisual_diffInteraction(combined_cellchat, weight.scale = T, measure = "weight")
dev.off()

netVisual_diffInteraction(combined_cellchat, weight.scale = T, measure = "weight",
                          sources.use = c("CAF", "Ductal 1", "Ductal 2", "Mono-Ly6c", "Mono/MDSC", "DC", "T/NK cells", "Macrophage"),
                          targets.use = c("CAF", "Ductal 1", "Ductal 2", "Mono-Ly6c", "Mono/MDSC", "DC", "T/NK cells", "Macrophage"))

gg1 <- netVisual_heatmap(combined_cellchat)
gg2 <- netVisual_heatmap(combined_cellchat, measure = "weight")
gg1 + gg2

tiff("diffInteraction_heatmap_no_endo.tiff", unit="in", width=7, height=6,res=500)
netVisual_heatmap(combined_cellchat, measure = "weight", 
                  row.show = c("CAF", "Ductal 1", "Ductal 2", "Mono-Ly6c", "Mono/MDSC", "T/NK cells", "Macrophage"),
                  col.show = c("CAF", "Ductal 1", "Ductal 2", "Mono-Ly6c", "Mono/MDSC", "T/NK cells", "Macrophage"))
dev.off()


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

tiff("interaction_weight_comparison.jpeg", unit="in", width=8, height=4.5, res=500)
patchwork::wrap_plots(plots = gg)
dev.off()

gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax,
                                               x.measure="outdeg_unweighted", y.measure = "indeg_unweighted")
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
tiff("interaction_number_comparison.jpeg", unit="in", width=8, height=4.5, res=500)
patchwork::wrap_plots(plots = gg)
dev.off()


# Hierarchy plot
pathways.show <- c("COLLAGEN")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
vertex.receiver = seq(1,10) # Left portion of hierarchy plot the shows signaling to dermal cells and right portion shows signaling to epidermal cells
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}


gg1 <- netAnalysis_signalingChanges_scatter(combined_cellchat, idents.use = "CAF")
gg2 <- netAnalysis_signalingChanges_scatter(combined_cellchat, idents.use = "Mono/MDSC")
gg3 <- netAnalysis_signalingChanges_scatter(combined_cellchat, idents.use = "Endothelial")
gg4 <- netAnalysis_signalingChanges_scatter(combined_cellchat, idents.use = "T/NK cells")
patchwork::wrap_plots(plots = list(gg1,gg2,gg3, gg4))



#Part II: Identify the conserved and context-specific signaling pathways
#functional
combined_cellchat <- computeNetSimilarityPairwise(combined_cellchat, type = "functional")
combined_cellchat <- netEmbedding(combined_cellchat, type = "functional", umap.method = "uwot")
combined_cellchat <- netClustering(combined_cellchat, type = "functional")
netVisual_embeddingPairwise(combined_cellchat, type = "functional", label.size=3, do.label=T)

 #structural ***** WARNING: TAKES TOO FUCKING LONG
# combined_cellchat <- computeNetSimilarityPairwise(combined_cellchat, type = "structural")
# combined_cellchat <- netEmbedding(combined_cellchat, type = "structural",umap.method = "uwot")
# combined_cellchat <- netClustering(combined_cellchat, type = "structural")
# netVisual_embeddingPairwise(combined_cellchat, type = "structural", label.size = 2)

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
gg1 <- netVisual_bubble(combined_cellchat, pairLR.use = pairLR.use.up, sources.use = 8, targets.use = c(7,8,9,10), comparison = c(1, 2), 
                        angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(combined_cellchat, pairLR.use = pairLR.use.down, sources.use = 8, targets.use = c(7,8,9,10), comparison = c(1, 2),
                        angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

tiff("interaction_weight_comparison.jpeg", unit="in", width=8, height=4.5, res=500)
patchwork::wrap_plots(plots = gg)
dev.off()


par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 8, targets.use = c(7:10), slot.name = 'net', net = net.up, 
                     lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 10, targets.use = c(7:10), slot.name = 'net', net = net.down, 
                     lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))



# Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("SPP1")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

tiff("SPP1_differential.jpeg", unit="in", width=10, height=6, res=500)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


pathways.show <- c("TGFb") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))


# Chord diagram
pathways.show <- c("FN1") 
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

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
load("immune.RData")
load("nonimmune.RData")

memory.limit(56000)
scKPC <- merge(immune, y=c(nonimmune))
scKPC <- DietSeurat(scKPC, assay="RNA")
scKPC <- NormalizeData(scKPC, normalization.method = "LogNormalize", scale.factor = 10000)
scKPC <- FindVariableFeatures(scKPC, selection.method = "vst", nfeatures = 2000)
features <- rownames(scKPC)
scKPC <- ScaleData(scKPC, features = features, vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
scKPC <- RunPCA(scKPC, features=features)
ElbowPlot(object = scKPC, ndims = 50)
scKPC <- RunTSNE(scKPC, dims = 1:30, verbose = FALSE)
scKPC <- RunUMAP(scKPC, dims = 1:30, verbose = FALSE)
scKPC <- FindNeighbors(scKPC, dims = 1:30,verbose = FALSE)
scKPC <- FindClusters(scKPC, resolution = 0.2)

DimPlot(scKPC, reduction="umap",label=TRUE,pt.size=1, group.by = "SingleR.label", repel = T)
temp <- scKPC
scKPC <- subset(scKPC, subset=(SingleR.label=="Acinar"| SingleR.label=="Ductal 0" | SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2" | 
                                 SingleR.label=="Ductal NM" | SingleR.label=="Endothelial" | SingleR.label=="Fibroblasts" | SingleR.label== "ADM" |
                                 SingleR.label== "iCAF" | SingleR.label== "myCAF" | SingleR.label== "myCAF-like PSC" | SingleR.label== "Endothelial" |
                                 SingleR.label=="Macro-Spp1" |SingleR.label=="Macro-Ly6c+Chil3+" | SingleR.label=="Macro-C1q" | 
                                 SingleR.label=="Macro-Ly6c+Isg15+" |SingleR.label=="Macro-Proliferating" | SingleR.label=="pDC" | 
                                 SingleR.label=="cDC1-Clec9a" | SingleR.label=="cDC1-Ccl22" | SingleR.label=="cDC2-Cd209a" |SingleR.label== "cDC2 Itgax"|
                                 SingleR.label=="T/NK cells" | SingleR.label=="Plasma cells" | SingleR.label=="Granulocytes" | SingleR.label=="B cells" 
))
scKPC[["SingleR.label"]] <- replace(scKPC[["SingleR.label"]], scKPC[["SingleR.label"]]=="myCAF-like PSC", "PSC")
tiff("treatment_tSNE_SingleR_labeled.jpeg", unit="in", width=15, height=8, res=300)
DimPlot(scKPC, reduction="umap",label=TRUE,pt.size=1,group.by = "SingleR.label", split.by = "treatment",repel = T)
dev.off()


treatment_split <- SplitObject(scKPC, split.by = "treatment")
control <- treatment_split$control
TAK981 <- treatment_split$TAK981

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
cellchat <- filterCommunication(cellchat, min.cells = 40)

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
par(mfrow = c(2,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

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
pathways.show <- c("MHC-I") #VCAM, TBFb, VEGF, PD-L1, MHCI, 


# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver =seq(1,7)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "hierarchy")

# Circle plot
pathways.show <- c("CXCL") #VCAM, TGFb, VEGF, PD-L1, MHCI, CD80, CD86, SELL
par(mfrow=c(2,2))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")


# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
# names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout="hierarchy")


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
control_cellchatDB <- control_cellchatDB.mouse 
control_cellchat@DB <- control_cellchatDB
control_cellchat <- subsetData(control_cellchat)
control_cellchat <- identifyOverExpressedGenes(control_cellchat)
control_cellchat <- identifyOverExpressedInteractions(control_cellchat)
control_cellchat <- projectData(control_cellchat, PPI.mouse)
control_cellchat <- computeCommunProb(control_cellchat)
control_cellchat <- filterCommunication(control_cellchat, min.cells = 40)
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
TAK981_cellchatDB <- TAK981_cellchatDB.mouse 
TAK981_cellchat@DB <- TAK981_cellchatDB
TAK981_cellchat <- subsetData(TAK981_cellchat)
TAK981_cellchat <- identifyOverExpressedGenes(TAK981_cellchat)
TAK981_cellchat <- identifyOverExpressedInteractions(TAK981_cellchat)
TAK981_cellchat <- projectData(TAK981_cellchat, PPI.mouse)
TAK981_cellchat <- computeCommunProb(TAK981_cellchat)
TAK981_cellchat <- filterCommunication(TAK981_cellchat, min.cells = 40)
TAK981_cellchat <- computeCommunProbPathway(TAK981_cellchat)
TAK981_cellchat <- aggregateNet(TAK981_cellchat)
groupSize <- as.numeric(table(TAK981_cellchat@idents))
par(mfrow = c(2,2), xpd=TRUE)
netVisual_circle(TAK981_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(TAK981_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

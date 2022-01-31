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
setwd("E:/TAK981_KPC/control_vs_TAK981")

mref <- ImmGenData() #mouse immune cells
sceM <- BaronPancreasData('mouse')
sceM <- sceM[,!is.na(sceM$label)]
library(scuttle)
sceM <- logNormCounts(sceM)
## import all necessary files
## use the following command to unzip tar.gz: tar -xvzf
load("immune.RData")

load("TAK981_only.RData")

C1.data <- Read10X(data.dir="E:/TAK981_KPC/cellranger/c1_55M/filtered_feature_bc_matrix")
C2.data <- Read10X(data.dir="E:/TAK981_KPC/cellranger/c2_55M/filtered_feature_bc_matrix")
T1.data <- Read10X(data.dir="E:/TAK981_KPC/cellranger/TAK981_1/filtered_feature_bc_matrix")
T2.data <- Read10X(data.dir="E:/TAK981_KPC/cellranger/TAK981_2/filtered_feature_bc_matrix")

C1 <- CreateSeuratObject(counts=C1.data, project= "C1")
C1@meta.data$treatment <- "control"
#C1 <- NormalizeData(C1)
C2 <- CreateSeuratObject(counts=C2.data, project= "C2")
C2@meta.data$treatment <- "control"
#C2 <- NormalizeData(C2)
T1 <- CreateSeuratObject(counts=T1.data, project= "T1")
T1@meta.data$treatment <- "TAK981"
#T1 <- NormalizeData(T1)
T2 <- CreateSeuratObject(counts=T2.data, project= "T2")
T2@meta.data$treatment <- "TAK981"
#T2 <- NormalizeData(T2)

###############################################################################
#Sample by sample QC (https://www.nature.com/articles/s41388-021-02054-3#Sec22)
C1[["percent.mito"]] <- PercentageFeatureSet(C1, pattern="^mt-|Hsp")
C2[["percent.mito"]] <- PercentageFeatureSet(C2, pattern="^mt-|Hsp")
T1[["percent.mito"]] <- PercentageFeatureSet(T1, pattern="^mt-|Hsp")
T2[["percent.mito"]] <- PercentageFeatureSet(T2, pattern="^mt-|Hsp")

C1[["percent.ribo"]] <- PercentageFeatureSet(object = C1, pattern = "Rps|Rpl|Mrpl|Mrps")
C2[["percent.ribo"]] <- PercentageFeatureSet(object = C2, pattern = "Rps|Rpl|Mrpl|Mrps")
T1[["percent.ribo"]] <- PercentageFeatureSet(object = T1, pattern = "Rps|Rpl|Mrpl|Mrps")
T2[["percent.ribo"]] <- PercentageFeatureSet(object = T2, pattern = "Rps|Rpl|Mrpl|Mrps")


#nCount_RNA: number of UMIs per cell
#nFeature_RNA: number of genes detected per cell
VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA","percent.mito", "percent.ribo"), ncol = 4)
#VlnPlot(C2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(T1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(T2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(P2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#200-10,000 genes detected within cells, 1000-40,000 UMIs counted, less than 25% mitochondrial reads, less than 40% ribosombal genes
C1 <- subset(C1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA<40000 & percent.mito < 20 & percent.ribo < 40)
C2 <- subset(C2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA<40000 & percent.mito < 20 & percent.ribo < 40)
T1 <- subset(T1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA<40000 & percent.mito < 20 & percent.ribo < 40)
T2 <- subset(T2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA<40000 & percent.mito < 20 & percent.ribo < 40)

###############################################################################

#merge all samples (traditional seurat)
# KPC <- merge(C1, y=c(C2, T1, T2, P1, P2), add.cell.ids = c("C1","C2","T1","T2", "P1", "P2"), project="KPC")
# KPC <- NormalizeData(KPC)
# table(KPC$orig.ident)

#feature selection (highly variable features)

# KPC <- FindVariableFeatures(KPC, nfeatures = 3000)
# KPC.top10 <- head(VariableFeatures(KPC),10)
# KPC.plot1 <- VariableFeaturePlot(KPC)
# KPC.plot2 <- LabelPoints(plot=KPC.plot1, points=KPC.top10, repel=TRUE)
# KPC.plot1 + KPC.plot2

#data scaling

# all.genes <- rownames(KPC)
# gc()
# memory.limit(56000)
# plan("multiprocess", workers = 4)
# KPC <- ScaleData(KPC, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
# plan("multiprocess", workers = 1)

#Run the following line if they don't allow large vector allocation
#KPC <- ScaleData(KPC)


###############################################################################
#USING SCTRANSFORM INSTEAD. Better than the traditional approach above.
KPC <- merge(C1, y=c(C2, T1, T2), add.cell.ids = c("C1","C2","T1","T2"), project="KPC")
KPC.list <- SplitObject(KPC, split.by = "orig.ident")

#Pairwise merging
#control <- merge(C1, y=C2, add.cell.ids = c("C1","C2"), project="KPC")
#TAK <- merge(T1, y=T2, add.cell.ids = c("T1","T2"), project="KPC")
#TAKP <- merge(P1, y=P2, add.cell.ids = c("P1", "P2"), project="KPC")
#KPC.list <- list(control, TAK, TAKP)


#KPC.list <- SplitObject(KPC, split.by = "treatment")
#KPC.list <- lapply(X = KPC.list, FUN = SCTransform, vars.to.regress=c("nCount_RNA","percent.mt", "percent.ribo"))

for (i in 1:length(KPC.list)) {
  KPC.list[[i]] <- SCTransform(KPC.list[[i]], vars.to.regress=c("nCount_RNA","percent.mito", "percent.ribo"))
}

features <- SelectIntegrationFeatures(object.list = KPC.list, nfeatures = 2500)
KPC.list <- PrepSCTIntegration(object.list = KPC.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = KPC.list, normalization.method = "SCT",
                                         anchor.features = features)
memory.limit(56000)
combined.KPC <- IntegrateData(anchorset = anchors, normalization.method = "SCT")


###############################################################################

#PCA

#control <- RunPCA(control, feautures=VariableFeatures(object=control))
#VizDimLoadings(control, dims=1:2, reduction="pca")
#DimPlot(control, reduction="pca")

#TAK <- RunPCA(TAK, feautures=VariableFeatures(object=TAK))
# #TP <- RunPCA(TP, feautures=VariableFeatures(object=TP))
# 
# KPC <- RunPCA(KPC, feautures=VariableFeatures(object=KPC))
# #VizDimLoadings(KPC, dims=1:2, reduction="pca")
# #DimPlot(KPC, reduction="pca")
# #DimHeatmap(KPC, dims = 1:15, cells = 500, balanced = TRUE)
# 
# ###############################################################################
# #Clustering
# KPC <- FindNeighbors(KPC, dims=1:15)
# KPC <- FindClusters(KPC, resolution = 1)
# 
# ###############################################################################
# #t-SNE
# KPC <- RunTSNE(KPC, dims=1:15,label=TRUE, check_duplicates = FALSE)
# jpeg("tSNE.jpeg", width=1305, height=1092)
# DimPlot(KPC, reduction="tsne")
# dev.off()
# 
# ###############################################################################
# #UMAP
# KPC <- RunUMAP(KPC, dims=1:15)
# jpeg("UMAP.jpeg", width=700, height=500)
# DimPlot(KPC, reduction="umap", label=TRUE)
# dev.off()
# 
# DimPlot(KPC, reduction="umap", label=TRUE, split.by = "treatment")
# ###############################################################################
# #differential expressed features (cluster biomarkers)
# 
# #find markers for specific cluster
# #control.cluster9 <-FindMarkers(control, ident.1 = 2, min.pct = 0.25)
# #head(control.cluster9, n=5)
combined.KPC <- ScaleData(object = combined.KPC, verbose = FALSE)
combined.KPC <- RunPCA(combined.KPC, verbose = FALSE)
combined.KPC <- RunTSNE(combined.KPC, dims = 1:40, verbose = FALSE)
combined.KPC <- RunUMAP(combined.KPC, dims = 1:40, verbose = FALSE)

combined.KPC <- FindNeighbors(combined.KPC, dims = 1:40, verbose = FALSE)
combined.KPC <- FindClusters(combined.KPC, resolution = 1.2)
DimPlot(combined.KPC, reduction="tsne",label=TRUE,pt.size=1)

# find markers for all clusters
# MUST USE "RNA" for assay when using samples with different conditions (i.e. control vs treated)
DefaultAssay(combined.KPC) <- "RNA"
#Must normalize data with Lognormalize for further DE analyses
combined.KPC <- NormalizeData(object =combined.KPC, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(combined.KPC)
combined.KPC <- ScaleData(object = combined.KPC, features = all.genes)
KPC.markers <- FindAllMarkers(combined.KPC, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
KPC.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

####################################################################################
#Heatmap for top 10 markers for every cluster
tiff("cluster_heatmap2.jpeg", unit="in", width=40, height=24, res=500)
DoHeatmap(combined.KPC, features = top10$gene) + NoLegend()
dev.off()

###############################################################################
tiff("tSNE_unlabeled.jpeg", unit="in", width=8, height=7, res=500)
DimPlot(combined.KPC, reduction="tsne",label=TRUE,pt.size=1)
dev.off()

DimPlot(combined.KPC, reduction="tsne",label=TRUE,pt.size=1)
VlnPlot(combined.KPC, "Ptprc", pt.size=0.5)+NoLegend()
immune <- subset(combined.KPC, idents=c(0:23)[-c(1,4,5,8,9,10,12,16,17,19,20,21)])
# cluster 5(or #6) is mt-genes. exclding...
nonimmune <- subset(combined.KPC, idents=c(0:23)[c(1,4,5,8,9,10,12,16,17,19,20,21)])

#Feature Plot (reference to https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00776-9#ethics)
FeaturePlot(combined.KPC, features=c("Ptprc", "Cd3d", "Cd19","Adgre1", "Itgam", "Ly6g", "Itgax"), reduction = "tsne")

####################################################################################################################
#automated cell type identification

#convert Seurat to SingleR-compatible dataset
scImmune <- as.SingleCellExperiment(immune)
scImmune <- SingleR(test=scImmune, ref=mref, assay.type.test = 1, labels = mref$label.main)
immune[["SingleR.label"]] <- scImmune$labels
rm(scImmune)
immune[["SingleR.label"]] <- replace(immune[["SingleR.label"]], immune[["SingleR.label"]]=="Microglia", "Macrophages")
immune[["SingleR.label"]] <- replace(immune[["SingleR.label"]], immune[["SingleR.label"]]=="Tgd", "T cells")
immune[["SingleR.label"]] <- replace(immune[["SingleR.label"]], immune[["SingleR.label"]]=="ILC", "NK cells")

#eliminate epithelial cells, stem cells, etc
# immune <- subset(immune, subset=(SingleR.label=="Monocytes" |SingleR.label == "Macrophages"|SingleR.label == "Basophils" |
#                                    SingleR.label == "DC" | SingleR.label == "T cells" | SingleR.label == "NK cells" | 
#                                    SingleR.label == "B cells" | SingleR.label == "ILC" | SingleR.label == "NKT" |
#                                    SingleR.label == "Tgd" | SingleR.label == "Mast cells" | SingleR.label == "Neutrophils"))
DimPlot(immune, reduction="tsne", label=TRUE,pt.size=1, repel = TRUE)
DimPlot(immune, reduction="tsne", label=TRUE,pt.size=1, group.by = "SingleR.label", repel = TRUE)


#combined.KPC <- subset(combined.KPC, subset= Ptprc > 0)
DimPlot(nonimmune, reduction="tsne", label=TRUE,pt.size=1)

scNonImmune <- as.SingleCellExperiment(nonimmune)
scNonImmune <- SingleR(test=scNonImmune, ref=sceM, assay.type.test = 1, labels = sceM$label)
nonimmune[["SingleR.label"]] <- scNonImmune$labels
rm(scNonImmune)
# nonimmune[["SingleR.label"]] <- replace(nonimmune[["SingleR.label"]], nonimmune[["SingleR.label"]]=="ductal", "Ductal 1")
# #nonimmune[["SingleR.label"]] <- replace(nonimmune[["SingleR.label"]], nonimmune[["SingleR.label"]]=="B_cell", "Ductal 1")
# nonimmune[["SingleR.label"]] <- replace(nonimmune[["SingleR.label"]], nonimmune[["SingleR.label"]]=="schwann", "Ductal 2")
# nonimmune[["SingleR.label"]] <- replace(nonimmune[["SingleR.label"]], nonimmune[["SingleR.label"]]=="endothelial", "Endothelial")
# nonimmune[["SingleR.label"]] <- replace(nonimmune[["SingleR.label"]], nonimmune[["SingleR.label"]]=="quiescent_stellate", "Fibroblasts")
# nonimmune[["SingleR.label"]] <- replace(nonimmune[["SingleR.label"]], nonimmune[["SingleR.label"]]=="activated_stellate", "Fibroblasts")

# tiff("nonimmune_temp_tsne.jpeg", unit="in", width=6, height=5, res=500)
# DimPlot(nonimmune, reduction="tsne", label=TRUE,pt.size=0.5, group.by = "SingleR.label", repel = TRUE)
# dev.off()

# nonimmune <- subset(nonimmune, subset=(SingleR.label=="Ductal 1" |SingleR.label == "Ductal 2"|
#                                    SingleR.label == "Fibroblasts 1"|SingleR.label == "Fibroblasts 2"|SingleR.label == "endothelial"))
DimPlot(nonimmune, reduction="tsne", label=TRUE,pt.size=1, group.by = "SingleR.label", repel = TRUE)


# nonimmune <- subset(combined.KPC, idents=c(0:24)[c(1,4,7,9,10,11,12,16,18,19,21,22,23)])
nonimmune.markers <- FindAllMarkers(nonimmune, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
nonimmune.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> nonimmune.top10

####################################################################################
#Heatmap for top 10 markers for every cluster
tiff("nonimmune_initial_cluster_heatmap.jpeg", unit="in", width=40, height=24, res=500)
DoHeatmap(nonimmune, features = nonimmune.top10$gene) + NoLegend()
dev.off()
DimPlot(nonimmune, reduction="tsne", label=TRUE,pt.size=1, repel = TRUE)
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==9, "Fibroblasts")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==18, "Fibroblasts") #SAA3 and C3 are highly expressed compared to Fibroblasts 1 
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==16, "Endothelial")

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==19, "Acinar")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==20, "Ductal NM") #the only EpCAM+
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==0, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==8, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==7, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==3, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==15, "Ductal 1")

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==4, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==11, "Ductal 2")


DimPlot(immune, reduction="umap", label=TRUE,pt.size=1, group.by = "SingleR.label", repel = TRUE)


save(immune, file="immune.RData")
save(nonimmune, file="nonimmune.RData")

# rm(list=setdiff(ls(), "kpc"))

#Merge labled immune and nonimmune cells
scKPC<-merge(immune, y=c(nonimmune))
scKPC <- DietSeurat(scKPC, assay="RNA")
scKPC.list <- SplitObject(scKPC, split.by = "orig.ident")

for (i in 1:length(scKPC.list)) {
  scKPC.list[[i]] <- SCTransform(scKPC.list[[i]], vars.to.regress=c("nCount_RNA","percent.mito", "percent.ribo"))
}
features <- SelectIntegrationFeatures(object.list = scKPC.list, nfeatures = 2000)
scKPC.list <- PrepSCTIntegration(object.list = scKPC.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = scKPC.list, normalization.method = "SCT",
                                  anchor.features = features)
memory.limit(56000)
scKPC <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
scKPC <- ScaleData(object = scKPC, verbose = FALSE)
scKPC <- RunPCA(scKPC, verbose = FALSE)
scKPC <- RunTSNE(scKPC, dims = 1:40, verbose = FALSE)
scKPC <- RunUMAP(scKPC, dims = 1:40, verbose = FALSE)
scKPC <- FindNeighbors(scKPC, dims = 1:40, verbose = FALSE)
scKPC <- FindClusters(scKPC, resolution = 1.4)

DimPlot(scKPC, reduction="tsne",label=TRUE,pt.size=1,group.by = "SingleR.label", repel = T)

# #Must normalize data with Lognormalize for further DE analyses
# DefaultAssay(scKPC) <- "RNA"
# scKPC <- NormalizeData(scKPC) #LogNormalize
# all.genes <- rownames(scKPC)
# scKPC <- ScaleData(object = scKPC, features = all.genes)
# scKPC.markers <- FindAllMarkers(scKPC, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# scKPC.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC) -> scKPC.top10
# 
# tiff("cluster_heatmap.jpeg", unit="in", width=40, height=24, res=500)
# DoHeatmap(scKPC, features = scKPC.top10$gene) + NoLegend()
# dev.off()
# tiff("cluster_heatmap_SingleR.jpeg", unit="in", width=40, height=24, res=500)
# DoHeatmap(scKPC, features = scKPC.top10$gene, group.by = "SingleR.label") + NoLegend()
# dev.off()
# DimPlot(scKPC, reduction="tsne",label=TRUE,pt.size=1, repel = T)


tiff("tSNE_SingleR_labeled.jpeg", unit="in", width=10, height=10, res=500)
DimPlot(scKPC, reduction="tsne", label=TRUE,pt.size=1, group.by = "SingleR.label", repel = TRUE)
dev.off()
tiff("treatment_tSNE_SingleR_labeled.jpeg", unit="in", width=11, height=5, res=500)
DimPlot(scKPC, reduction="tsne", label=TRUE,pt.size=0.8, group.by = "SingleR.label", split.by = "treatment", repel = TRUE)
dev.off()
tiff("individual_tSNE_SingleR_labeled.jpeg", unit="in", width=20, height=5, res=500)
DimPlot(scKPC, reduction="tsne", label=TRUE,pt.size=1, group.by = "SingleR.label", split.by = "orig.ident", repel = TRUE)
dev.off()


DimPlot(scKPC, reduction="tsne", label=TRUE,pt.size=0.8, group.by = "SingleR.label", repel = TRUE)

nonimmune <- subset(scKPC, subset=(SingleR.label=="Acinar"| SingleR.label=="Ductal 0" | SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2" | 
                                     SingleR.label=="Ductal NM" | SingleR.label=="Endothelial" | SingleR.label=="Fibroblasts"))
save(nonimmune, file="nonimmune.RData")

##
#using scRNA dataset as reference
#library(scRNAseq)
#scref <- StoeckiusHashingData(mode='mouse')
#library(scuttle)
#scref <- logNormCounts(scref)
#scpred.KPC <- SingleR(test=sceG, ref=sceM, labels=sceM$label, de.method="wilcox")
#table(scpred.KPC$labels)
#plotScoreHeatmap(scpred.KPC)z
##
rm(scKPC2, sceM, scKPC.list, anchors, nonimmune, immune.temp, immune.list, temp)
rm(immune, nonimmune.markers, nonimmune.top10, immune.markers, mref, immune.temp)
save.image("E:/TAK981_KPC/control_vs_TAK981/TAK981_only.RData")

####################################################################################################################
load("immune.RData")
mref <- ImmGenData() #mouse immune cells
sceM <- BaronPancreasData('mouse')
sceM <- sceM[,!is.na(sceM$label)]
library(scuttle)
sceM <- logNormCounts(sceM)

immune <- DietSeurat(immune, assay="RNA")
immune.list <- SplitObject(immune, split.by = "orig.ident")

for (i in 1:length(immune.list)) {
  immune.list[[i]] <- SCTransform(immune.list[[i]], vars.to.regress=c("nCount_RNA","percent.mito", "percent.ribo"))
}
features <- SelectIntegrationFeatures(object.list = immune.list, nfeatures = 3000)
immune.list <- PrepSCTIntegration(object.list = immune.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = immune.list, normalization.method = "SCT",
                                  anchor.features = features)
memory.limit(56000)
immune <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

immune <- ScaleData(object = immune, vars.to.regress=c("nCount_RNA","percent.mito", "percent.ribo"))
immune <- RunPCA(immune, verbose = FALSE)
immune <- RunTSNE(immune, dims = 1:40, verbose = FALSE)
immune <- RunUMAP(immune, dims = 1:40, verbose = FALSE)
immune <- FindNeighbors(immune, reduction = "umap",dims = 1:2, verbose = FALSE)
immune <- FindClusters(immune, resolution = 0.2)
DimPlot(immune, reduction="umap",label=TRUE,pt.size=1)

# immune_saved <- immune

DefaultAssay(immune) <- "RNA"
immune <- NormalizeData(immune) #LogNormalize
all.genes <- rownames(immune)
immune <- ScaleData(object = immune, features = all.genes)
immune.markers <- FindAllMarkers(immune, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immune.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10
write.csv(immune.markers, "immune_markers2.csv", row.names = T)
jpeg("immuune_cluster_heatmap.jpeg", width=3000, height=2000)
DoHeatmap(immune, features = top10$gene) + NoLegend()
dev.off()

tiff("immune_unlabeled.jpeg", unit="in", width=6, height=5, res=500)
DimPlot(immune, reduction="umap",label=TRUE,pt.size=1)
dev.off()

DimPlot(immune, reduction="umap",label=TRUE,pt.size=1, group.by = "SingleR.label", repel=T)

S100a9 <- VlnPlot(immune, features = c("S100a9"), group.by = "SingleR.label") +
  ylab("S100a9") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

Clec9a <- VlnPlot(immune, features = c("Clec9a"), group.by = "SingleR.label") +
  ylab("Clec9a") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

Isg15 <- VlnPlot(immune, features = c("Isg15"), group.by = "SingleR.label") +
  ylab("Isg15") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

Mki67 <- VlnPlot(immune, features = c("Mki67"), group.by = "SingleR.label") +
  ylab("Mki67") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

C1qa <- VlnPlot(immune, features = c("C1qa"), group.by = "SingleR.label") +
  ylab("C1qa") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

Spp1 <- VlnPlot(immune, features = c("Spp1"), group.by = "SingleR.label") +
  ylab("Spp1") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

Ly6c2 <- VlnPlot(immune, features = c("Ly6c2"), group.by = "SingleR.label") +
  ylab("Ly6c2") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

Ccl22 <- VlnPlot(immune, features = c("Ccl22"), group.by = "SingleR.label") +
  ylab("Ccl22") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

Cd209a <- VlnPlot(immune, features = c("Cd209a"), group.by = "SingleR.label") +
  ylab("Cd209a") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Itgax <- VlnPlot(immune, features = c("Itgax"), group.by = "SingleR.label") +
  ylab("Itgax") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

Cd79a <- VlnPlot(immune, features = c("Cd79a"), group.by = "SingleR.label") +
  ylab("Cd79a") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Jchain <- VlnPlot(immune, features = c("Jchain"), group.by = "SingleR.label") +
  ylab("Jchain") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Cd3e <- VlnPlot(immune, features = c("Cd3e"), group.by = "SingleR.label") +
          ylab("Cd3e") +
          scale_y_continuous(position="left")+
          theme(axis.title.x =element_blank(),
                axis.text.x = element_text(size = 25),
                axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
                legend.position = 'none',
                plot.title = element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, size=1))

tiff("immune_vlnplot.jpeg", unit="in", width=24, height=20, res=500)
Cd79a/Ccl22/Clec9a/Cd209a/Itgax/S100a9/C1qa/Isg15/Ly6c2/Mki67/Spp1/Jchain/Cd3e
dev.off()

VlnPlot(Macrophages, features = c("Cx3cr1"), group.by = "SingleR.label")
VlnPlot(Macrophages, features = c("Hspa1a"), group.by = "SingleR.label", split.by = "treatment")
VlnPlot(immune, features = c("Ly6c2"))

VlnPlot(Macrophages, features = c("Isg15"))
VlnPlot(immune, features = c("Ccl3"))
VlnPlot(immune, features = c("Arg2"))

#macrophage
Macrophages <- subset(immune, subset=(SingleR.label=="Macro-Spp1" |SingleR.label=="Macro-Ly6c+Chil3+" |
                                        SingleR.label=="Macro-Ly6c+Isg15+" |SingleR.label=="Macro-Proliferating" | SingleR.label=="Macro-C1q"))
save(Macrophages, file="Macrophages.RData")
rm(Macrophages)
Spp1 <- subset(immune, subset=(SingleR.label=="Macro-Spp1"))
C1q <- subset(immune, subset=(SingleR.label=="Macro-C1q"))


#cIAP12 paper and drug response Cell paper
VlnPlot(Macrophages, features = c("Adgre1"))
VlnPlot(Macrophages, features = c("H2-Aa"))
VlnPlot(Macrophages, features = c("H2-Eb1"))
VlnPlot(Macrophages, features = c("Cd74"))

VlnPlot(Macrophages, features = c("Il1b"))

VlnPlot(Macrophages, features = c("C1qc"))

VlnPlot(Macrophages, features = c("Maf"))

VlnPlot(Macrophages, features = c("Mafb"))

dot_features <- c("Nfkb1", "Nfkb2",
                  "Chuk","Ikbkb", "Ikbkg", "Rel")
dot_features <- c("Atf1", "Atf2","Atf3","Atf4","Atf5","Atf6", "Atf7", "Batf", "Batf2", "Batf3", "Jdp2",
                  "Jun", "Junb", "Jund",
                  "Fos", "Fosb", "Fosl1", "Fosl2",
                  "Maf", "Mafa", "Mafb", "Maff", "Mafg", "Mafk")
DotPlot(Macrophages, features = dot_features, group.by="treatment") + scale_size(range = c(1,15)) + RotatedAxis() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

dot_features <- c("Spp1","Arg1","Vegfa", "Mki67","Top2a","Ube2c","Ly6c2","Isg15", "Plac8" ,"Ccr2","Chil3","Vcan","C1qb","Mafb", "Maf") 
dittoDotPlot(Macrophages, vars = dot_features, group.by = "SingleR.label",size=10, max.color="red", min.color = "yellow") +
  theme(axis.title.y=element_blank(),
        axis.text.x = element_text(angle=30)) 
# dittoDotPlot(Macrophages, vars = dot_features, group.by = "SingleR.label", size=10, max.color="red", min.color = "yellow") + coord_flip() +
#   theme(axis.title.x=element_blank(),
# axis.text.x = element_text(angle=30))


dittoHeatmap(Macrophages, dot_features,
             annot.by = c("SingleR.label", "treatment"))
tiff("macrophage_dittodotplot.jpeg", unit="in", width=10, height= 5, res=500)
dittoDotPlot(Macrophages, vars = dot_features, group.by = "SingleR.label",size=12, max.color="red", min.color = "yellow") +
  theme(axis.title.y=element_blank(),
        axis.text.y = element_text(face="bold"))
dev.off()
tiff("macrophage_dittodotplot_flipped.jpeg", unit="in", width=6, height= 5.5, res=500)
dittoDotPlot(Macrophages, vars = dot_features, group.by = "SingleR.label", size=8, max.color="red", min.color = "yellow") + coord_flip() +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(size=10,angle=30,face="bold"))
dev.off()
tiff("immune_dittoHeatmap.jpeg", unit="in", width=7, height=5, res=500)
dittoHeatmap(Macrophages, dot_features,
             annot.by = c("SingleR.label", "treatment"))
dev.off()

dot_features <- c("Cd80", "Cd86", "Sod2", "Tnf", "Nos2", "Stat1","Irf3","Irf5","Il12b","Tlr4","Il18", "Ccl2","Ccl3", "Ccl4", "Cxcl9", "Cxcl10", "Cxcl16",
                  "Mrc1", "Chil3", "Retnla", "Arg1", "Stat3",
                  "Il10ra", "Tnfsf14", "Il6", "Sphk1", 
                  "Cd163", "Tgfb1","Mertk", "Tlr8",
                  "Vegfa")

tiff("M1 M2 Macrophages.jpeg", unit="in", width=17, height=3.5, res=500)
DotPlot(Macrophages, features = dot_features, group.by="treatment") + scale_size(range = c(1,15)) + RotatedAxis() +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
dev.off()

# DotPlot(Macrophages, features = dot_features, group.by="treatment") + scale_size(range = c(1,10)) + 
#   geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
#   scale_colour_viridis(option="magma") +
#   guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
#   RotatedAxis() +
#   theme(axis.title = element_blank())

#Macro-Vegfa
VlnPlot(Macrophages, features = c("Vegfa"))
VlnPlot(Macrophages, features = c("Arg1"))
VlnPlot(Macrophages, features = c("Spp1"))
VlnPlot(Macrophages, features = c("Cxcl3")) # pro-angiogenic chemokines https://www.frontiersin.org/articles/10.3389/fphys.2013.00159/full

#Macro-Phagocytic
VlnPlot(Macrophages, features = c("Ly6c2"))

#Macro-ISG
VlnPlot(Macrophages, features = c("Isg15"))
VlnPlot(Macrophages, features = c("Stat1"))
VlnPlot(Macrophages, features = c("Cd86"))

# mref <- ImmGenData() #mouse immune cells

#cluster 13 is not immune cell. some fibroblasts and fibrocytes express cd45. 
temp <- immune

SRimmune <- as.SingleCellExperiment(temp)
SRimmune <- SingleR(test=SRimmune, ref=mref, assay.type.test = 1, labels = mref$label.main)
temp[["SingleR.label"]] <- SRimmune$labels
rm(SRimmune)
DimPlot(immune, reduction="umap", label=TRUE,pt.size=1, group.by = "SingleR.label", repel = TRUE)

VlnPlot(Macrophages, features = c("Mafb"))

# temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==13, "Epithelial")
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==2, "B cells")
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==27, "Plasma cells") #
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==10, "T/NK cells")
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==24, "T/NK cells")
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==28, "T/NK cells")

temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==15, "Granulocytes") #Il1b+ Arg2+  https://www.science.org/doi/10.1126/sciimmunol.aay6017?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==22, "Granulocytes")
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==30, "Granulocytes")
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==31, "Granulocytes")

# https://www.cell.com/pb-assets/products/nucleus/nucleus-phagocytes/rnd-systems-dendritic-cells-br.pdf
# https://www.science.org/doi/10.1126/scitranslmed.abf5058
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==20, "cDC2-Cd209a") #Cd11b+ Cd11c(int) Cdh1+ Clec9a- Xcr1- CD209a+
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==25, "cDC1-Ccl22") #Relb+ Itgae(Cd103)- CD8a-  Xcr1- Cd209a- Clec9a- Ly75(Cd205)+ 
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==23, "cDC2-Itgax") #expressing Itgax and Mgl2, Ear2, 
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==26, "cDC1-Clec9a") #Itgae(Cd103)+ CD8a-  Xcr1+  Clec9a+ Batf3+
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==29, "pDC") #Cd11b- Siglech+ Bst2+ Tlr7+ Clec9a+ Ly6c2+ 

# M1 high tumor also proliferative
# https://www.nature.com/articles/s41598-020-73624-w
# https://www.nature.com/articles/s41467-017-01711-0
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==3, "Macro-C1q") #MHC Class II high macropahge (cd74 is invariant MHCii) https://www.frontiersin.org/articles/10.3389/fimmu.2018.01132/full
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==7, "Macro-C1q")
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==8, "Macro-C1q")

temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==4, "Macro-Ly6c+Chil3+") #high Il1b, Ly6c2, Ccr2. pro-tumorogensis genes such as Lyz2, Ifitm3, Vim, S100a6, https://www.jimmunol.org/content/206/1_Supplement/101.01
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==17, "Macro-Ly6c+Chil3+") #high Il1b, Ly6c2, Ccr2. pro-tumorogensis genes such as Lyz2, Ifitm3, Vim, S100a6, https://www.jimmunol.org/content/206/1_Supplement/101.01

temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==5, "Macro-Ly6c+Isg15+") #Ly6C+CCR2+ indicates BMDM. highly expressed genes are related to IFN response (IFIT, ISG, OAS). 


temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==11, "Macro-Proliferating") #Ki67 high
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==12, "Macro-Proliferating") #Ki67 high

temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==9, "Macro-Spp1") #high in Arg1. may be MDSC. determine after InferCNV
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==14, "Macro-Spp1") 
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==16, "Macro-Spp1") 
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==21, "Macro-Spp1") 


temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==0, "Macro-C1q") #Heat shock protein (HSP) high macrophages
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==1, "Macro-C1q")
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==6, "Macro-C1q")
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==13, "Macro-C1q")
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==18, "Macro-C1q") #HSP70 cluster
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==19, "Macro-C1q")

immune <- temp
# rm(temp)

tiff("immune_umap_SingleR.jpeg", unit="in", width=7, height=6, res=300)
DimPlot(immune, reduction="umap", label=TRUE,pt.size=0.8, group.by = "SingleR.label", repel = TRUE) + theme(legend.position = "none")
dev.off()
tiff("immune_umap_SingleR_treatment.jpeg", unit="in", width=18, height=6, res=300)
DimPlot(immune, reduction="umap", label=TRUE,pt.size=1, group.by = "SingleR.label", split.by = "treatment", repel = TRUE) + theme(legend.position = "none")
dev.off()

tiff("immune_celltype_proportion_treatment.jpeg", unit="in", width=9, height=6, res=500)
plot_group_proportions(immune, graph.type = "dodge")
dev.off()

tiff("immune_celltype_proportion_treatment_stacked.jpeg", unit="in", width=3.5, height=6, res=500)
plot_group_proportions(immune, graph.type = "stacked")
dev.off()

########################################################################################
#Proportion comparison - montecarlo
library("scProportionTest")
immune_prop <- sc_utils(immune)

prop_test <- permutation_test(
  immune_prop, cluster_identity = "SingleR.label",
  sample_1 = "TAK981", sample_2 = "control",
  sample_identity = "treatment"
)
permutation_plot(prop_test)+ scale_colour_discrete(labels = function(x) str_wrap(x, width = 7)) 

tiff("immune_celltype_proportion_test.jpeg", unit="in", width=7, height=5, res=500)
permutation_plot(prop_test)+ scale_colour_discrete(labels = function(x) str_wrap(x, width = 7))
dev.off()

########################################################################################
#Slingshot map psuedotime
# library(scater)
# scimmune <- as.SingleCellExperiment(immune, assay="RNA")
# scimmune <- runPCA(scimmune)
# scimmune <- slingshot(scimmune)  # no clusters
scimmune <- immune
Idents(scimmune) <- "SingleR.label"
scimmune <- slingshot(Embeddings(scimmune, "umap"), clusterLabels = scimmune@active.ident, 
                      start.clus = "Monocytes/MDSCs", stretch = 0)
lnes <- getLineages(reducedDim(scimmune,"PCA"),
                    scimmune@clusterLabels, start.clus = "Monocytes/MDSCs")

library(Polychrome)
library(ggbeeswarm)
library(ggthemes)

# this define the cluster color. You can change it with different color scheme.
my_color <- createPalette(length(levels(colnames(scimmune@clusterLabels))), c("red", "blue"), M=1000)
names(my_color) <- unique(as.character(colnames(scimmune@clusterLabels)))

colData(scimmune)
slingshot_df <- data.frame(colData(scimmune))

plot(reducedDims(scimmune), col = my_color[as.character(colnames(scimmune@clusterLabels))], 
     pch=16, 
     asp = 1)
legend("bottomleft",legend = names(my_color[levels(colnames(scimmune@clusterLabels))]),  
       fill = my_color[levels(colnames(scimmune@clusterLabels))])
lines(SlingshotDataSet(lnes), lwd=2, type = 'lineages', col = c("black"))

plot(reducedDim(scimmune), col=c("black"),pch = 16, cex = 0.5)
lines(scimmune, lwd = 2, type = 'lineages', col = 'black')



########################################################################################
# Trajectory analysis
# https://github.com/cole-trapnell-lab/monocle-release/issues/388
# https://github.com/satijalab/seurat/issues/1658
myeloid <- subset(immune, subset=(SingleR.label=="Macro-Spp1" |SingleR.label=="Macro-Ly6c+Chil3+" |
                                    SingleR.label=="Macro-Ly6c+Isg15+" |SingleR.label=="Macro-Proliferating" | SingleR.label=="Macro-C1q" |
                                    SingleR.label=="cDC1-Clec9a" | SingleR.label== "cDC1-Ccl22" | SingleR.label=="cDC2-Cd209a" | 
                                    SingleR.label== "cDC2-Itgax"))
DimPlot(myeloid, reduction="umap", label=TRUE,pt.size=1, group.by = "SingleR.label", repel = TRUE)
myeloid <- ProjectDim(myeloid, reduction = "pca")
myeloid_matrix <- myeloid@assays$RNA@counts
myeloid_metadata <- myeloid@meta.data
myeloid_annot <- data.frame(gene_short_name = rownames(myeloid@assays$RNA), row.names = rownames(myeloid@assays$RNA))

myeloid.cds <- new_cell_data_set(
                myeloid_matrix,
                cell_metadata = myeloid_metadata,
                gene_metadata = myeloid_annot)

# rowData(myeloid.cds)$gene_short_name <- row.names(rowData(myeloid.cds))
# option 1: preprocess cds with monocle
myeloid.cds <- preprocess_cds(myeloid.cds, num_dim = 100)
plot_pc_variance_explained(myeloid.cds)
myeloid.cds = align_cds(myeloid.cds, num_dim = 100, alignment_group = "orig.ident")
myeloid.cds <- reduce_dimension(myeloid.cds, preprocess_method="PCA")
myeloid.cds <- cluster_cells(myeloid.cds, reduction_method = "UMAP", resolution = 0.1)
plot_cells(myeloid.cds, color_cells_by = "SingleR.label", group_label_size = 3.5, label_cell_groups=T, label_groups_by_cluster=F)


# Option 2: Transfer seurat embedding
reducedDim(myeloid.cds, type = "PCA") <- myeloid@reductions$pca@cell.embeddings 
myeloid.cds@preprocess_aux$prop_var_expl <- myeloid@reductions$pca@stdev
plot_pc_variance_explained(myeloid.cds)
#transfer umap
myeloid.cds@int_colData@listData$reducedDims$UMAP <- myeloid@reductions$umap@cell.embeddings
#copy cluster info from Seruat
myeloid.cds@clusters$UMAP_so$clusters <- myeloid@meta.data$integrated_snn_res.1
myeloid.cds <- cluster_cells(myeloid.cds, reduction_method = "UMAP", resolution = 0.01)
rownames(myeloid.cds@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(myeloid.cds@int_colData@listData$reducedDims$UMAP) <- NULL

#    plot_cells(my.cds)
# cds <- align_cds(cds, alignment_group = "batch")
# myeloid.cds = align_cds(myeloid.cds, num_dim = 100, alignment_group = "orig.ident")
# myeloid.cds <- reduce_dimension(myeloid.cds, preprocess_method="PCA",reduction_method="UMAP")
# myeloid.cds <- cluster_cells(myeloid.cds, reduction_method="UMAP",resolution=1)
myeloid.cds <- learn_graph(myeloid.cds)
plot_cells(myeloid.cds,color_cells_by = "SingleR.label", label_cell_groups=T, label_groups_by_cluster=F,
           group_label_size=4, graph_label_size = 4,cell_size = 1.2)

tiff("Macrophage_trajectory.jpeg", unit="in", width=6, height=6, res=300)
plot_cells(myeloid.cds,color_cells_by = "SingleR.label", label_cell_groups=T, label_groups_by_cluster=F,
           group_label_size=4, graph_label_size = 4,cell_size = 1)
dev.off()



myeloid.cds3d <- reduce_dimension(myeloid.cds,max_components = 3)
myeloid.cds3d <- cluster_cells(myeloid.cds3d,reduction_method = "UMAP", resolution = 0.1)
myeloid.cds3d <- learn_graph(myeloid.cds3d)
plot_cells_3d(myeloid.cds3d,color_palette = c("red","blue","orange","grey","yellow"),color_cells_by = "SingleR.label",cell_size=30)

marker_test_res <- top_markers(myeloid.cds, group_cells_by="SingleR.label", 
                               reference_cells=1000, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
plot_genes_by_group(myeloid.cds,
                    top_specific_marker_ids,
                    group_cells_by="SingleR.label",
                    ordering_type="maximal_on_diag",
                    max.size=5)

#Single cell Trajectory 
myeloid.cds <- order_cells(myeloid.cds)
plot_cells(myeloid.cds,color_cells_by = "pseudotime", label_cell_groups=T, label_groups_by_cluster=F,
           label_leaves=TRUE, label_branch_points=TRUE,
           group_label_size=4, graph_label_size = 4,cell_size = 1, scale_to_range = T)
tiff("Macrophage_pseudotime.jpeg", unit="in", width=7, height=6, res=300)
plot_cells(myeloid.cds,color_cells_by = "pseudotime", label_cell_groups=T, label_groups_by_cluster=F,
           label_leaves=TRUE, label_branch_points=TRUE,
           group_label_size=4, graph_label_size = 4,cell_size = 1, scale_to_range = T)
dev.off()


#######################################################################################################################
#Pathway enrichment
immune.cerebro <- immune

Idents(immune.cerebro) <- "seurat_clusters"
immune.cerebro <- BuildClusterTree(
  immune.cerebro,
  dims = 1:30,
  reorder = FALSE,
  reorder.numeric = FALSE
)
immune.cerebro@misc$trees$seurat_clusters <- immune.cerebro@tools$BuildClusterTree

Idents(immune.cerebro) <- "SingleR.label"
immune.cerebro <- BuildClusterTree(
  immune.cerebro,
  dims = 1:30,
  reorder = FALSE,
  reorder.numeric = FALSE
)
immune.cerebro@misc$trees$SingleR.label <- immune.cerebro@tools$BuildClusterTree

library(cerebroApp)
immune.cerebro <- addPercentMtRibo(
  immune.cerebro,
  organism = 'mm',
  gene_nomenclature = 'name'
)
immune.cerebro <- getMostExpressedGenes(
  immune.cerebro,
  assay = 'RNA',
  groups = c('seurat_clusters','SingleR.label')
)


immune.cerebro <- getMarkerGenes(
  immune.cerebro,
  assay = 'RNA',
  organism = 'mm',
  groups = c('seurat_clusters','SingleR.label'),
  name = 'immune.cerebro',
  only_pos = TRUE,
  min_pct = 0.25,
  thresh_logFC = 0.25,
  thres_p_val = 0.01
)

immune.cerebro <- getEnrichedPathways(
  immune.cerebro,
  databases = c("GO_Biological_Process_2018", "GO_Cellular_Component_2018",
                "GO_Molecular_Function_2018", "KEGG_2016", "WikiPathways_2016", "Reactome_2016",
                "Mouse_Gene_Atlas"),
  marker_genes_input = 'immune.cerebro',
  adj_p_cutoff = 0.01,
  max_terms = 100
)


#GSVA
# hallmark_pathway <- gmtPathways("h.all.v7.4.symbols.gmt.txt")
# launchCerebroV1.3()

immune.cerebro <- performGeneSetEnrichmentAnalysis(
  immune.cerebro,
  GMT_file = "h.all.v7.5.1.symbols.gmt",
  groups = c('seurat_clusters','SingleR.label'),
  thres_p_val = 0.05,
  thres_q_val = 0.1
)


library(monocle)
monocle <- newCellDataSet(
  immune.cerebro@assays$RNA@counts,
  phenoData = new('AnnotatedDataFrame', data = immune.cerebro@meta.data),
  featureData = new('AnnotatedDataFrame', data = data.frame(
    gene_short_name = rownames(immune.cerebro@assays$RNA@counts),
    row.names = rownames(immune.cerebro@assays$RNA@counts))
  )
)

monocle <- estimateSizeFactors(monocle)
monocle <- setOrderingFilter(monocle, immune.cerebro@assays$RNA@var.features)
monocle <- reduceDimension(monocle, max_components = 2, method = 'DDRTree')
monocle <- orderCells(monocle)

immune.cerebro <- extractMonocleTrajectory(monocle, immune.cerebro, 'highly_variable_genes')

exportFromSeurat(
  immune.cerebro,
  assay = 'RNA',
  slot = 'data',
  file = paste0('cerebro_pbmc_seurat_', Sys.Date(), '.crb'),
  experiment_name = 'pbmc',
  organism = 'hg',
  groups = c('seurat_clusters','SingleR.label'),
  nUMI = 'nCount_RNA',
  nGene = 'nFeature_RNA',
  add_all_meta_data = TRUE,
  verbose = FALSE
)

launchCerebro()

#######################################################################################################################
source("scFunctions.R")
plot_group_proportions(immune, graph.type = "dodge")
plot_group_proportions(immune, graph.type = "stacked")

tiff("immune_celltype_proportion_treatment.jpeg", unit="in", width=9, height=6, res=500)
plot_group_proportions(immune, graph.type = "dodge")
dev.off()

tiff("immune_celltype_proportion_treatment_stacked.jpeg", unit="in", width=3.5, height=6, res=500)
plot_group_proportions(immune, graph.type = "stacked")
dev.off()


plot_heatmap_proportions(immune, graph.type = "by.cell")

tiff("immune_proportion_heatmap_treatment.jpeg", unit="in", width=9, height=3, res=300)
plot_heatmap_proportions(immune, graph.type = "by.cell")
dev.off()

save(immune, file="immune.RData")
Macrophages <- subset(immune, subset=(SingleR.label=="Macrophages" |SingleR.label=="Macrophages-cycling"))
save(Macrophages, file="Macrophages.RData")
TNK <- subset(immune, subset=(SingleR.label=="T/NK cells"))
save(TNK, file="TNK.RData")
DC <- subset(immune, subset=(SingleR.label=="cDC1-Ccl22" |SingleR.label=="cDC1-Clec9a" |SingleR.label=="cDC2-Cd209a" |
                             SingleR.label=="pDC" ))
save(DC, file="DC.RData")
Granulocytes <- subset(immune, subset=(SingleR.label=="Granulocytes"))
save(Granulocytes, file="Granulocytes.RData")
B <- subset(immune, subset=(SingleR.label=="B cells" | SingleR.label=="Plasma cells"))
save(B, file="B.RData")

##############################################################################################################################
#Macrophages
load("Macrophages.RData")
# 
# SRMacrophages <- as.SingleCellExperiment(Macrophages)
# SRMacrophages <- SingleR(test=SRMacrophages, ref=mref, assay.type.test = 1, labels = mref$label.main)
# Macrophages[["SingleR.label"]] <- SRMacrophages$labels
# rm(SRMacrophages)
# DimPlot(Macrophages, reduction="tsne", label=TRUE,pt.size=3, group.by = "SingleR.label", split.by = "treatment")
# Macrophages[["SingleR.label"]] <- replace(Macrophages[["SingleR.label"]], Macrophages[["SingleR.label"]]=="Microglia", "Macrophages")

#new.cluster.ids <- c("Myeloids", "Myeloids")
#names(new.cluster.ids) <- levels(Myeloids)
#names(new.cluster.ids)
#Myeloids <- RenameIdents(Myeloids,new.cluster.ids)

#https://www.researchgate.net/figure/Expression-of-known-M1-M2-and-TAM-marker-genes-in-our-GAMs-data-set_fig4_272092097
#Fcgr3(Cd16), 
#Mariottoni, Paula; Jiang Simon, et al., Frontiers in Medicine 2021 
#M1 markers: Gbp5, H2-Ab1(MHCii), Cd80, Cd86, Nos2(iNOS), Stat1, TNFa
#M2 markers: Mrc1(Cd206), Msr1(Cd204), Cd163, Chil3 (Ym1), Tgm2



tiff("M1 M2 Macrophages.jpeg", unit="in", width=17, height=3.5, res=300)
DotPlot(Macrophages, features = dot_features, group.by="treatment") + scale_size(range = c(1,15)) + RotatedAxis() + 
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F") + theme(axis.title = element_blank())
dev.off()

tiff("M1 M2 Macrophages individual.jpeg", unit="in", width=17, height=5, res=300)
DotPlot(Macrophages, features = dot_features, group.by="orig.ident") + scale_size(range = c(1,15)) + RotatedAxis() + 
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F") + theme(axis.title = element_blank())
dev.off()


dittoDotPlot(Macrophages, vars = dot_features, group.by = "treatment", size=10)

dittoHeatmap(Macrophages, dot_features,
             annot.by = c("SingleR.label", "treatment"))


#Various M2 population
#Chil3/4 = Ym1/2
#Tnfsf14 = LIGHT
dot_features <- c("Mrc1", "Chil3", "Chil4", "Retnla", "Arg1", "Il1r1",  "Il1r2",
                  "Il10ra", "Tnf", "Tnfsf14", "Il6", "Sphk1",
                  "Cd163", "Tgfb1", "Mertk", "Tlr8", "Tlr1",
                  "Vegfa")

tiff("M2 subset dotplot.jpeg", unit="in", width=12, height=5, res=300)
DotPlot(Macrophages, features = dot_features, group.by="treatment") + scale_size(range = c(1,15)) + RotatedAxis() + 
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()


#M1 polarization pathway
#NFkBcytokines:  "Il1b",, "Il6", "TNF" 
dot_features <- c("Stat5a","Irf5", "Stat1", "Socs3", "Irf3", "Ifnb1","Tlr4","Nfkb1", "Nfkb2",
                  "Il1b", "Il6", "Il15", "Tnf",
                  "Il10ra", "Stat3")

tiff("M1 polarization genes.jpeg", unit="in", width=7, height=5, res=300)
DotPlot(Macrophages, features = dot_features, group.by="treatment") + scale_size(range = c(3,10)) + RotatedAxis() + 
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()

##############################################################################################################
#simple GSEA
library(escape)
GS.hallmark <- getGeneSets(library = "H") #Hallmark gene sets
Macrophages.ES <- enrichIt(obj = Macrophages, gene.sets = GS.hallmark, groups = 1000, cores = 2)
Macrophages <- Seurat::AddMetaData(Macrophages, Macrophages.ES)

colors <- colorRampPalette(c("#0348A6", "#7AC5FF", "#C6FDEC", "#FFB433", "#FF4B20"))
dittoHeatmap(Macrophages, genes = NULL, metas = names(Macrophages.ES), 
             order.by = "treatment",
             annot.by = "treatment", 
             fontsize = 7, 
             cluster_cols = TRUE,
             heatmap.colors = colors(50))

#############################################################################################################
# 
# #EnrichR
# library(enrichR)
# dbs <- listEnrichrDbs()
# GO <- gmtPathways("GO_Biological_Process_2021.txt")
# Mac <- Macrophages
# Mac[["orig.ident"]] <- replace(Mac[["orig.ident"]], Mac[["orig.ident"]]=="C1", "control")
# Mac[["orig.ident"]] <- replace(Mac[["orig.ident"]], Mac[["orig.ident"]]=="C2", "control")
# Mac[["orig.ident"]] <- replace(Mac[["orig.ident"]], Mac[["orig.ident"]]=="T1", "TAK981")
# Mac[["orig.ident"]] <- replace(Mac[["orig.ident"]], Mac[["orig.ident"]]=="T2", "TAK981")
# Idents(Mac) <- "orig.ident" #set all idents to orig.ident i.e. control, TAK981, etc
# DEenrichRPlot(Mac, ident.1 = "TAK981", ident.2 =  "control", enrich.database = "WikiPathways_2019_Mouse", max.genes=10000)
# 

#############################################################################################################
# gseGO
library("org.Mm.eg.db", character.only = TRUE)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(ggnewscale)
organism <-org.Mm.eg.db

Macrophages_ranked_list <- FindMarkers(Mac, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)
# order list, pull out gene name and log2fc, and convert genes to uppercase
Macrophages_ranked_list <- Macrophages_ranked_list[order(Macrophages_ranked_list$avg_log2FC, decreasing = T),]
Macrophages_ranked_list$Gene.name <- rownames(Macrophages_ranked_list)
Macrophages_ranked_list <- Macrophages_ranked_list[,c("Gene.name", "avg_log2FC")]
rownames(Macrophages_ranked_list) <- NULL
head(Macrophages_ranked_list)

Macrophages_ranked_list2 <- prepare_ranked_list(Macrophages_ranked_list)
CAP_Macrophages_ranked_list <-Macrophages_ranked_list
CAP_Macrophages_ranked_list$Gene.name <- str_to_upper((Macrophages_ranked_list$Gene.name))

gse <- gseGO(geneList=Macrophages_ranked_list2, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 100, 
             maxGSSize = 500, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
goplot(gse)
require(DOSE)
enrichplot::dotplot(gse, showCategory=15, split=".sign", title="GO Biological Process Enrichment - Macrophage TAK981 treatment", font.size=10) + facet_grid(.~.sign) +
  scale_y_discrete(labels=function(gse) str_wrap(gse, width=20))
enrichplot::dotplot(gse, showCategory=20)


tiff("Macrophaegs_gseGO.jpeg", unit="in", width=10, height=12, res=300)
enrichplot::dotplot(gse, showCategory=15, split=".sign", title="GO Biological Process Enrichment - Macrophage TAK981 treatment", font.size=8) + facet_grid(.~.sign) +
  scale_y_discrete(labels=function(gse) str_wrap(gse, width=40))
dev.off()

#############################################################################################################
Macrophages.temp <- Macrophages
Macrophages.temp[["orig.ident"]] <- replace(Macrophages.temp[["orig.ident"]], Macrophages.temp[["orig.ident"]]=="C1", "control")
Macrophages.temp[["orig.ident"]] <- replace(Macrophages.temp[["orig.ident"]], Macrophages.temp[["orig.ident"]]=="C2", "control")
Macrophages.temp[["orig.ident"]] <- replace(Macrophages.temp[["orig.ident"]], Macrophages.temp[["orig.ident"]]=="T1", "TAK981")
Macrophages.temp[["orig.ident"]] <- replace(Macrophages.temp[["orig.ident"]], Macrophages.temp[["orig.ident"]]=="T2", "TAK981")

Idents(Macrophages.temp) <- "orig.ident" #set all idents to orig.ident i.e. control, TAK981, etc

Macrophages_ranked_list <- FindMarkers(Macrophages.temp, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)

# order list, pull out gene name and log2fc, and convert genes to uppercase
Macrophages_ranked_list <- Macrophages_ranked_list[order(Macrophages_ranked_list$avg_log2FC, decreasing = T),]
Macrophages_ranked_list$Gene.name <- str_to_upper(rownames(Macrophages_ranked_list))
Macrophages_ranked_list <- Macrophages_ranked_list[,c("Gene.name", "avg_log2FC")]
rownames(Macrophages_ranked_list) <- NULL
Macrophages_ranked_list <- prepare_ranked_list(Macrophages_ranked_list)
head(Macrophages_ranked_list)

Macrophages_fgsea_results <- fgseaMultilevel(pathways = hallmark_pathway,
                                 stats = Macrophages_ranked_list,
                                 minSize = 15,
                                 maxSize = Inf)

Macrophages_fgsea_df<- Macrophages_fgsea_results %>% arrange (desc(NES)) %>% select (pathway, pval, padj, NES, size)

Macrophages_fgsea_results2 <- Macrophages_fgsea_results[1:(length(Macrophages_fgsea_results)-1)]
Macrophages_fgsea_results2 <- apply(Macrophages_fgsea_results2,7,as.character)

write.csv(Macrophages_fgsea_results2, "Macrophages_fgsea.csv",row.names = FALSE)

topPathwaysUp <- Macrophages_fgsea_results[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- Macrophages_fgsea_results[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

macrophage_gseatable <- plotGseaTable(hallmark_pathway[topPathways],
                          stats=Macrophages_ranked_list,
                          fgseaRes=Macrophages_fgsea_results, 
                          gseaParam = 0.5)

tiff("macrophage_gseatable.jpeg", unit="in", width=20, height=5, res=300)
plotGseaTable(hallmark_pathway[topPathways],
              stats=Macrophages_ranked_list,
              fgseaRes=Macrophages_fgsea_results, 
              gseaParam = 0.5)
dev.off()


waterfall_plot(fgsea_results, "Pathways enriched in TAK981 treated vs control Macrophages")

tiff("waterfall_control_vs_TAK981_macrophages.jpeg", unit="in", width=9, height=6, res=300)
waterfall_plot(fgsea_results, "Pathways enriched in TAK981-treated vs control macrophages")
dev.off()

# example of pathway highly enriched in treated
tiff("Macrophages_gsea_Ifna.jpeg", unit="in", width=5, height=3, res=300)
plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_ALPHA_RESPONSE" , Macrophages_ranked_list)
dev.off()
tiff("Macrophages_gsea_Ifng.jpeg", unit="in", width=5, height=3, res=300)
plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_GAMMA_RESPONSE" , Macrophages_ranked_list)
dev.off()
tiff("Macrophages_gsea_TNFa.jpeg", unit="in", width=5, height=3, res=300)
plot_enrichment(hallmark_pathway, "HALLMARK_TNFA_SIGNALING_VIA_NFKB" , Macrophages_ranked_list)
dev.off()


############################################################################################################
############################################################################################################
load("DC.RData")
DC <- subset(Myeloid, subset=SingleR.label=="DC")

dot_features <- c("Itgam","Itgax","H2-Ab1","Cd40", "Tlr7", "Tlr9")

DotPlot(DC, features = dot_features, group.by="orig.ident") + scale_size(range = c(1,15)) + RotatedAxis() + 
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")


############################################################################################################
#TNK analysis
load("TNK.RData")

mref <- ImmGenData()
setwd("E:/TAK981_KPC/control_vs_TAK981")
hallmark_pathway <- gmtPathways("h.all.v7.4.symbols.gmt.txt")

TNK <- DietSeurat(TNK, assay="RNA")
DefaultAssay(TNK) <- "RNA"
TNK <- NormalizeData(TNK) #LogNormalize
TNK <- FindVariableFeatures(TNK, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(TNK)
TNK <- ScaleData(object = TNK, features = all.genes)
TNK <- RunPCA (TNK, features = VariableFeatures(object = TNK))
TNK <- FindNeighbors(TNK, dims = 1:40)
TNK <- FindClusters(TNK, resolution = c(0.6))
TNK <- RunTSNE(TNK, dims = 1:40)

TNK.markers <- FindAllMarkers(TNK, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TNK.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10
jpeg("TNK_cluster_heatmap.jpeg", width=2000, height=900)
DoHeatmap(TNK, features = top10$gene) + NoLegend()
dev.off()
DimPlot(TNK, reduction="tsne",label=TRUE,pt.size=1, group.by = "SingleR.label",split.by = "treatment")

DimPlot(TNK, reduction="umap",label=TRUE,pt.size=3, split.by = "treatment")

SRTNK <- as.SingleCellExperiment(TNK)
SRTNK <- SingleR(test=SRTNK, ref=mref, assay.type.test = 1, labels = mref$label.main)
TNK[["SingleR.label"]] <- SRTNK$labels
rm(SRTNK)

TNK[["SingleR.label"]] <- replace(TNK[["SingleR.label"]], TNK[["SingleR.label"]]=="Tgd", "T cells")
TNK[["SingleR.label"]] <- replace(TNK[["SingleR.label"]], TNK[["SingleR.label"]]=="ILC", "NK cells")
#TNK <- subset(TNK, subset=(SingleR.label=="T cells" |SingleR.label == "NK cells"))
DimPlot(TNK, reduction="tsne", label=TRUE,pt.size=3, group.by = "SingleR.label", split.by = "treatment")

TNK <- addPercentMtRibo(
  TNK,
  organism = 'mm',
  gene_nomenclature = 'name'
)


############################################################################################################

tlymph <- subset(TNK, subset=(SingleR.label=="T cells"))
DimPlot(tlymph, reduction="tsne", label=TRUE,pt.size=3, group.by = "SingleR.label", split.by = "treatment")

####################################################################################################################
#CD8 T cell analysis
cd8t <- subset(tlymph, subset= Cd8a > 0 | Cd8b1 > 0)
DimPlot(cd8t, reduction="tsne", label=TRUE,pt.size=3, group.by = "SingleR.label", split.by = "treatment")

#Il2ra(Cd25), Tnfrsf4(Ox40), Tnfrsf9(41bb)
dot_features <- c("Cd28", "Ifng", "Il2","Il2ra","Cd69","Cd44","Tnfrsf4", "Tnfrsf9", "Gzmb", "Ctla4","Pdcd1", "Havcr2", "Lag3")
tiff("Cd8 actexh individual.jpeg", unit="in", width=15, height=4, res=300)
DotPlot(cd8t, features = dot_features, group.by="orig.ident") + scale_size(range = c(1,20)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()

dot_features <- c("Atf1", "Atf2","Atf3","Atf4","Atf5","Atf6", "Atf7", "Batf", "Batf2", "Batf3", "Jdp2",
                  "Jun", "Junb", "Jund",
                  "Fos", "Fosb", "Fosl1", "Fosl2",
                  "Maf", "Mafa", "Mafb", "Maff", "Mafg", "Mafk")
tiff("cd8t Ap-1.jpeg", unit="in", width=15, height=4, res=300)
DotPlot(cd8t, features = dot_features, group.by="treatment") + scale_size(range = c(1,15)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()

#NFAT pathway (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3694398/)
#Ppp3cc (Calcineurin)
#NFAT1 = NFATc2
#NFAT2= NFATc1
#NFAT4= NFATc3
#NFAT5


dot_features <- c("Cracr2a","Nfatc1", "Nfatc2", "Ppp3ca","Ppp3cc")
tiff("Cd8 NFAT.jpeg", unit="in", width=8, height=5, res=300)
DotPlot(cd8t, features = dot_features, group.by="treatment") + scale_size(range = c(1,20)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()

#NF-kb pathway
#IKK-alpha = CHUK
#IKK-beta = Ikbkb
#NEMO = Ikbkg

dot_features <- c("Nfkb1", "Nfkb2",
                  "Chuk","Ikbkb", "Ikbkg", "Rel")

tiff("Cd8 NFKB.jpeg", unit="in", width=8, height=5, res=300)
DotPlot(cd8t, features = dot_features, group.by="treatment") + scale_size(range = c(1,20)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()

# Trac - TCR receptor alpha constant "Trac","Trbc1","Trbc2",
# Cd28 PI3k signaling
# Pik4ca - PI3k
# PRkcq PKC-theta
dot_features <- c("Cd28","Icos","Pik3ca", "Prkcq", "Mapk1",
                  "Plcg1","Grb2", "Vav1", "Card9")

tiff("PI3K CD28 signaling NfKB.jpeg", unit="in", width=8, height=5, res=300)
DotPlot(cd8t, features = dot_features, group.by="treatment") + scale_size(range = c(1,20)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()

#NKT-like CD8t cells
#Cd49b = Itga2
#NKG2D = Klrk1
#CD94 = Klrd1
#NK1.1 = Klrb1
#iNKT1 genes - T-bet(Tbx21), Ifng, 
#iNKT2 genes: Il4, PLZF(Zbtb16)
#iNKT17 genes: Il17, ROR??t(Rorc),
#CD69 for activated NK cells
dot_features <- c("Ncam1", "Itga2","Klrk1","Klrd1", "Klrb1", "Klrb1a","Klrb1b", "Klrb1c", "Klrb1d",
                  "Trbc1", "Trac",
                  "Tbx21", "Ifng",
                  "Il4", "Zbtb16",
                  "Il17a","Rorc",
                  "Cd69", "Gzmb")

#tiff("PI3K CD28 signaling NfKB.jpeg", unit="in", width=8, height=5, res=300)
DotPlot(cd8t, features = dot_features, group.by="treatment") + scale_size(range = c(1,20)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
#dev.off()
####################################################################################################################
#GSEA

#####################################################################################################################
gsea_TNK <- TNK
#Macrophages <- SetIdent(Macrophages, value = "Macrophages") #set all clusters to macrophages

gsea_TNK[["orig.ident"]] <- replace(gsea_TNK[["orig.ident"]], gsea_TNK[["orig.ident"]]=="C1", "control")
gsea_TNK[["orig.ident"]] <- replace(gsea_TNK[["orig.ident"]], gsea_TNK[["orig.ident"]]=="C2", "control")
gsea_TNK[["orig.ident"]] <- replace(gsea_TNK[["orig.ident"]], gsea_TNK[["orig.ident"]]=="T1", "TAK981")
gsea_TNK[["orig.ident"]] <- replace(gsea_TNK[["orig.ident"]], gsea_TNK[["orig.ident"]]=="T2", "TAK981")

Idents(gsea_TNK) <- "orig.ident" #set all idents to orig.ident i.e. control, TAK981, etc

TNK_ranked_list <- FindMarkers(gsea_TNK, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)

TNK_ranked_list <- TNK_ranked_list[order(TNK_ranked_list$avg_log2FC, decreasing = T),]
TNK_ranked_list$Gene.name <- str_to_upper(rownames(TNK_ranked_list))
TNK_ranked_list <- TNK_ranked_list[,c("Gene.name", "avg_log2FC")]
rownames(TNK_ranked_list) <- NULL

TNK_ranked_list <- prepare_ranked_list(TNK_ranked_list)

TNK_fgsea_results <- fgseaMultilevel(pathways = msigdb.all,
                                 stats = TNK_ranked_list,
                                 minSize = 15,
                                 maxSize = Inf)

TNK_fgsea_df<- TNK_fgsea_results %>% arrange (desc(NES)) %>% select (pathway, pval, padj, NES, size)
TNK_fgsea_df <- as.data.frame(TNK_fgsea_results)
TNK_fgsea_df = data.frame(lapply(TNK_fgsea_df, as.character), stringsAsFactors=FALSE)

write.csv(TNK_fgsea_df, "TNK_fgsea.csv",row.names = FALSE)




topPathwaysUp <- TNK_fgsea_results[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- TNK_fgsea_results[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(hallmark_pathway[topPathways],
              stats=TNK_ranked_list,
              fgseaRes=TNK_fgsea_results, 
              gseaParam = 0.5)


tiff("TNK_gseatable.jpeg", unit="in", width=20, height=5, res=300)
plotGseaTable(hallmark_pathway[topPathways],
              stats=TNK_ranked_list,
              fgseaRes=Macrophages_fgsea_results, 
              gseaParam = 0.5)
dev.off()
plot_enrichment <- function (geneset, pathway, ranked_list) {
  plotEnrichment(geneset[[pathway]], ranked_list)+labs (title = pathway)
}

# example of pathway highly enriched in treated

waterfall_plot <- function (fsgea_results, graph_title) {
  fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot(aes(reorder(short_name,NES), NES)) +
    geom_bar(stat= "identity", aes(fill = padj<0.05))+
    coord_flip()+
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
    theme(axis.text.y = element_text(size = 7), 
          plot.title = element_text(hjust = 0))
}
# example of pathway highly enriched in treated

waterfall_plot(fgsea_results, "Pathways enriched in TAK981-treated vs control T cells")

tiff("waterfall_control_vs_TAK981_Cd8.jpeg", unit="in", width=9, height=6, res=300)
waterfall_plot(fgsea_results, "Pathways enriched in TAK981-treated vs control CD8+ T cells")
dev.off()

tiff("cd8_gsea_TNFa.jpeg", unit="in", width=9, height=5, res=300)
plot_enrichment(hallmark_pathway, "HALLMARK_TNFA_SIGNALING_VIA_NFKB" , TNK_ranked_list)
dev.off()

# 
####################################################################################################################
#NKT cell analysis
NKT <- subset(TNK, subset=(SingleR.label=="NKT"))
DimPlot(NKT, reduction="tsne", label=TRUE,pt.size=3, group.by = "SingleR.label", split.by = "orig.ident")

dot_features <- c("Cd28", "Ifng","Il2ra", "Foxp3","Cd69","Cd44","Tnfrsf4", "Tnfrsf9", "Gzmb", "Ctla4","Pdcd1", "Havcr2", "Lag3")

#CD56 = Ncam1
#Cd49b = Itga2
#NKG2D = Klrk1
#CD94 = Klrd1
#NK1.1 = Klrb1
#iNKT1 genes - T-bet(Tbx21), Ifng, 
#iNKT2 genes: Il4, PLZF(Zbtb16)
#iNKT17 genes: Il17, ROR??t(Rorc),
#CD69 for activated NK cells

# dot_features <- c("Ncam1", "Itga2","Klrk1","Klrd1", "Klrb1f", "Cd44",
#                   "Trbc1", "Trac",
#                   "Tbx21", "Ifng",
#                   "Il4", "Zbtb16",
#                   "Il17a","Rorc",
#                   "Cd69", "Gzmb")
tiff("NKT actexh.jpeg", unit="in", width=15, height=5, res=300)
DotPlot(NKT, features = dot_features, group.by="treatment") + scale_size(range = c(1,15)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()
tiff("NKT actexh individual.jpeg", unit="in", width=15, height=5, res=300)
DotPlot(NKT, features = dot_features, group.by="orig.ident") + scale_size(range = c(1,15)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()
dot_features <- c("Atf1", "Atf2","Atf3","Atf4","Atf5","Atf6", "Atf7", "Batf", "Batf2", "Batf3", "Jdp2",
                  "Jun", "Junb", "Jund",
                  "Fos", "Fosb", "Fosl1", "Fosl2",
                  "Maf", "Mafa", "Mafb", "Maff", "Mafg", "Mafk")
tiff("NKT Ap-1.jpeg", unit="in", width=15, height=4, res=300)
DotPlot(cd8t, features = dot_features, group.by="treatment") + scale_size(range = c(1,15)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()

dot_features <- c("Cracr2a","Nfatc1", "Nfatc2", "Ppp3ca","Ppp3cc")
tiff("NKT NFAT.jpeg", unit="in", width=8, height=5, res=300)
DotPlot(cd8t, features = dot_features, group.by="treatment") + scale_size(range = c(1,20)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()

dot_features <- c("Nfkb1", "Nfkb2",
                  "Chuk","Ikbkb", "Ikbkg", "Rel")

tiff("NKT NFKB.jpeg", unit="in", width=8, height=5, res=300)
DotPlot(cd8t, features = dot_features, group.by="treatment") + scale_size(range = c(1,20)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()
#################################################################################################################
gsea_TNK <- NKT

gsea_TNK[["orig.ident"]] <- replace(gsea_TNK[["orig.ident"]], gsea_TNK[["orig.ident"]]=="C1", "control")
gsea_TNK[["orig.ident"]] <- replace(gsea_TNK[["orig.ident"]], gsea_TNK[["orig.ident"]]=="C2", "control")
gsea_TNK[["orig.ident"]] <- replace(gsea_TNK[["orig.ident"]], gsea_TNK[["orig.ident"]]=="T1", "TAK981")
gsea_TNK[["orig.ident"]] <- replace(gsea_TNK[["orig.ident"]], gsea_TNK[["orig.ident"]]=="T2", "TAK981")

Idents(gsea_TNK) <- "orig.ident" #set all idents to orig.ident i.e. control, TAK981, etc

TNK_ranked_list <- FindMarkers(gsea_TNK, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)

TNK_ranked_list <- TNK_ranked_list[order(TNK_ranked_list$avg_log2FC, decreasing = T),]
TNK_ranked_list$Gene.name <- str_to_upper(rownames(TNK_ranked_list))
TNK_ranked_list <- TNK_ranked_list[,c("Gene.name", "avg_log2FC")]
rownames(TNK_ranked_list) <- NULL

TNK_ranked_list <- prepare_ranked_list(TNK_ranked_list)

fgsea_results <- fgseaMultilevel(pathways = hallmark_pathway,
                                 stats = TNK_ranked_list,
                                 minSize = 15,
                                 maxSize = Inf)

fgsea_df<- fgsea_results %>% arrange (desc(NES)) %>% select (pathway, pval, padj, NES, size)
write.csv(fgsea_df, "fgsea.csv",row.names = FALSE)

# example of pathway highly enriched in treated
plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_ALPHA_RESPONSE" , TNK_ranked_list)

waterfall_plot(fgsea_results, "Pathways enriched in TAK981-treated vs control NKT cells")

tiff("waterfall_control_vs_TAK981_NKT.jpeg", unit="in", width=9, height=6, res=300)
waterfall_plot(fgsea_results, "Pathways enriched in TAK981-treated vs control NKT cells")
dev.off()


tiff("NKT_gsea_TNFa.jpeg", unit="in", width=5, height=3, res=300)
plot_enrichment(hallmark_pathway, "HALLMARK_TNFA_SIGNALING_VIA_NFKB" , TNK_ranked_list)
dev.off()


#################################################################################################################
#################################################################################################################
Blymph <- subset(TNK, subset=(SingleR.label=="B cells"))

#tiff("Blypmh_tsne.jpeg", unit="in", width=15, height=6, res=300)
DimPlot(Blymph, reduction="tsne",label=TRUE,pt.size=3, split.by = "treatment")
#dev.off()

Blymph <- NormalizeData(Blymph, normalization.method = "LogNormalize", scale.factor = 10000)
Blymph <- FindVariableFeatures(Blymph, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Blymph)
Blymph <- ScaleData(Blymph, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo"))
Blymph <- RunPCA (Blymph, features = VariableFeatures(object = Macrophages))
Blymph <- FindNeighbors(Blymph, dims = 1:40)
Blymph <- FindClusters(Blymph, resolution = c(0.6))
Blymph <- RunTSNE(Blymph, dims = 1:40)
Blymph <- RunUMAP(Blymph, dims = 1:40)
Blymph.markers <- FindAllMarkers(Blymph, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Blymph.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> Blymph.top10


SRBlymph <- as.SingleCellExperiment(Blymph)
SRBlymph <- SingleR(test=SRBlymph, ref=mref, assay.type.test = 1, labels = mref$label.main)
Blymph[["SingleR.label"]] <- SRBlymph$labels
rm(SRBlymph)
DimPlot(Blymph, reduction="umap", label=TRUE,pt.size=3, group.by = "SingleR.label", split.by = "treatment")


#naive b cell: Cd20+, Cd38-, CD27-, IgD+
#follicular b cell: cd1d low, Cd23(Fcer2) high, Cd24 low, 
#marginal zone b cell: Slc22a2, Cd1d high, Cd21 high, IgM high,
#memory b cell: Obf1(Pou2af1)+, cd21(cr2)+, PD-L2(PDCD1LG2)+
#plasma cells: Cd27 high, Cd138(Sdc1) high, IRF4+, Ly6k
#Breg: CD24+,, Tgfb1+, IgM high
dot_features <- c("Fcer2a","Ighd", "H2-Ab1",
                  "Cd1d1", "Cr2","Ighm",
                  "Pou2af1", "Pdcd1lg2", 
                  "Cd27","Sdc1","Prdm1", "Irf4", "Ly6k",
                  "Cd24a", "Tgfb1")
tiff("Blymph.jpeg", unit="in", width=12, height=5, res=300)
DotPlot(Blymph, features = dot_features, group.by="treatment") + scale_size(range = c(1,15)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()

dot_features <- c("Nfkbia","Nfkbib", "Nfkbie", "Chuk","Ikbkb", "Ikbkg",
                  "Nfkb1", "Rel", "Rela", "Nfkb2", "Relb")
tiff("Blymph_NFkB.jpeg", unit="in", width=12, height=5, res=300)
DotPlot(Blymph, features = dot_features, group.by="treatment") + scale_size(range = c(1,15)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()



gsea_Blymph <- Blymph

gsea_Blymph[["orig.ident"]] <- replace(gsea_Blymph[["orig.ident"]], gsea_Blymph[["orig.ident"]]=="C1", "control")
gsea_Blymph[["orig.ident"]] <- replace(gsea_Blymph[["orig.ident"]], gsea_Blymph[["orig.ident"]]=="C2", "control")
gsea_Blymph[["orig.ident"]] <- replace(gsea_Blymph[["orig.ident"]], gsea_Blymph[["orig.ident"]]=="T1", "TAK981")
gsea_Blymph[["orig.ident"]] <- replace(gsea_Blymph[["orig.ident"]], gsea_Blymph[["orig.ident"]]=="T2", "TAK981")

Idents(gsea_Blymph) <- "orig.ident" #set all idents to orig.ident i.e. control, TAK981, etc

Blymph_ranked_list <- FindMarkers(gsea_Blymph, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)

Blymph_ranked_list <- Blymph_ranked_list[order(Blymph_ranked_list$avg_log2FC, decreasing = T),]
Blymph_ranked_list$Gene.name <- str_to_upper(rownames(Blymph_ranked_list))
Blymph_ranked_list <- Blymph_ranked_list[,c("Gene.name", "avg_log2FC")]
rownames(Blymph_ranked_list) <- NULL

Blymph_ranked_list <- prepare_ranked_list(Blymph_ranked_list)

fgsea_results <- fgseaMultilevel(pathways = hallmark_pathway,
                                 stats = Blymph_ranked_list,
                                 minSize = 15,
                                 maxSize = Inf)

fgsea_df<- fgsea_results %>% arrange (desc(NES)) %>% select (pathway, pval, padj, NES, size)
write.csv(fgsea_df, "fgsea.csv",row.names = FALSE)

waterfall_plot(fgsea_results, "Pathways enriched in TAK981-treated vs control B cells")

tiff("waterfall_control_vs_TAK981_B_cells.jpeg", unit="in", width=9, height=6, res=300)
waterfall_plot(fgsea_results, "Pathways enriched in TAK981-treated vs control B cells")
dev.off()

Blymph_ranked_list <- FindMarkers(gsea_Blymph, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)

Blymph_ranked_list <- Blymph_ranked_list[order(Blymph_ranked_list$avg_log2FC, decreasing = T),]
Blymph_ranked_list$Gene.name <- str_to_upper(rownames(Blymph_ranked_list))
Blymph_ranked_list <- Blymph_ranked_list[,c("Gene.name", "avg_log2FC")]
rownames(Blymph_ranked_list) <- NULL

Blymph_ranked_list <- prepare_ranked_list(Blymph_ranked_list)

fgsea_results <- fgseaMultilevel(pathways = hallmark_pathway,
                                 stats = Blymph_ranked_list,
                                 minSize = 15,
                                 maxSize = 500)

fgsea_results %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()

waterfall_plot(fgsea_results, "Pathways enriched in TAK981-treated vs control B cells")

tiff("waterfall_control_vs_TAK981_B_cells.jpeg", unit="in", width=9, height=6, res=300)
waterfall_plot(fgsea_results, "Pathways enriched in TAK981-treated vs control B cells")
dev.off()

tiff("B_gsea_apoptosis.jpeg", unit="in", width=5, height=3, res=300)
plot_enrichment(hallmark_pathway, "HALLMARK_APOPTOSIS" , TNK_ranked_list)
dev.off()
############################################################################################################

############################################################################################################
#nonimmune cell analysis
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

mref <- ImmGenData() #mouse immune cells
sceM <- BaronPancreasData('mouse')
sceM <- sceM[,!is.na(sceM$label)]
library(scuttle)
sceM <- logNormCounts(sceM)
## import all necessary files
## use the following command to unzip tar.gz: tar -xvzf
setwd("E:/TAK981_KPC/control_vs_TAK981")
hallmark_pathway <- gmtPathways("h.all.v7.4.symbols.gmt.txt")
load('nonimmune.RData')
original <- nonimmune
#perform scTransform
scKPC <- nonimmune
scKPC <- DietSeurat(scKPC, assay="RNA")
scKPC.list <- SplitObject(scKPC, split.by = "orig.ident")

for (i in 1:length(scKPC.list)) {
  scKPC.list[[i]] <- SCTransform(scKPC.list[[i]], vars.to.regress=c("nCount_RNA","percent.mito", "percent.ribo"))
}
features <- SelectIntegrationFeatures(object.list = scKPC.list, nfeatures = 2000)
scKPC.list <- PrepSCTIntegration(object.list = scKPC.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = scKPC.list, normalization.method = "SCT",
                                  anchor.features = features)
memory.limit(56000)
scKPC <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
#DefaultAssay(scKPC) <- "integrated"
scKPC <- ScaleData(object = scKPC, features = features)
scKPC <- RunPCA(scKPC, verbose = FALSE)
scKPC <- RunTSNE(scKPC, dims = 1:40, verbose = FALSE)
scKPC <- RunUMAP(scKPC, dims = 1:40, verbose = FALSE)
scKPC <- FindNeighbors(scKPC, dims = 1:40, verbose = FALSE)
scKPC <- FindClusters(scKPC, resolution = 1.4)
DimPlot(scKPC, reduction="tsne",label=TRUE,pt.size=1)

#Must normalize data with Lognormalize for further DE analyses
DefaultAssay(scKPC) <- "RNA"
scKPC <- NormalizeData(scKPC)
all.genes <- rownames(scKPC)
scKPC <- ScaleData(object = scKPC, features = all.genes)
scKPC.markers <- FindAllMarkers(scKPC, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scKPC.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10
write.csv(scKPC.markers, "nonimmune_markers.csv", row.names = T)

tiff("nonimmune_cluster_heatmap.jpeg", unit="in", width=30, height=18, res=500)
DoHeatmap(scKPC, features = top10$gene) + NoLegend()
dev.off()

DimPlot(scKPC, reduction="tsne",label=TRUE,pt.size=1)

tiff("nonimmune_tsne_unlabeled.jpeg", unit="in", width=6, height=5, res=300)
DimPlot(scKPC, reduction="tsne",label=TRUE,pt.size=1, repel = T)
dev.off()

DimPlot(scKPC, reduction="tsne",label=TRUE,pt.size=1, group.by = "SingleR.label")

temp <- scKPC
SRKPC <- as.SingleCellExperiment(temp)
SRKPC <- SingleR(test=SRKPC, ref=sceM, assay.type.test = 1, labels = sceM$label)
temp[["SingleR.label"]] <- SRKPC$labels
rm(SRKPC)
DimPlot(temp, reduction="tsne",label=TRUE,pt.size=1, group.by = "SingleR.label")

nonimmune <- scKPC
VlnPlot(nonimmune, features = c("Ptprc"))





nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==8, "myCAF")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==15, "iCAF")  
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==19, "myCAF-like PSC") #Acta2 (a-SMA) increase  

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==12, "Endothelial")

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==14, "ADM") #Acinar-to-ductal metaplasia area

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==13, "Ductal 0") #Ptprc+  Epcam+
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==17, "Ductal 0") #Ptprc+  Epcam+

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==18, "Ductal NM") # EpCAM+

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==0, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==1, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==3, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==5, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==6, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==7, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==10, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==11, "mt") #mt and ribo genes. exclude

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==2, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==4, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==9, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==16, "Ductal 2")

kpc <- subset(nonimmune, subset=(SingleR.label=="ADM" | SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2" | 
                                         SingleR.label=="Ductal 0" | SingleR.label=="Ductal NM" | SingleR.label=="Endothelial" |
                                         SingleR.label=="iCAF"| SingleR.label=="myCAF-like PSC" | SingleR.label=="myCAF"))

DimPlot(kpc, reduction="tsne",label=TRUE,pt.size=1, group.by = "SingleR.label")

tiff("nonimmune_cluster_heatmap_SingleR.jpeg", unit="in", width=35, height=18, res=500)
DoHeatmap(kpc, features = top10$gene, group.by = "SingleR.label") + NoLegend()
dev.off()


tiff("nonimmune_tsne.jpeg", unit="in", width=6, height=5, res=300)
DimPlot(kpc, reduction="tsne", label=TRUE,pt.size=1, group.by = "SingleR.label",
        repel = TRUE)
dev.off()

tiff("nonimmune_tsne_individual.jpeg", unit="in", width=15, height=5, res=300)
DimPlot(kpc, reduction="tsne", label=TRUE,pt.size=1, group.by = "SingleR.label", split.by = "treatment", 
        repel = TRUE)
dev.off()

###########################################################################################################################
ductal <- subset(kpc, subset=(SingleR.label=="ADM" |SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2" | SingleR.label=="Ductal 0" | SingleR.label=="Ductal NM"))
CAF <- subset(kpc, subset=(SingleR.label=="myCAF" | SingleR.label=="iCAF" | SingleR.label=="myCAF-like PSC"))

tiff("kpc_Krt19.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(kpc, features = c("Krt19"), group.by = "SingleR.label")
dev.off()

tiff("kpc_Tspan8.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(kpc, features = c("Tspan8"), group.by = "SingleR.label")
dev.off()

tiff("kpc_Sox9.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(kpc, features = c("Sox9"), group.by = "SingleR.label")
dev.off()

tiff("kpc_Pecam1.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(kpc, features = c("Pecam1"), group.by = "SingleR.label")
dev.off()

tiff("kpc_Prss3.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(kpc, features = c("Prss3"), group.by = "SingleR.label")
dev.off()

tiff("kpc_Rgs5.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(kpc, features = c("Rgs5"), group.by = "SingleR.label")
dev.off()

tiff("kpc_Col1a1.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(kpc, features = c("Col1a1"), group.by = "SingleR.label")
dev.off()

################################################################### 

#common fibroblast marker
tiff("kpc_Col1a1.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(kpc, features = c("Col1a1"),group.by = "SingleR.label")
dev.off()
#common CAF marker https://cancerdiscovery.aacrjournals.org/content/9/8/1102
tiff("kpc_Dcn.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(kpc, features = c("Onecut2"),group.by = "SingleR.label") #DCN is higher in iCAF
dev.off()

Prss3 <- VlnPlot(kpc, features = c("Prss3"), group.by = "SingleR.label") +
  ylab("Prss3") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Krt19 <- VlnPlot(kpc, features = c("Krt19"), group.by = "SingleR.label") +
  ylab("Krt19") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Epcam <- VlnPlot(kpc, features = c("Epcam"), group.by = "SingleR.label") +
  ylab("Epcam") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Pecam1 <- VlnPlot(kpc, features = c("Pecam1"), group.by = "SingleR.label") +
  ylab("Pecam1") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Col1a1 <- VlnPlot(kpc, features = c("Col1a1"), group.by = "SingleR.label") +
  ylab("Col1a1") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Dcn <- VlnPlot(kpc, features = c("Dcn"), group.by = "SingleR.label") +
  ylab("Dcn") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Acta2 <- VlnPlot(kpc, features = c("Acta2"), group.by = "SingleR.label") +
  ylab("Acta2") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Postn <- VlnPlot(kpc, features = c("Postn"), group.by = "SingleR.label") +
  ylab("Postn") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Rgs5 <- VlnPlot(kpc, features = c("Rgs5"), group.by = "SingleR.label") +
  ylab("Rgs5") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_text(size = 25),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

tiff("kpc_vlnplot.jpeg", unit="in", width=24, height=20, res=500)
Prss3/Krt19/Epcam/Pecam1/Col1a1/Dcn/Acta2/Postn/Rgs5
dev.off()

#myCAF - indced through TGFb/SMAD2/3 
tiff("CAF_Acta2(a-SMA).jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(CAF, features = c("Acta2"),group.by = "SingleR.label") #a-SMA for activated myCAFs
dev.off()
tiff("CAF_Postn.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(CAF, features = c("Postn"),group.by = "SingleR.label") #iCAF show high Il6 but lower aSMA
dev.off()
tiff("CAF_Tagln.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(CAF, features = c("Tagln"),group.by = "SingleR.label") #iCAF show high Il6 but lower aSMA
dev.off()

#iCAF - induced through IL1/JAK-STAT3
tiff("CAF_Il6.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(CAF, features = c("Il6"),group.by = "SingleR.label") #iCAF show high Il6 but lower aSMA
dev.off()
tiff("CAF_Dpt.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(CAF, features = c("Dpt"),group.by = "SingleR.label") #matrix protein
dev.off()
tiff("CAF_Has1.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(CAF, features = c("Has1"),group.by = "SingleR.label")
dev.off()
tiff("CAF_Pdgfra.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(CAF, features = c("Pdgfra"),group.by = "SingleR.label")
dev.off()

tiff("CAF_S100a4.jpeg", unit="in", width=5, height=4, res=300)
VlnPlot(CAF, features = c("S100a4"),group.by = "SingleR.label",split.by="treatment") #iCAF show high Il6 but lower aSMA
dev.off()



##################################################################################################################################### 

FeaturePlot(ductal, features=c("Foxq1"),pt.size=1.5,reduction="tsne", split.by = "treatment")
VlnPlot(nonimmune, features = c("Cd47"), group.by = "SingleR.label",  split.by = "treatment")
VlnPlot(ductal, features = c("Aldh1a1", "Pou2f1", "Nes"), group.by = "SingleR.label", split.by = "treatment")
# Epithelial vs Mesenchymal markers 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6777805/#sd
# epithelial markers: Ocln, Gjb1, and Tjp1, Cldn4, Cldn3, Epcam

# cancer stem cell markers: Cldn3, Epcam, Cdh1 (E-cadherin)
tiff("ductal_Aldh1a1.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(ductal, features = c("Aldh1a1"), group.by = "SingleR.label", split.by = "treatment")
dev.off()
tiff("ductal_Pou2f1.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(ductal, features = c("Pou2f1"), group.by = "SingleR.label", split.by = "treatment")
dev.off()
tiff("ductal_Nes.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(ductal, features = c("Nes"), group.by = "SingleR.label", split.by = "treatment")
dev.off()

# Epithelial markers: Cldn3, Epcam, Cdh1 (E-cadherin)
tiff("ductal_Cldn3.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(ductal, features = c("Cldn3"), group.by = "SingleR.label", split.by = "treatment")
dev.off()
tiff("ductal_Epcam.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(ductal, features = c("Epcam"), group.by = "SingleR.label", split.by = "treatment")
dev.off()
tiff("ductal_Cdh1.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(ductal, features = c("Cdh1"), group.by = "SingleR.label", split.by = "treatment")
dev.off()


# mesenchymal markers: S100a4, Vim, Cdh2 (N-cadherin)
tiff("ductal_S100a4.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(ductal, features = c("S100a4"), group.by = "SingleR.label", split.by = "treatment")
dev.off()
tiff("ductal_Vim.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(ductal, features = c("Vim"), group.by = "SingleR.label", split.by = "treatment")
dev.off()
tiff("ductal_Cdh2.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(ductal, features = c("Cdh2"), group.by = "SingleR.label", split.by = "treatment")
dev.off()


#Proliferation markers: Mki67, Top2a
tiff("ductal_Mki67.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(ductal, features = c("Mki67"), group.by = "SingleR.label", split.by = "treatment")
dev.off()
tiff("ductal_Top2a.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(ductal, features = c("Top2a"), group.by = "SingleR.label", split.by = "treatment")
dev.off()

#Cell Cycle marker
tiff("ductal_Ccnb2.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(ductal, features = c("Ccnb2"), group.by = "SingleR.label", split.by = "treatment")
dev.off()


########################################
#tiff("kpc.jpeg", unit="in", width=15, height=5, res=300)
DotPlot(kpc, features = dot_features, group.by="treatment") + scale_size(range = c(1,15)) + RotatedAxis() + 
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
#dev.off()


VlnPlot(ductal, features = c("Setdb1"), group.by = "SingleR.label", split.by = "treatment")

dot_features <- c("Setdb1","Ifnar1", "Ifnar2", "Ifnb1", "Ifit1", "Ifit2", "Ifit3")
#dot_features <- c("Atf1", "Atf2","Atf3","Atf4","Atf5","Atf6", "Atf7", "Batf", "Batf2", "Batf3", "Jdp2",
#                  "Jun", "Junb", "Jund",
 #                 "Fos", "Fosb", "Fosl1", "Fosl2",
  #                "Maf", "Mafa", "Mafb", "Maff", "Mafg", "Mafk")
tiff("ductal_Setdb1Ifn.jpeg", unit="in", width=8, height=4, res=300)
DotPlot(ductal, features = dot_features, group.by="orig.ident") + scale_size(range = c(1,15)) + RotatedAxis() + 
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()
####################################################################################################################
# Group-wise cell proportion plot (https://erilu.github.io/single-cell-rnaseq-analysis/#group-wise_analysis)
source("scFunctions.R")

tiff("nonimmune_celltype_proportion.jpeg", unit="in", width=9, height=6, res=300)
plot_group_proportions(kpc, graph.type = "dodge")
dev.off()

tiff("nonimmune_celltype_proportion_stacked.jpeg", unit="in", width=3.5, height=6, res=300)
plot_group_proportions(kpc, graph.type = "stacked")
dev.off()

tiff("immune_proportion_heatmap.jpeg", unit="in", width=8, height=4, res=300)
plot_heatmap_proportions(kpc, graph.type = "by.cell")
dev.off()

####################################################################################################################
#GSEA
kpc.subset <- ductal
kpc.subset[["orig.ident"]] <- replace(kpc.subset[["orig.ident"]], kpc.subset[["orig.ident"]]=="C1", "control")
kpc.subset[["orig.ident"]] <- replace(kpc.subset[["orig.ident"]], kpc.subset[["orig.ident"]]=="C2", "control")
kpc.subset[["orig.ident"]] <- replace(kpc.subset[["orig.ident"]], kpc.subset[["orig.ident"]]=="T1", "TAK981")
kpc.subset[["orig.ident"]] <- replace(kpc.subset[["orig.ident"]], kpc.subset[["orig.ident"]]=="T2", "TAK981")
Idents(kpc.subset) <- "orig.ident" #set all idents to orig.ident i.e. control, TAK981, etc

kpc_ranked_list <- FindMarkers(kpc.subset, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)
temp_list <- kpc_ranked_list
hallmark_pathway <- gmtPathways("h.all.v7.4.symbols.gmt.txt")
# order list, pull out gene name and log2fc, and convert genes to uppercase
kpc_ranked_list <- kpc_ranked_list[order(kpc_ranked_list$avg_log2FC, decreasing = T),]
kpc_ranked_list$Gene.name <- str_to_upper(rownames(kpc_ranked_list))
kpc_ranked_list <- kpc_ranked_list[,c("Gene.name", "avg_log2FC")]
rownames(kpc_ranked_list) <- NULL
head(kpc_ranked_list)

kpc_ranked_list <- prepare_ranked_list(kpc_ranked_list)
head(kpc_ranked_list)

fgsea_results <- fgseaMultilevel(pathways = hallmark_pathway,
                                 stats = kpc_ranked_list,
                                 minSize = 15,
                                 maxSize = Inf)

fgsea_df<- fgsea_results %>% arrange (desc(NES)) %>% select (pathway, pval, padj, NES, size)
write.csv(fgsea_df, "ductal_fgsea.csv",row.names = FALSE)

# example of pathway highly enriched in treated
tiff("ductal_TNFa.jpeg", unit="in", width=5, height=3, res=300)
plot_enrichment(hallmark_pathway, "HALLMARK_TNFA_SIGNALING_VIA_NFKB" , kpc_ranked_list)
dev.off()

tiff("ductal_apoptosis.jpeg", unit="in", width=5, height=3, res=300)
plot_enrichment(hallmark_pathway, "HALLMARK_APOPTOSIS" , kpc_ranked_list)
dev.off()
tiff("ductal_oxidative.jpeg", unit="in", width=5, height=3, res=300)
plot_enrichment(hallmark_pathway, "HALLMARK_OXIDATIVE_PHOSPHORYLATION" , kpc_ranked_list)
dev.off()

tiff("ductal_gsea.jpeg", unit="in", width=10, height=7, res=300)
waterfall_plot(fgsea_results, "Pathways enriched in TAK981-treated vs control ductal cells")
dev.off()


#############################################################################################################
library("org.Mm.eg.db", character.only = TRUE)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(ggnewscale)
organism <-org.Mm.eg.db

ductal_ranked_list <- FindMarkers(ductal, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)
# order list, pull out gene name and log2fc, and convert genes to uppercase
ductal_ranked_list <- ductal_ranked_list[order(ductal_ranked_list$avg_log2FC, decreasing = T),]
ductal_ranked_list$Gene.name <- rownames(ductal_ranked_list)
ductal_ranked_list <- ductal_ranked_list[,c("Gene.name", "avg_log2FC")]
rownames(ductal_ranked_list) <- NULL

ductal_ranked_list2 <- prepare_ranked_list(ductal_ranked_list)

gse <- gseGO(geneList=ductal_ranked_list2, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 100, 
             maxGSSize = 500, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
goplot(gse)
require(DOSE)
enrichplot::dotplot(gse, showCategory=15, split=".sign", title="GO Biological Process Enrichment - KPC TAK981 treatment", font.size=10) + facet_grid(.~.sign) +
  scale_y_discrete(labels=function(gse) str_wrap(gse, width=20))
enrichplot::dotplot(gse, showCategory=20)

tiff("KPC_gseGO.jpeg", unit="in", width=10, height=12, res=300)
enrichplot::dotplot(gse, showCategory=15, split=".sign", title="GO Biological Process Enrichment - KPC TAK981 treatment", font.size=8) + facet_grid(.~.sign) +
  scale_y_discrete(labels=function(gse) str_wrap(gse, width=40))
dev.off()

#############################################################################################################
######### ductal cell (cancer cell) only GSEA
ductal <- subset(kpc, subset=(SingleR.label=="ADM" |SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2" | SingleR.label=="Ductal 0" | SingleR.label=="Ductal NM"))

ductal[["orig.ident"]] <- replace(ductal[["orig.ident"]], ductal[["orig.ident"]]=="C1", "control")
ductal[["orig.ident"]] <- replace(ductal[["orig.ident"]], ductal[["orig.ident"]]=="C2", "control")
ductal[["orig.ident"]] <- replace(ductal[["orig.ident"]], ductal[["orig.ident"]]=="T1", "TAK981")
ductal[["orig.ident"]] <- replace(ductal[["orig.ident"]], ductal[["orig.ident"]]=="T2", "TAK981")
Idents(ductal) <- "orig.ident" #set all idents to orig.ident i.e. control, TAK981, etc

ductal_ranked_list <- FindMarkers(ductal, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)
# order list, pull out gene name and log2fc, and convert genes to uppercase
ductal_ranked_list <- ductal_ranked_list[order(ductal_ranked_list$avg_log2FC, decreasing = T),]
ductal_ranked_list$Gene.name <- str_to_upper(rownames(ductal_ranked_list))
ductal_ranked_list <- ductal_ranked_list[,c("Gene.name", "avg_log2FC")]
rownames(ductal_ranked_list) <- NULL
ductal_ranked_list <- prepare_ranked_list(ductal_ranked_list)
head(ductal_ranked_list)

ductal_fgsea_results <- fgseaMultilevel(pathways = hallmark_pathway,
                                 stats = ductal_ranked_list,
                                 minSize = 15,
                                 maxSize = Inf)

ductal_fgsea_df<- ductal_fgsea_results %>% arrange (desc(NES)) %>% select (pathway, pval, padj, NES, size)
write.csv(ductal_fgsea_df, "ductal_fgsea.csv",row.names = FALSE)

waterfall_plot(ductal_fgsea_results, "Pathways enriched in TAK981-treated vs control ductal cells")
# example of pathway highly enriched in treated
plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_ALPHA_RESPONSE" , ductal_ranked_list)
plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_GAMMA_RESPONSE" , ductal_ranked_list)

tiff("ductal_gsea.jpeg", unit="in", width=10, height=7, res=300)
waterfall_plot(ductal_fgsea_results, "Pathways enriched in TAK981-treated vs control ductal cells")
dev.off()


################## fibroblasts only GSEA
fibroblasts <- subset(kpc, subset=(SingleR.label == "fibroblasts"))

fibroblasts[["orig.ident"]] <- replace(fibroblasts[["orig.ident"]], fibroblasts[["orig.ident"]]=="C1", "control")
fibroblasts[["orig.ident"]] <- replace(fibroblasts[["orig.ident"]], fibroblasts[["orig.ident"]]=="C2", "control")
fibroblasts[["orig.ident"]] <- replace(fibroblasts[["orig.ident"]], fibroblasts[["orig.ident"]]=="T1", "TAK981")
fibroblasts[["orig.ident"]] <- replace(fibroblasts[["orig.ident"]], fibroblasts[["orig.ident"]]=="T2", "TAK981")
Idents(fibroblasts) <- "orig.ident" #set all idents to orig.ident i.e. control, TAK981, etc

fibroblasts_ranked_list <- FindMarkers(fibroblasts, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)
# order list, pull out gene name and log2fc, and convert genes to uppercase
fibroblasts_ranked_list <- fibroblasts_ranked_list[order(fibroblasts_ranked_list$avg_log2FC, decreasing = T),]
fibroblasts_ranked_list$Gene.name <- str_to_upper(rownames(fibroblasts_ranked_list))
fibroblasts_ranked_list <- fibroblasts_ranked_list[,c("Gene.name", "avg_log2FC")]
rownames(fibroblasts_ranked_list) <- NULL
fibroblasts_ranked_list <- prepare_ranked_list(fibroblasts_ranked_list)
head(fibroblasts_ranked_list)

fibroblasts_fgsea_results <- fgseaMultilevel(pathways = hallmark_pathway,
                                        stats = fibroblasts_ranked_list,
                                        minSize = 15,
                                        maxSize = Inf)

fibroblasts_fgsea_df<- fibroblasts_fgsea_results %>% arrange (desc(NES)) %>% select (pathway, pval, padj, NES, size)
write.csv(fibroblasts_fgsea_df, "fibroblasts_fgsea.csv",row.names = FALSE)

waterfall_plot(fibroblasts_fgsea_results, "Pathways enriched in TAK981-treated vs control fibroblasts cells")
# example of pathway highly enriched in treated
plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_ALPHA_RESPONSE" , fibroblasts_ranked_list)
plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_GAMMA_RESPONSE" , fibroblasts_ranked_list)

tiff("fibroblasts_gsea.jpeg", unit="in", width=10, height=7, res=300)
waterfall_plot(fibroblasts_fgsea_results, "Pathways enriched in TAK981-treated vs control fibroblasts cells")
dev.off()


################## endothelial only GSEA
endothelial <- subset(kpc, subset=(SingleR.label == "endothelial"))

endothelial[["orig.ident"]] <- replace(endothelial[["orig.ident"]], endothelial[["orig.ident"]]=="C1", "control")
endothelial[["orig.ident"]] <- replace(endothelial[["orig.ident"]], endothelial[["orig.ident"]]=="C2", "control")
endothelial[["orig.ident"]] <- replace(endothelial[["orig.ident"]], endothelial[["orig.ident"]]=="T1", "TAK981")
endothelial[["orig.ident"]] <- replace(endothelial[["orig.ident"]], endothelial[["orig.ident"]]=="T2", "TAK981")
Idents(endothelial) <- "orig.ident" #set all idents to orig.ident i.e. control, TAK981, etc

endothelial_ranked_list <- FindMarkers(endothelial, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)
# order list, pull out gene name and log2fc, and convert genes to uppercase
endothelial_ranked_list <- endothelial_ranked_list[order(endothelial_ranked_list$avg_log2FC, decreasing = T),]
endothelial_ranked_list$Gene.name <- str_to_upper(rownames(endothelial_ranked_list))
endothelial_ranked_list <- endothelial_ranked_list[,c("Gene.name", "avg_log2FC")]
rownames(endothelial_ranked_list) <- NULL
endothelial_ranked_list <- prepare_ranked_list(endothelial_ranked_list)
head(endothelial_ranked_list)

endothelial_fgsea_results <- fgseaMultilevel(pathways = hallmark_pathway,
                                             stats = endothelial_ranked_list,
                                             minSize = 15,
                                             maxSize = Inf)

endothelial_fgsea_df<- endothelial_fgsea_results %>% arrange (desc(NES)) %>% select (pathway, pval, padj, NES, size)
write.csv(endothelial_fgsea_df, "endothelial_fgsea.csv",row.names = FALSE)

waterfall_plot(endothelial_fgsea_results, "Pathways enriched in TAK981-treated vs control endothelial cells")
# example of pathway highly enriched in treated
plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_ALPHA_RESPONSE" , endothelial_ranked_list)
plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_GAMMA_RESPONSE" , endothelial_ranked_list)

tiff("endothelial_gsea.jpeg", unit="in", width=10, height=7, res=300)
waterfall_plot(endothelial_fgsea_results, "Pathways enriched in TAK981-treated vs control endothelial cells")
dev.off()



#############################################################################################################################
# CellChat
load("immune.RData")
load("nonimmune.RData")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

#Set ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

memory.limit(56000)
scKPC <- merge(immune, y=c(nonimmune))
scKPC <- DietSeurat(scKPC, assay="RNA")
scKPC <- NormalizeData(scKPC,normalization.method = "LogNormalize", scale.factor = 10000)
scKPC <- FindVariableFeatures(scKPC, selection.method = "vst", nfeatures = 2000)
scKPC <- ScaleData(object = scKPC, verbose = FALSE)
scKPC <- RunPCA(scKPC, verbose = FALSE)
scKPC <- RunUMAP(scKPC, dims = 1:40, verbose = FALSE)
scKPC <- FindNeighbors(scKPC, verbose = FALSE)
scKPC <- FindClusters(scKPC, resolution = 0.4)
scKPC <- subset(scKPC, subset=(SingleR.label=="Ductal 0" | SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2" | 
                                 SingleR.label=="Ductal NM" | SingleR.label=="Endothelial" | SingleR.label=="Fibroblasts" | SingleR.label== "ADM" |
                                 SingleR.label== "iCAF" | SingleR.label== "myCAF" | SingleR.label== "myCAF-like PSC" | SingleR.label== "Endothelial" |
                                 SingleR.label=="Macro-Spp1" |SingleR.label=="Macro-Ly6c+Chil3+" | SingleR.label=="Macro-C1q" | 
                                 SingleR.label=="Macro-Ly6c+Isg15+" |SingleR.label=="Macro-Proliferating" | SingleR.label=="pDC" | 
                                 SingleR.label=="cDC1-Clec9a" | SingleR.label=="cDC1-Ccl22" | SingleR.label=="cDC2-Cd209a" |SingleR.label== "cDC2-Itgax"|
                                 SingleR.label=="T/NK cells" | SingleR.label=="Plasma cells" | SingleR.label=="Granulocytes" | SingleR.label=="B cells"))

DimPlot(scKPC, reduction="umap",label=TRUE,pt.size=1, repel = T)
DimPlot(scKPC, reduction="umap",label=TRUE,pt.size=1,group.by = "SingleR.label", repel = T)

treatment_split <- SplitObject(scKPC, split.by = "treatment")
control <- treatment_split$control
TAK981 <- treatment_split$TAK981

# Part I: Data input & processing and initialization of CellChat object
control_matrix <- GetAssayData(control, assay = "RNA", slot = "data") # normalized data matrix
Idents(control) <- "SingleR.label"
control_labels <- Idents(control)
control_meta <- data.frame(group = control_labels, row.names = names(control_labels)) # create a dataframe of the cell labels

control_cellchat <- createCellChat(object = control_matrix, meta = control_meta, group.by = "group")
control_cellchat <- addMeta(control_cellchat, meta = control_meta, meta.name = "labels")
control_cellchat <- setIdent(control_cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(control_cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(control_cellchat@idents)) # number of cells in each cell group

#Preprocessing the expression data for cell-cell communication analysis
control_cellchat@DB <- CellChatDB
control_cellchat <- subsetData(control_cellchat) # This step is necessary even if using the whole database
# future::plan("multiprocess", workers = 4)

control_cellchat <- identifyOverExpressedGenes(control_cellchat)
control_cellchat <- identifyOverExpressedInteractions(control_cellchat)
control_cellchat <- projectData(control_cellchat, PPI.mouse)


# Part II: Inference of cell-cell communication network
# Compute the communication probability and infer cellular communication network
control_cellchat <- computeCommunProb(control_cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
control_cellchat <- filterCommunication(control_cellchat, min.cells = 40)

# Extract the inferred cellular communication network as a data frame

# Option 1
# returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
# Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
df.net <- subsetCommunication(control_cellchat)

# # Option 2
# # gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# df.net <- subsetCommunication(control_cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
# 
# # Option 3
# # gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
# df.net <- subsetCommunication(control_cellchat, signaling = c("WNT", "TGFb")) 


# Infer the cell-cell communication at a signaling pathway level
control_cellchat <- computeCommunProbPathway(control_cellchat)

# Calculate the aggregated cell-cell communication network
# USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.
control_cellchat <- aggregateNet(control_cellchat)
groupSize <- as.numeric(table(control_cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(control_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(control_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# examine the signaling sent from each cell group
mat <- control_cellchat@net$weight
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
levels(control_cellchat@idents)
vertex.receiver = c(13,14,15,16,19,20,21)
netVisual_aggregate(control_cellchat, signaling = pathways.show,  vertex.receiver = c(13,14,15,16,19,20,21),layout = "hierarchy")

# Circle plot
pathways.show <- c("VEGF") #VCAM, TGFb, VEGF, PD-L1, MHCI, CD80, CD86, SELL
par(mfrow=c(1,1))
netVisual_aggregate(control_cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(control_cellchat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(control_cellchat, signaling = pathways.show, color.heatmap = "Reds")


# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
# names(group.cellType) <- levels(cellchat@idents)
# netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))


#Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(control_cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(control_cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(control_cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout="hierarchy")
par(mfrow=c(2,3))
netVisual_individual(control_cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
netVisual_individual(control_cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")


# Automatically save the plots of the all inferred network for quick exploration
# Access all the signaling pathways showing significant communications
pathways.show.all <- control_cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(control_cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(control_cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 12, height = 10, units = 'in', dpi = 300)
}





###################### Repeat for TAK981
TAK981_matrix <- GetAssayData(TAK981, assay = "RNA", slot = "data") # normalized data matrix
Idents(TAK981) <- "SingleR.label"
TAK981_labels <- Idents(TAK981)
TAK981_meta <- data.frame(group = TAK981_labels, row.names = names(TAK981_labels)) # create a dataframe of the cell labels
TAK981_cellchat <- createCellChat(object = TAK981_matrix, meta = TAK981_meta, group.by = "group")
TAK981_cellchat <- addMeta(TAK981_cellchat, meta = TAK981_meta, meta.name = "labels")
TAK981_cellchat <- setIdent(TAK981_cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(TAK981_cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(TAK981_cellchat@idents)) # number of cells in each cell group
TAK981_cellchat@DB <- CellChatDB
TAK981_cellchat <- subsetData(TAK981_cellchat) # This step is necessary even if using the whole database
TAK981_cellchat <- identifyOverExpressedGenes(TAK981_cellchat)
TAK981_cellchat <- identifyOverExpressedInteractions(TAK981_cellchat)
TAK981_cellchat <- projectData(TAK981_cellchat, PPI.mouse)
TAK981_cellchat <- computeCommunProb(TAK981_cellchat)
TAK981_cellchat <- filterCommunication(TAK981_cellchat, min.cells = 5)
TAK981_cellchat <- computeCommunProbPathway(TAK981_cellchat)
TAK981_cellchat <- aggregateNet(TAK981_cellchat)
groupSize <- as.numeric(table(TAK981_cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(TAK981_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(TAK981_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


#merge control and TAK981
object.list <- list(control = control_cellchat, TAK981 = TAK981_cellchat)

for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# cellchat <- filterCommunication(cellchat, min.cells = 40)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
par(mfrow = c(2,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macro-Spp1")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macro-C1q")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, -1, 1
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))

##########################################################################################################################

setwd("E:/TAK981_KPC/control_vs_TAK981")
#load('nonimmune.RData')
 
#InferCNV
kpc_matrix <- GetAssayData(scKPC, slot="counts", assay = "RNA")
annotations <- as.data.frame(scKPC$SingleR.label)
write.table(annotations, file="cellAnnotations.txt", sep="\t", col.names = FALSE,quote = FALSE)

kpc_cnv <- CreateInfercnvObject(raw_counts_matrix=kpc_matrix,
                                annotations_file="cellAnnotations.txt",
                                delim="\t",
                                gene_order_file="vM28.annotation.txt",
                                ref_group_names=c("Endothelial"))
memory.limit(56000)
kpc_cnv = infercnv::run(kpc_cnv,
                         cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                         out_dir="InferCNV",  # dir is auto-created for storing outputs
                         cluster_by_groups=T,   # cluster
                         denoise=T,
                         HMM=T)


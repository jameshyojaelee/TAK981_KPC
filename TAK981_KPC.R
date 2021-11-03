library(dplyr)
library(Seurat)
library(patchwork)

## import all necessary files
## use the following command to unzip tar.gz: tar -xvzf

C1.data <- Read10X(data.dir="control_1")
C2.data <- Read10X(data.dir="control_2")
T1.data <- Read10X(data.dir="TAK981_1")
T2.data <- Read10X(data.dir="TAK981_2")
P1.data <- Read10X(data.dir="TAK981_Anti-PDL1_1")
P2.data <- Read10X(data.dir="TAK981_Anti-PDL1_2")
aggr.data <-Read10X(data.dir = "aggr")


C1 <- CreateSeuratObject(counts=C1.data, project= "C1")
#C1 <- NormalizeData(C1)
C2 <- CreateSeuratObject(counts=C2.data, project= "C2")
#C2 <- NormalizeData(C2)
T1 <- CreateSeuratObject(counts=T1.data, project= "T1")
#T1 <- NormalizeData(T1)
T2 <- CreateSeuratObject(counts=T2.data, project= "T2")
#T2 <- NormalizeData(T2)
P1 <- CreateSeuratObject(counts=P1.data, project= "P1")
#P1 <- NormalizeData(P1)
P2 <- CreateSeuratObject(counts=P2.data, project= "P2")
#P2 <- NormalizeData(P2)

###############################################################################
#merge biological replicates
control <- merge(C1, y= C2, add.cell.ids = c("C1", "C2"), project="control")
control <- NormalizeData(control)
table(control$orig.ident)

TAK <- merge(T1, y=T2, add.cell.ids = c("T1", "T2"), project="TAK")
TAK <- NormalizeData(TAK)
table(TAK$orig.ident)

TP <- merge(P1, y=P2, add.cell.ids = c("P1", "P2"), project="TP")
TP <- NormalizeData(TP)
table(TP$orig.ident)


###############################################################################
#feature selection (highly variable features)
control <- FindVariableFeatures(control, nfeatures = 2000)
c.top10 <- head(VariableFeatures(control),10)
c.plot1 <- VariableFeaturePlot(control)
c.plot2 <- LabelPoints(plot=c.plot1, points=c.top10, repel=TRUE)
c.plot1 + c.plot2 

#feature selection (highly variable features)
TAK <- FindVariableFeatures(TAK)
t.top10 <- head(VariableFeatures(TAK),10)
t.plot1 <- VariableFeaturePlot(TAK)
t.plot2 <- LabelPoints(plot=t.plot1, points=t.top10, repel=TRUE)
t.plot1 + t.plot2 

#feature selection (highly variable features)
TP <- FindVariableFeatures(TP, nfeatures = 2000)
p.top10 <- head(VariableFeatures(TP),10)
p.plot1 <- VariableFeaturePlot(TP)
p.plot2 <- LabelPoints(plot=p.plot1, points=top10, repel=TRUE)
p.plot1 + p.plot2 

###############################################################################

#data scaling
all.genes <- rownames(control)
control <- ScaleData(control, features = all.genes)

all.genes <- rownames(TAK)
TAK <- ScaleData(TAK, features = all.genes)

all.genes <- rownames(TP)
TP <- ScaleData(TP, features = all.genes)

###############################################################################

#PCA

control <- RunPCA(control, feautures=VariableFeatures(object=control))
#VizDimLoadings(control, dims=1:2, reduction="pca")
#DimPlot(control, reduction="pca")

TAK <- RunPCA(TAK, feautures=VariableFeatures(object=TAK))
TP <- RunPCA(TP, feautures=VariableFeatures(object=TP))

###############################################################################
#Clustering
control <- FindNeighbors(control, dims=1:10)
control <- FindClusters(control, resolution = 0.5)

TAK <- FindNeighbors(TAK, dims=1:10)
TAK <- FindClusters(TAK, resolution = 0.5)

TP <- FindNeighbors(TP, dims=1:10)
TP <- FindClusters(TP, resolution = 0.5)

###############################################################################
#UMAP
control <- RunUMAP(control, dims=1:10)

TAK <- RunUMAP(TAK, dims=1:10)

TP <- RunUMAP(TP, dims=1:10)


DimPlot(control, reduction="umap")
DimPlot(TAK, reduction="umap")
DimPlot(TP, reduction="umap")


###############################################################################
#differential expressed features (cluster biomarkers)

#find markers for specific cluster
#control.cluster9 <-FindMarkers(control, ident.1 = 2, min.pct = 0.25)
#head(control.cluster9, n=5)

#find markers for all clusters
control.markers <- FindAllMarkers(control, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
control.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#find markers for all clusters
TAK.markers <- FindAllMarkers(TAK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TAK.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#find markers for all clusters
TP.markers <- FindAllMarkers(TP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TP.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


###############################################################################

#ViolinPlot
VlnPlot(control, features = c("Ptprc"))


#Feature Plot
FeaturePlot(control, features=c("Ptprc", "Cd3d", "Cd4", "Cd8b1", "Cd74", "Cd19"))


FeaturePlot(TAK, features=c("Ptprc", "Cd3d", "Cd4", "Cd8b1", "Cd74"))


FeaturePlot(TP, features=c("Ptprc", "Cd3d", "Cd4", "Cd8b1", "Cd74"))


###############################################################################
#Heatmap for top 20 markers for every cluster

control.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> c.top10
DoHeatmap(control, features = c.top10$gene) + NoLegend()

###############################################################################






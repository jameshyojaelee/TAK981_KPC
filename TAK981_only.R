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
setwd("X:/control_vs_TAK981")
source("scFunctions.R")
# 
# mref <- ImmGenData() #mouse immune cells
# sceM <- BaronPancreasData('mouse')
# sceM <- sceM[,!is.na(sceM$label)]
# library(scuttle)
# sceM <- logNormCounts(sceM)
# ## import all necessary files
# # load("immune.RData")
# load("nonimmune.RData")
load("new.RData")

## use the following command to unzip tar.gz: tar -xvzf
C1.data <- Read10X(data.dir="E:/TAK981_KPC/cellranger/c1_55M/filtered_feature_bc_matrix")
C2.data <- Read10X(data.dir="E:/TAK981_KPC/cellranger/c2_55M/filtered_feature_bc_matrix")
T1.data <- Read10X(data.dir="E:/TAK981_KPC/cellranger/TAK981_1/filtered_feature_bc_matrix")
T2.data <- Read10X(data.dir="E:/TAK981_KPC/cellranger/TAK981_2/filtered_feature_bc_matrix")

#############################################################################################
counts_per_cell <- Matrix::colSums(C1.data)
counts_per_gene <- Matrix::rowSums(C1.data)
genes_per_cell <- Matrix::colSums(C1.data>0) # count a gene only if it has non-zero reads mapped.
cells_per_gene <- Matrix::rowSums(C1.data>0) # only count cells where the gene is expressed
hist(log10(counts_per_cell+1),main='counts per cell',col='blue')
mito_genes <- grep("^mt-", rownames(C1.data) , ignore.case=T, value=T)
mito_gene_read_counts = Matrix::colSums(C1.data[mito_genes,])
pct_mito = mito_gene_read_counts / counts_per_cell * 100
plot(sort(pct_mito), xlab = "cells sorted by percentage mitochondrial counts", ylab = "percentage mitochondrial counts")

counts_per_cell <- Matrix::colSums(C2.data)
counts_per_gene <- Matrix::rowSums(C2.data)
genes_per_cell <- Matrix::colSums(C2.data>0) # count a gene only if it has non-zero reads mapped.
cells_per_gene <- Matrix::rowSums(C2.data>0) # only count cells where the gene is expressed
hist(log10(counts_per_cell+1),main='counts per cell',col='blue')
mito_genes <- grep("^mt-", rownames(C2.data) , ignore.case=T, value=T)
mito_gene_read_counts = Matrix::colSums(C2.data[mito_genes,])
pct_mito = mito_gene_read_counts / counts_per_cell * 100
plot(sort(pct_mito), xlab = "cells sorted by percentage mitochondrial counts", ylab = "percentage mitochondrial counts")

counts_per_cell <- Matrix::colSums(T1.data)
counts_per_gene <- Matrix::rowSums(T1.data)
genes_per_cell <- Matrix::colSums(T1.data>0) # count a gene only if it has non-zero reads mapped.
cells_per_gene <- Matrix::rowSums(T1.data>0) # only count cells where the gene is expressed
hist(log10(counts_per_cell+1),main='counts per cell',col='blue')
mito_genes <- grep("^mt-", rownames(T1.data) , ignore.case=T, value=T)
mito_gene_read_counts = Matrix::colSums(T1.data[mito_genes,])
pct_mito = mito_gene_read_counts / counts_per_cell * 100
plot(sort(pct_mito), xlab = "cells sorted by percentage mitochondrial counts", ylab = "percentage mitochondrial counts")

counts_per_cell <- Matrix::colSums(T2.data)
counts_per_gene <- Matrix::rowSums(T2.data)
genes_per_cell <- Matrix::colSums(T2.data>0) # count a gene only if it has non-zero reads mapped.
cells_per_gene <- Matrix::rowSums(T2.data>0) # only count cells where the gene is expressed
hist(log10(counts_per_cell+1),main='counts per cell',col='blue')
mito_genes <- grep("^mt-", rownames(T2.data) , ignore.case=T, value=T)
mito_gene_read_counts = Matrix::colSums(T2.data[mito_genes,])
pct_mito = mito_gene_read_counts / counts_per_cell * 100
plot(sort(pct_mito), xlab = "cells sorted by percentage mitochondrial counts", ylab = "percentage mitochondrial counts")

#############################################################################################

C1 <- CreateSeuratObject(counts=C1.data, project= "C1")
C1@meta.data$treatment <- "control"
C2 <- CreateSeuratObject(counts=C2.data, project= "C2")
C2@meta.data$treatment <- "control"
T1 <- CreateSeuratObject(counts=T1.data, project= "T1")
T1@meta.data$treatment <- "TAK981"
T2 <- CreateSeuratObject(counts=T2.data, project= "T2")
T2@meta.data$treatment <- "TAK981"

C1@assays$RNA #5556 cells
C2@assays$RNA #4596
T1@assays$RNA #2872
T2@assays$RNA #3866
#############################################################################################

#############################################################################################
#Sample by sample QC (https://www.nature.com/articles/s41388-021-02054-3#Sec22)
C1[["percent.mito"]] <- PercentageFeatureSet(C1, pattern="^mt-")
C2[["percent.mito"]] <- PercentageFeatureSet(C2, pattern="^mt-")
T1[["percent.mito"]] <- PercentageFeatureSet(T1, pattern="^mt-")
T2[["percent.mito"]] <- PercentageFeatureSet(T2, pattern="^mt-")

C1[["percent.ribo"]] <- PercentageFeatureSet(object = C1, pattern = "Rps|Rpl|Mrpl|Mrps")
C2[["percent.ribo"]] <- PercentageFeatureSet(object = C2, pattern = "Rps|Rpl|Mrpl|Mrps")
T1[["percent.ribo"]] <- PercentageFeatureSet(object = T1, pattern = "Rps|Rpl|Mrpl|Mrps")
T2[["percent.ribo"]] <- PercentageFeatureSet(object = T2, pattern = "Rps|Rpl|Mrpl|Mrps")

#nFeature_RNA: number of genes detected per cell
#nCount_RNA: number of UMIs per cell
VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA","percent.mito", "percent.ribo"),ncol = 4)
plot(C1$nCount_RNA,C1$nFeature_RNA, log='xy', col='blue')
plot(C1$nCount_RNA,C1$percent.mito, col='blue')

VlnPlot(C2, features = c("nFeature_RNA", "nCount_RNA","percent.mito", "percent.ribo"),ncol = 4)
plot(C2$nCount_RNA,C2$nFeature_RNA, log='xy', col='blue')
plot(C2$nCount_RNA,C2$percent.mito, col='blue') 

VlnPlot(T1, features = c("nFeature_RNA", "nCount_RNA","percent.mito", "percent.ribo"),ncol = 4)
plot(T1$nCount_RNA,T1$nFeature_RNA, log='xy', col='blue')
plot(T1$nCount_RNA,T1$percent.mito, col='blue') 

VlnPlot(T2, features = c("nFeature_RNA", "nCount_RNA","percent.mito", "percent.ribo"),ncol = 4)
plot(T2$nCount_RNA,T2$nFeature_RNA, log='xy', col='blue')
plot(T2$nCount_RNA,T2$percent.mito, col='blue')

#200-10,000 genes detected within cells, 1000-40,000 UMIs counted, less than 20% mitochondrial reads, less than 40% ribosombal genes
C1 <- subset(C1, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA<40000 & percent.mito < 10 & percent.ribo < 40)
C2 <- subset(C2, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA<40000 & percent.mito < 10 & percent.ribo < 40)
T1 <- subset(T1, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA<40000 & percent.mito < 10 & percent.ribo < 40)
T2 <- subset(T2, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA<40000 & percent.mito < 10 & percent.ribo < 40)

C1@assays$RNA
C2@assays$RNA
T1@assays$RNA
T2@assays$RNA

###################################################################################################################
#USING SCTRANSFORM INSTEAD. Better than the traditional approach above.
memory.limit(56000)
KPC <- merge(C1, y=c(C2, T1, T2), add.cell.ids = c("C1","C2","T1","T2"), project="KPC")
KPC.list <- SplitObject(KPC, split.by = "orig.ident")
for (i in 1:length(KPC.list)) {
  KPC.list[[i]] <- SCTransform(KPC.list[[i]], vars.to.regress=c("percent.mito", "percent.ribo"))
}
features <- SelectIntegrationFeatures(object.list = KPC.list, nfeatures = 2000)
KPC.list <- PrepSCTIntegration(object.list = KPC.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = KPC.list, normalization.method = "SCT",
                                         anchor.features = features)
combined.KPC <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
combined.KPC <- RunPCA(combined.KPC, verbose = FALSE)
combined.KPC <- RunTSNE(combined.KPC, dims = 1:40, verbose = FALSE)
combined.KPC <- RunUMAP(combined.KPC, dims = 1:40, verbose = FALSE)
combined.KPC <- FindNeighbors(combined.KPC, reduction = "umap",dims = 1:2, verbose = FALSE)
combined.KPC <- FindClusters(combined.KPC, resolution = 0.2)
DimPlot(combined.KPC, reduction="umap",label=TRUE,pt.size=1)

# find markers for all clusters
# MUST USE "RNA" for assay when using samples with different conditions (i.e. control vs treated)
DefaultAssay(combined.KPC) <- "RNA"
Idents(combined.KPC) <- "SingleR.label"
#Must normalize data with Lognormalize for further DE analyses
combined.KPC <- NormalizeData(object =combined.KPC, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(combined.KPC)
combined.KPC <- ScaleData(object = combined.KPC, features = all.genes)
KPC.markers <- FindAllMarkers(combined.KPC, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
KPC.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10
write.csv(KPC.markers, "kpc_markers.csv", row.names = T)
#Heatmap for top 10 markers for every cluster

tiff("cluster_heatmap.jpeg", unit="in", width=20, height=16, res=500)
DoHeatmap(combined.KPC, features = top10$gene) + NoLegend()
dev.off()
# 
# #Conserved markers between two treatment groups
# KPC.conserved.markers <- FindConservedMarkers(combined.KPC, assay="RNA",only.pos = TRUE, idents.1=0,grouping.var = "treatment")
# KPC.conserved.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC) -> conserved.top10

###############################################################################
tiff("tSNE_unlabeled2.jpeg", unit="in", width=8, height=7, res=500)
DimPlot(combined.KPC, reduction="umap",label=TRUE,pt.size=1, split.by = "treatment")
dev.off()

VlnPlot(combined.KPC, "Ptprc", pt.size=0.5)+NoLegend()
# immune <- subset(combined.KPC, idents=c(0:23)[-c(1,4,5,8,9,10,12,16,17,19,20,21)])
# # cluster 5(or #6) is mt-genes. exclding...
# nonimmune <- subset(combined.KPC, idents=c(0:23)[c(1,4,5,8,9,10,12,16,17,19,20,21)])

sccombined.KPC <- as.SingleCellExperiment(combined.KPC)
sccombined.KPC <- SingleR(test=sccombined.KPC, ref=mref, assay.type.test = 1, labels = mref$label.main)
combined.KPC[["SingleR.label"]] <- sccombined.KPC$labels
#Feature Plot (reference to https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00776-9#ethics)

# temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==20, "cDC2-Cd209a") #Cd11b+ Cd11c(int) Cdh1+ Clec9a- Xcr1- CD209a+
# temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==25, "cDC1-Ccl22") #Relb+ Itgae(Cd103)- CD8a-  Xcr1- Cd209a- Clec9a- Ly75(Cd205)+ 
# temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==23, "cDC2-Itgax") #expressing Itgax and Mgl2, Ear2, 
# temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==26, "cDC1-Clec9a") #Itgae(Cd103)+ CD8a-  Xcr1+  Clec9a+ Batf3+
# temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==29, "pDC") #Cd11b- Siglech+ Bst2+ Tlr7+ Clec9a+ Ly6c2+ 


FeaturePlot(combined.KPC, features=c("Col1a1", "Col1a2","Col4a1", "Col6a1", "Sdc4"), split.by = "treatment")
FeaturePlot(combined.KPC, features=c("Cd47"), split.by = "treatment")
VlnPlot(combined.KPC, features=c("Sirpa"), group.by="SingleR.label", split.by = "treatment")


combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==0, "Ductal")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==1, "Ductal")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==2, "Myeloid")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==3, "Myeloid")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==4, "Ductal")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==5, "Myeloid")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==6, "Myeloid")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==7, "CAF")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==8, "Myeloid")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==9, "Ductal")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==10, "T/NK cells")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==11, "DC")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==12, "B cells")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==13, "Ductal")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==14, "Granulocytes") 
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==15, "Endothelial")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==16, "Ductal")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==17, "DC")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==18, "Ductal")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==19, "Granulocytes")

tiff("umap_split.jpeg", unit="in", width=10, height=5, res=500)
DimPlot(combined.KPC, reduction="umap",label=F,pt.size=0.4, group.by = "SingleR.label")
dev.off()

tiff("KPC_CD40.jpeg", unit="in", width=8, height=6, res=500)
FeaturePlot(combined.KPC, feature=c("Cd40"),pt.size=0.7)
dev.off()

DimPlot(combined.KPC, reduction="umap",label=TRUE,pt.size=1)

VlnPlot(combined.KPC, features=c("Vegfa"),group.by = "SingleR.label", split.by = "treatment")

tiff("umap.jpeg", unit="in", width=7, height=5, res=500)
DimPlot(combined.KPC, reduction="umap",label=F,pt.size=0.3, group.by = "SingleR.label")
dev.off()

tiff("Sumo1_total_vlnplot.jpeg", unit="in", width=7, height=5, res=500)
VlnPlot(combined.KPC, features=c('Sumo1'),group.by = "SingleR.label", split.by="treatment")
dev.off()
tiff("Sumo2_total_vlnplot.jpeg", unit="in", width=7, height=5, res=500)
VlnPlot(combined.KPC, features=c('Sumo2'),group.by = "SingleR.label", split.by="treatment")
dev.off()
tiff("Sumo3_total_vlnplot.jpeg", unit="in", width=7, height=5, res=500)
VlnPlot(combined.KPC, features=c('Sumo3'),group.by = "SingleR.label", split.by="treatment")
dev.off()
tiff("Ube2i_total_vlnplot.jpeg", unit="in", width=7, height=5, res=500)
VlnPlot(combined.KPC, features=c('Ube2i'),group.by = "SingleR.label", split.by="treatment")
dev.off()
tiff("Uba2_total_vlnplot.jpeg", unit="in", width=7, height=5, res=500)
VlnPlot(combined.KPC, features=c('Uba2'),group.by = "SingleR.label", split.by="treatment")
dev.off()

FeaturePlot(combined.KPC, features=c("H2-K1"))
FeaturePlot(combined.KPC, features=c("H2-D1"))
VlnPlot(combined.KPC, features=c('H2-K1'),group.by = "SingleR.label", split.by="treatment")
VlnPlot(combined.KPC, features=c('H2-D1'),group.by = "SingleR.label", split.by="treatment")
VlnPlot(combined.KPC, features=c('Gapdh'),group.by = "SingleR.label", split.by="treatment")
VlnPlot(combined.KPC, features=c('Actb'),group.by = "SingleR.label", split.by="treatment")


tiff("SUMO_total.jpeg", unit="in", width=12, height=12, res=500)
FeaturePlot(combined.KPC, features=c("Sumo1", "Sumo2","Sumo3", "Sae1","Uba2", "Ube2i"))
dev.off()

tiff("SUMO_total_vlnplot.jpeg", unit="in", width=17, height=8, res=500)
VlnPlot(combined.KPC, features=c("Sumo1", "Sumo2","Sumo3", "Sae1","Uba2", "Ube2i"),group.by = "SingleR.label", split.by = "treatment",pt.size=0)
dev.off()

dot_features <- c("Cd3e","Gzmb","Nkg7","Apoe", "C1qc","Lyz2", "S100a8", "S100a9", "G0s2", "Pecam1", "Cd34",
                  "Krt18", "Krt19", "Itgax", "Ccl22","Dcn", "Col1a1", "Col3a1", "Ms4a1", "Cd79a")

tiff("cellmarkers_dotplot.jpeg", unit="in", width=10, height=4, res=500)
DotPlot(combined.KPC, features = dot_features, group.by="SingleR.label") +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(1,10)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

tiff("combined.KPC_celltype_proportion_treatment.jpeg", unit="in", width=9, height=6, res=500)
plot_group_proportions(combined.KPC, graph.type = "dodge")
dev.off()

tiff("combined.KPC_celltype_proportion_treatment_stacked.jpeg", unit="in", width=3, height=6, res=500)
plot_group_proportions(combined.KPC, graph.type = "stacked")
dev.off()

combined.KPC_prop <- sc_utils(combined.KPC)

prop_test <- permutation_test(
  combined.KPC_prop, cluster_identity = "SingleR.label",
  sample_1 = "control", sample_2 = "TAK981",
  sample_identity = "treatment"
)
permutation_plot(prop_test)+ scale_colour_discrete(labels = function(x) str_wrap(x, width = 7))

tiff("kpc_Ifi27.jpeg", unit="in", width=5, height=5, res=500)
FeaturePlot(combined.KPC, "Ifi27")
dev.off()

tiff("combined.KPC_celltype_proportion_test.jpeg", unit="in", width=5, height=1.8, res=500)
permutation_plot(prop_test)+ scale_colour_discrete(labels = function(x) str_wrap(x, width = 7))
dev.off()

VlnPlot(combined.KPC, features=c("Cd74"),group.by = "SingleR.label", split.by = "treatment")
FeaturePlot(combined.KPC, features=c("Cd74"))


dot_features <- c("Cd40", "Cd47", "H2-Aa", "H2-D1")
VlnPlot(combined.KPC, features = c("Cd40"), group.by = "SingleR.label", split.by = "treatment", pt.size=0.01)
dittoDotPlot(combined.KPC, vars = dot_features, group.by= "treatment", split.by = "SingleR.label", split.nrow = 8,
             size=10, max.color="red", min.color = "yellow") +
  theme(axis.title.y=element_blank(), axis.text.x = element_text(angle=30))


immune <- subset(combined.KPC, idents=c(0:19)[-c(1,2,5,8,10,14,16,17,19)])
nonimmune <- subset(combined.KPC, idents=c(0:19)[c(1,2,5,8,10,14,16,17,19)])
VlnPlot(immune, features=c("Krt18"))

#save data
save(immune, file="immune.RData")
save(nonimmune, file="nonimmune.RData")

combined.KPC@assays$RNA@data
combined.KPC@assays$SCT@data


df <- as.data.frame(combined.KPC@assays$RNA@data)
write.csv(df, file )

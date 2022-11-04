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
library(DoubletFinder)
library("scProportionTest")
setwd("E:/TAK981_KPC/control_vs_TAK981")

mref <- ImmGenData() #mouse immune cells
sceM <- BaronPancreasData('mouse')
sceM <- sceM[,!is.na(sceM$label)]
library(scuttle)
sceM <- logNormCounts(sceM)
## import all necessary files
# load("immune.RData")
# load("nonimmune.RData")
load("TAK981_only.RData")

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
C1 <- subset(C1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA<40000 & percent.mito < 20 & percent.ribo < 40)
C2 <- subset(C2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA<40000 & percent.mito < 20 & percent.ribo < 40)
T1 <- subset(T1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA<40000 & percent.mito < 20 & percent.ribo < 40)
T2 <- subset(T2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 1000 & nCount_RNA<40000 & percent.mito < 20 & percent.ribo < 40)

C1@assays$RNA
C2@assays$RNA
T1@assays$RNA
T2@assays$RNA

###################################################################################################################
#USING SCTRANSFORM INSTEAD. Better than the traditional approach above.

memory.limit(56000)
C1 <- NormalizeData(C1)
C1 <- FindVariableFeatures(C1, selection.method = "vst", nfeatures = 2000)
C1 <- ScaleData(C1)
C1 <- RunPCA(C1)
C1 <- RunUMAP(C1, dims = 1:10)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(C1, PCs = 1:10, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# # ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
# sweep.res.list <- paramSweep_v3(control, PCs = 1:10, sct = F)
# gt.calls <- control@meta.data[rownames(sweep.res.list[[1]]), "GT"].   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results
# sweep.stats <- summarizeSweep(sweep.res.list, GT = TRUE, GT.calls = gt.calls)
# bcmvn <- find.pK(sweep.stats)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- control@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- control@meta.data$ClusteringResults
nExp_poi <- round(0.035*nrow(C1@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
# nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
C1 <- doubletFinder_v3(C1, PCs = 1:10, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
DF.name = colnames(C1@meta.data)[grepl("DF.classification", colnames(C1@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(C1, group.by = "orig.ident") + NoAxes(),
                   DimPlot(C1, group.by = DF.name) + NoAxes())
C1 = C1[, C1@meta.data[, DF.name] == "Singlet"]


C2 <- NormalizeData(C2)
C2 <- FindVariableFeatures(C2, selection.method = "vst", nfeatures = 2000)
C2 <- ScaleData(C2)
C2 <- RunPCA(C2)
C2 <- RunUMAP(C2, dims = 1:10)
sweep.res.list <- paramSweep_v3(C2, PCs = 1:10, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
nExp_poi <- round(0.035*nrow(C2@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
C2 <- doubletFinder_v3(C2, PCs = 1:10, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
DF.name = colnames(C2@meta.data)[grepl("DF.classification", colnames(C2@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(C2, group.by = "orig.ident") + NoAxes(),
                   DimPlot(C2, group.by = DF.name) + NoAxes())
C2 = C2[, C2@meta.data[, DF.name] == "Singlet"]

T1 <- NormalizeData(T1)
T1 <- FindVariableFeatures(T1, selection.method = "vst", nfeatures = 2000)
T1 <- ScaleData(T1)
T1 <- RunPCA(T1)
T1 <- RunUMAP(T1, dims = 1:10)
sweep.res.list <- paramSweep_v3(T1, PCs = 1:10, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
nExp_poi <- round(0.023*nrow(T1@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
T1 <- doubletFinder_v3(T1, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
DF.name = colnames(T1@meta.data)[grepl("DF.classification", colnames(T1@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(T1, group.by = "orig.ident") + NoAxes(),
                   DimPlot(T1, group.by = DF.name) + NoAxes())
T1 = T1[, T1@meta.data[, DF.name] == "Singlet"]


T2 <- NormalizeData(T2)
T2 <- FindVariableFeatures(T2, selection.method = "vst", nfeatures = 2000)
T2 <- ScaleData(T2)
T2 <- RunPCA(T2)
T2 <- RunUMAP(T2, dims = 1:10)
sweep.res.list <- paramSweep_v3(T2, PCs = 1:10, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
nExp_poi <- round(0.031*nrow(T2@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
T2 <- doubletFinder_v3(T2, PCs = 1:10, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
DF.name = colnames(T2@meta.data)[grepl("DF.classification", colnames(T2@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(T2, group.by = "orig.ident") + NoAxes(),
                   DimPlot(T2, group.by = DF.name) + NoAxes())
T2 = T2[, T2@meta.data[, DF.name] == "Singlet"]


###################################################################################################################
KPC <- merge(C1, y=c(C2,T1,T2), add.cell.ids = c("C1","C2","T1","T2"), project="KPC")
KPC.list <- SplitObject(KPC, split.by = "orig.ident")
memory.limit(56000)
for (i in 1:length(KPC.list)) {KPC.list[[i]] <- SCTransform(KPC.list[[i]], vars.to.regress=c("nCount_RNA","percent.mito", "percent.ribo"))}
KPC.features <- SelectIntegrationFeatures(object.list = KPC.list, nfeatures = 2000)
KPC.list <- PrepSCTIntegration(object.list = KPC.list, anchor.features = KPC.features)
KPC.anchors <- FindIntegrationAnchors(object.list = KPC.list, normalization.method = "SCT",
                                          anchor.features = KPC.features)
combined.KPC <- IntegrateData(anchorset = KPC.anchors, normalization.method = "SCT")
combined.KPC <- ScaleData(object = combined.KPC, verbose = FALSE)
combined.KPC <- RunPCA(combined.KPC, verbose = FALSE)
combined.KPC <- RunTSNE(combined.KPC, dims = 1:40, verbose = FALSE)
combined.KPC <- RunUMAP(combined.KPC, dims = 1:40, verbose = FALSE)
combined.KPC <- FindNeighbors(combined.KPC, reduction = "umap",dims = 1:2, verbose = FALSE)
combined.KPC <- FindClusters(combined.KPC, resolution = 0.1)
DimPlot(combined.KPC, reduction="umap",label=TRUE,pt.size=1)

# find markers for all clusters
# MUST USE "RNA" for assay when using samples with different conditions (i.e. KPC vs treated)
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
tiff("cluster_heatmap_temp.jpeg", unit="in", width=40, height=24, res=500)
DoHeatmap(combined.KPC, features = top10$gene) + NoLegend()
dev.off()

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


# 
# 
# temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==20, "cDC2-Cd209a") #Cd11b+ Cd11c(int) Cdh1+ Clec9a- Xcr1- CD209a+
# temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==25, "cDC1-Ccl22") #Relb+ Itgae(Cd103)- CD8a-  Xcr1- Cd209a- Clec9a- Ly75(Cd205)+ 
# temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==23, "cDC2-Itgax") #expressing Itgax and Mgl2, Ear2, 
# temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==26, "cDC1-Clec9a") #Itgae(Cd103)+ CD8a-  Xcr1+  Clec9a+ Batf3+
# temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==29, "pDC") #Cd11b- Siglech+ Bst2+ Tlr7+ Clec9a+ Ly6c2+ 
DimPlot(combined.KPC, reduction="umap",label=TRUE,pt.size=1)
VlnPlot(combined.KPC, features=c("Cd3e"))
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==0, "Ductal")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==1, "Ductal")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==2, "Myeloid")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==3, "Myeloid")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==4, "Myeloid")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==5, "CAF")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==6, "Granulocytes")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==7, "T/NK cells")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==8, "B cells")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==9, "Myeloid")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==10, "Endothelial")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==11, "ADM")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==12, "Myeloid")
combined.KPC[["SingleR.label"]]<- replace(combined.KPC[["SingleR.label"]], combined.KPC[["seurat_clusters"]]==13, "Ductal")
DimPlot(combined.KPC, reduction="umap",label=TRUE,pt.size=1, group.by = "SingleR.label")

tiff("tSNE.jpeg", unit="in", width=8, height=4, res=500)
DimPlot(combined.KPC, reduction="tsne",label=TRUE,pt.size=1, group.by = "SingleR.label")
dev.off()

dot_features <- c("Cd3e","Gzmb","Nkg7","Siglech","Ly6c2","Apoe", "C1qc","Lyz2", "S100a8", "S100a9", "G0s2", "Dcn", "Col1a1", "Col3a1", "Pecam1", "Cd34",
                  "Krt18", "Krt19", "Itgax", "Ccl22","Ms4a1", "Cd79a","Ctrb1", "Prss2", "Try5") 

tiff("cellmarkers_dotplot.jpeg", unit="in", width=12, height=5, res=500)
DotPlot(combined.KPC, features = dot_features, group.by="SingleR.label") +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(1,10)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

immune <- subset(combined.KPC, idents=c(0:13)[-c(1,2,6,11,12,14)])
nonimmune <- subset(combined.KPC, idents=c(0:23)[c(1,2,6,11,12,14)])
VlnPlot(nonimmune, features=c("Ptprc"))

save(immune, file="immune.RData")
save(nonimmune, file="nonimmune.RData")



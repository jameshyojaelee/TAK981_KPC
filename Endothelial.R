#nonendo cell analysis
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
setwd("D:/TAK981_KPC/control_vs_TAK981")
library("scProportionTest")
source("scFunctions.R")

load("endo.RData")

DefaultAssay(endo) <- "RNA"
endo <- NormalizeData(endo)
all.genes <- rownames(endo)
endo <- ScaleData(endo, features = all.genes,verbose = FALSE)
endo <- FindVariableFeatures(endo, selection.method = "vst", nfeatures = 2500,verbose = FALSE)
endo <- RunPCA(endo, verbose = FALSE)
endo <- RunUMAP(endo, dims = 1:40,verbose = FALSE)
endo <- FindNeighbors(endo, dims = 1:2, reduction = "umap", verbose = FALSE)
endo <- FindClusters(endo, resolution = 1)
DimPlot(endo, reduction="umap",label=TRUE,pt.size=1, group.by = "SingleR.label", repel=T)

dot_features <- c("Vegfa","Des","Kdr","Hif1a")
DotPlot(endo, features = dot_features, group.by = "SingleR.label", split.by = "treatment")
















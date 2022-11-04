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
library(reticulate)
library("scProportionTest")
setwd("D:/TAK981_KPC/control_vs_TAK981")
source("scFunctions.R")
load("TAN.RData")

TAN <- DietSeurat(TAN, assay="RNA")
DefaultAssay(TAN) <- "RNA"
TAN <- NormalizeData(TAN,verbose = FALSE)
TAN <- FindVariableFeatures(TAN, selection.method = "vst", nfeatures = 2000)
features <- rownames(TAN)
TAN <- ScaleData(object = TAN, features = features,verbose = FALSE)
TAN <- RunPCA(TAN, verbose = FALSE)
TAN <- RunTSNE(TAN, dims = 1:40, verbose = FALSE)
TAN <- RunUMAP(TAN, dims = 1:40, verbose = FALSE)
TAN <- FindNeighbors(TAN, dims = 1:40, verbose = FALSE)
TAN <- FindClusters(TAN, resolution = 1,verbose = FALSE)
DimPlot(TAN, reduction="umap",label=TRUE,pt.size=2, group.by = "SingleR.label" )

mref <- ImmGenData()
SRTNK <- as.SingleCellExperiment(TAN)
SRTNK <- SingleR(test=SRTNK, ref=mref, assay.type.test = 1, labels = mref$label.main)
TAN[["SingleR.label"]] <- SRTNK$labels

TAN_subset <- subset(TAN, subset=(SingleR.label=="Neutrophils"))

VlnPlot(TAN, "Spp1")

dot_features <- c("Tnf", "Ccl3", "Icam1", "Fas",
                  "Cxcr4","Mmp9", "Il1b","Hmgb1",
                  "Osm", "Vegfa", 
                  "Arg1", "Cd274", "Vsir","Tgfbr1", "Tgfbr2")

tiff("N1 N2 TAN.jpeg", unit="in", width=9, height=3, res=500)
DotPlot(TAN_subset, features = dot_features, group.by="treatment") +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(0,15)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

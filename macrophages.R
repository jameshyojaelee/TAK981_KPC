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
setwd("X:/control_vs_TAK981")
source("scFunctions.R")
load("Macrophages.RData")

Macrophages <- DietSeurat(Macrophages, assay="RNA")
Macrophages.list <- SplitObject(Macrophages, split.by = "orig.ident")

for (i in 1:length(Macrophages.list)) {
  Macrophages.list[[i]] <- SCTransform(Macrophages.list[[i]], vars.to.regress=c("percent.mito", "percent.ribo"))
}
features <- SelectIntegrationFeatures(object.list = Macrophages.list, nfeatures = 2000)
Macrophages.list <- PrepSCTIntegration(object.list = Macrophages.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = Macrophages.list, normalization.method = "SCT",
                                  anchor.features = features)
Macrophages <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
DefaultAssay(Macrophages) <- "integrated"
Macrophages <- ScaleData(object = Macrophages, features=features)
Macrophages <- RunPCA(Macrophages, verbose = FALSE)
# Macrophages <- RunTSNE(Macrophages, dims = 1:40, verbose = FALSE)
Macrophages <- RunUMAP(Macrophages, dims = 1:40,verbose = FALSE)
Macrophages <- FindNeighbors(Macrophages,dims = 1:40, verbose = FALSE)
Macrophages <- FindClusters(Macrophages, resolution = 0.8)
DimPlot(Macrophages, reduction="umap",label=TRUE,pt.size=1)
DimPlot(Macrophages, reduction="umap",label=TRUE,pt.size=1, group.by = "SingleR.label")

# lognormalize for marker genes
DefaultAssay(Macrophages) <- "RNA"
Macrophages <- NormalizeData(Macrophages) #LogNormalize
all.genes <- rownames(Macrophages)
Macrophages <- ScaleData(object = Macrophages, features = all.genes)
Macrophages.markers <- FindAllMarkers(Macrophages, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Macrophages.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> Macrophages.top10
write.csv(Macrophages.markers, "Macrophages.markers.csv", row.names = T)
tiff("Macrophages.markers_heatmap.jpeg", unit="in", width=14, height=10, res=500)
DoHeatmap(Macrophages, features = Macrophages.top10) + NoLegend()
dev.off()

FeaturePlot(Macrophages, "Col1a2")
VlnPlot(Macrophages, "Cd40")
VlnPlot(Macrophages, "Cd40", group.by = "SingleR.label", split.by = "treatment")

Macrophages[["SingleR.label"]] <- replace(Macrophages[["SingleR.label"]], Macrophages[["seurat_clusters"]]==0, "Macro-MHCII") #MHC Class II high macropahge (cd74 is invariant MHCii) https://www.frontiersin.org/articles/10.3389/fimmu.2018.01132/full
Macrophages[["SingleR.label"]] <- replace(Macrophages[["SingleR.label"]], Macrophages[["seurat_clusters"]]==1, "Macro-Spp1")
Macrophages[["SingleR.label"]] <- replace(Macrophages[["SingleR.label"]], Macrophages[["seurat_clusters"]]==2, "Macro-C1q") 
Macrophages[["SingleR.label"]] <- replace(Macrophages[["SingleR.label"]], Macrophages[["seurat_clusters"]]==3, "Macro-Proliferating") #Ki67 high
Macrophages[["SingleR.label"]] <- replace(Macrophages[["SingleR.label"]], Macrophages[["seurat_clusters"]]==4, "Mono/MDSC") #high Il1b, Ly6c2, Ccr2. pro-tumorogensis genes such as Lyz2, Ifitm3, Vim, S100a6, https://www.jimmunol.org/content/206/1_Supplement/101.01
Macrophages[["SingleR.label"]] <- replace(Macrophages[["SingleR.label"]], Macrophages[["seurat_clusters"]]==5, "Macro-Spp1") 
Macrophages[["SingleR.label"]] <- replace(Macrophages[["SingleR.label"]], Macrophages[["seurat_clusters"]]==6, "Macro-C1q")  #FCN1 (Fcna+)
Macrophages[["SingleR.label"]] <- replace(Macrophages[["SingleR.label"]], Macrophages[["seurat_clusters"]]==7, "Mono-Ly6c") #Ly6C+CCR2+ indicates BMDM. highly expressed genes are related to IFN response (IFIT, ISG, OAS).  CXCL10 high
Macrophages[["SingleR.label"]] <- replace(Macrophages[["SingleR.label"]], Macrophages[["seurat_clusters"]]==8, "Macro-Spp1") #MMP9+ CAF-related 

DimPlot(Macrophages, reduction="umap",label=TRUE,pt.size=1, group.by = "SingleR.label")
VlnPlot(Macrophages, features=c("Cd47"), group.by = "SingleR.label", split.by = "treatment")
VlnPlot(Macrophages, features=c("Sirpa"), group.by = "SingleR.label", split.by = "treatment")
VlnPlot(Macrophages, features=c("Sirpa"), group.by = "SingleR.label", split.by = "treatment")

########################################################################################
library(EnhancedVolcano)

DefaultAssay(Macrophages) <- "RNA"

Macrophages.temp <- Macrophages
Idents(Macrophages.temp) <- "treatment" #set all idents to orig.ident i.e. control, TAK981, etc
Macrophages_ranked_list <- FindMarkers(Macrophages.temp, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)
EnhancedVolcano(Macrophages_ranked_list,
                lab = rownames(Macrophages_ranked_list),
                x = 'avg_log2FC',
                y = 'p_val')

Macrophages.temp <- subset(Macrophages)
Idents(Macrophages.temp) <- "treatment" #set all idents to orig.ident i.e. control, TAK981, etc

Macrophages_ranked_list <- FindMarkers(Macrophages.temp, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0, logfc.threshold = 0)
EnhancedVolcano(Macrophages_ranked_list,
                lab = rownames(Macrophages_ranked_list),
                x = 'avg_log2FC',
                y = 'p_val')

tiff("control_vs_TAK981_TAM.jpeg", unit="in", width=15, height=10, res=500)
EnhancedVolcano(Macrophages_ranked_list, lab = rownames(Macrophages_ranked_list),
                x = "avg_log2FC", y = "p_val",
                pCutoff = 0.01, FCcutoff = 0, labSize = 5, pointSize = 2,
                selectLab = c("Cd40","Cd80", "Cd86", "Sod2", "Nos2", "Stat1", "Il1b","Il18","Cxcl9", "Cxcl10", "Cxcl16",
                              "Cd68","Mrc1", "Msr1", "Chil3", "Retnla", "Arg1", "Stat3", 
                              "Il10ra", "Il6",
                              "Cd163", "Tgfb1",
                              "Vegfa",
                              "Ccl5", "Ccl7", "Ly6c2","Isg15"),
                labCol = 'black',
                labFace = 'bold',
                parseLabels = TRUE,
                legendPosition = "none", 
                title = "Control vs TAK981 TAM",
                drawConnectors = TRUE,
                widthConnectors = 0.6,
                arrowheads = FALSE,
                colConnectors = 'black',
                max.overlaps = 30)
dev.off()


########################################################################################

# Macrophages
save(Macrophages,file="Macrophages.RData")

#cIAP12 paper and drug response Cell paper

FeaturePlot(Macrophages,features=c("Cxcl9", "Cxcl10"), pt.size = 1.2)
FeaturePlot(Macrophages,features=c("Ccl5"), pt.size = 1.2)

VlnPlot(Macrophages, features = c("Cxcl9", "Cxcl10"), group.by = "SingleR.label")
VlnPlot(Macrophages, features = c("Ccl5"), group.by = "SingleR.label", split.by = "treatment")

VlnPlot(Macrophages, features = c("Spp1"), group.by = "SingleR.label")
VlnPlot(Macrophages, features = c("H2-Aa"))
VlnPlot(Macrophages, features = c("H2-Eb1"))
VlnPlot(Macrophages, features = c("Cd74"))
VlnPlot(Macrophages, features = c("Il1b"))
VlnPlot(Macrophages, features = c("C1qc"))
VlnPlot(Macrophages, features = c("Maf"))
VlnPlot(Macrophages, features = c("Mafb"))

dot_features <- c("Cd40","Cd80", "Cd86", "Sod2", "Nos2", "Stat1", "Il1b","Irf3","Irf5","Tlr4","Il18","Cxcl9", "Cxcl10", "Cxcl16",
                  "Cd68","Mrc1", "Chil3", "Retnla", "Arg1", "Stat3",
                  "Il10ra", "Il6",
                  "Cd163", "Tgfb1","Mertk",
                  "Vegfa")
Idents(Macrophages) <- "SingleR.label"
DotPlot(Macrophages, features = dot_features, split.by="treatment", cols ="YlOrRd") +
  scale_size(range = c(1,15)) + RotatedAxis() + theme(axis.title = element_blank())



FeaturePlot(Macrophages,features=c("Cxcl9", "Cxcl10"), pt.size = 1.2)
VlnPlot(Macrophages,features=c("Cxcl10"), group.by="SingleR.label", split.by = "treatment")

tiff("Sumo1_macrophage_vlnplot.jpeg", unit="in", width=7, height=5, res=500)
VlnPlot(Macrophages, features=c('Sumo1'),group.by = "SingleR.label", split.by="treatment")
dev.off()
tiff("Sumo2_macrophage_vlnplot.jpeg", unit="in", width=7, height=5, res=500)
VlnPlot(Macrophages, features=c('Sumo2'),group.by = "SingleR.label", split.by="treatment")
dev.off()
tiff("Sumo3_macrophage_vlnplot.jpeg", unit="in", width=7, height=5, res=500)
VlnPlot(Macrophages, features=c('Sumo3'),group.by = "SingleR.label", split.by="treatment")
dev.off()
tiff("Ube2i_macrophage_vlnplot.jpeg", unit="in", width=7, height=5, res=500)
VlnPlot(Macrophages, features=c('Ube2i'),group.by = "SingleR.label", split.by="treatment")
dev.off()
tiff("Uba2_macrophage_vlnplot.jpeg", unit="in", width=7, height=5, res=500)
VlnPlot(Macrophages, features=c('Uba2'),group.by = "SingleR.label", split.by="treatment")
dev.off()


tiff("Itgam_macrophage_vlnplot.jpeg", unit="in", width=7, height=5, res=500)
VlnPlot(Macrophages, features=c('Itgam'),group.by = "SingleR.label", split.by="treatment")
dev.off()


source("DotPlot_custom.R")
DotPlot(Macrophages, features = dot_features, group.by="orig.ident") +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(1,10)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

dittoDotPlot(Macrophages2, vars = dot_features, group.by = "treatment", split.by="SingleR.label", 
             adjustment = "z-score", legend.show = F,
             split.nrow = 6,
             size=10, max.color="red", min.color = "yellow")+
  theme(axis.title.y=element_blank(),
        axis.text.x = element_text(angle=30))

dittoHeatmap(Macrophages, genes = dot_features, order.by = "SingleR.label", annot.by=c("SingleR.label","treatment"),
             scaled.to.max = TRUE)


dot_features <- c("Cd40","Cd80", "Cd86", "Sod2", "Nos2", "Stat1", "Il1b","Irf3","Irf5","Tlr4","Il18","Cxcl9", "Cxcl10",
                  "Mrc1", "Msr1", "Chil3", "Arg1", "Stat3",
                  "Il10ra", "Il6",
                  "Cd163", "Tgfb1","Mertk",
                  "Vegfa")

tiff("M1 M2 Macrophages2.jpeg", unit="in", width=12, height=2.5, res=500)
DotPlot(Macrophages, features = dot_features, group.by="treatment") +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(0,10)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
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

# mref <- ImmGenData() #mouse Macrophages cells


Isg15 <- VlnPlot(Macrophages, features = c("Isg15"), group.by = "SingleR.label") +
  ylab("Isg15") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 40, margin=margin(r = 30), angle=0,vjust = 0.55),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Isg15
Mki67 <- VlnPlot(Macrophages, features = c("Mki67"), group.by = "SingleR.label") +
  ylab("Mki67") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 40, margin=margin(r = 30), angle=0,vjust = 0.55),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

C1qa <- VlnPlot(Macrophages, features = c("C1qa"), group.by = "SingleR.label") +
  ylab("C1qa") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 40, margin=margin(r = 30), angle=0,vjust = 0.55),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

Spp1 <- VlnPlot(Macrophages, features = c("Spp1"), group.by = "SingleR.label") +
  ylab("Spp1") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 40, margin=margin(r = 30), angle=0,vjust = 0.55),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

Ly6c2 <- VlnPlot(Macrophages, features = c("Ly6c2"), group.by = "SingleR.label") +
  ylab("Ly6c2") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 40, margin=margin(r = 30), angle=0,vjust = 0.55),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

Chil3 <- VlnPlot(Macrophages, features = c("Chil3"), group.by = "SingleR.label") +
  ylab("Chil3") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 40, margin=margin(r = 30), angle=0,vjust = 0.55),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Ccr2 <- VlnPlot(Macrophages, features = c("Ccr2"), group.by = "SingleR.label") +
  ylab("Ccr2") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.title.y.left = element_text(size= 40, margin=margin(r = 30), angle=0,vjust = 0.55),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

tiff("Macrophages_vlnplot.jpeg", unit="in", width=20, height=20, res=500)
C1qa/Ly6c2/Isg15/Mki67/Spp1/Ccr2
dev.off()

# https://www.science.org/doi/10.1126/sciimmunol.aay6017?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed
dot_features <- c("C1qa","C1qc", "Apoe","Cd74", "H2-Ab1", "H2-Aa", "Mki67", "Top2a", "Ccnb1","Spp1", "Mmp9", "Vegfa", "Arg1", 
                  "Ly6c2", "Ifit2","Isg15",
                  "Ccr2", "Chil3", "Il1b")
tiff("Macrophages_dotplot.jpeg", unit="in", width=10, height=4, res=500)
DotPlot(Macrophages, features = dot_features, group.by="SingleR.label") +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(1,12)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

# M1 & M2 cytokines and chemokines
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4188125/#:~:text=Those%20macrophages%20express%20a%20series,IL%2D6%20and%20express%20iNOS.

temp <- Macrophages
Idents(temp) <- "SingleR.label"

dot_features <- c("Tnf", "Il1b", "Il6", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl8", "Ccl9",
                  "Il10", "Tgfbi", "Ccl17", "Ccl24",
                  "Igf1", "Vegfa", "Irf4", "Mrc1", "Cd163")

DotPlot(temp, features = dot_features, split.by="SingleR.label") +
  scale_size(range = c(1,12)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

dittoDotPlot(temp, vars = dot_features, group.by = "treatment", split.by="SingleR.label", split.nrow = 5,
             size=10, max.color="red", min.color = "yellow") + theme(axis.title.y=element_blank(), axis.text.x = element_text(angle=45))


tiff("Macrophages_celltype_proportion_treatment.jpeg", unit="in", width=9, height=6, res=500)
plot_group_proportions(Macrophages, graph.type = "dodge")
dev.off()

tiff("Macrophages_celltype_proportion_treatment_stacked.jpeg", unit="in", width=3, height=6, res=500)
plot_group_proportions(Macrophages, graph.type = "stacked")
dev.off()

Macrophages_prop <- sc_utils(Macrophages)

prop_test <- permutation_test(
  Macrophages_prop, cluster_identity = "SingleR.label",
  sample_1 = "control", sample_2 = "TAK981",
  sample_identity = "treatment"
)
permutation_plot(prop_test)+ scale_colour_discrete(labels = function(x) str_wrap(x, width = 7))

tiff("Macrophages_celltype_proportion_test.jpeg", unit="in", width=5, height=1.8, res=500)
permutation_plot(prop_test)+ scale_colour_discrete(labels = function(x) str_wrap(x, width = 7))
dev.off()

tiff("Macrophages_umap.jpeg", unit="in", width=5.5, height=3.8, res=500)
DimPlot(Macrophages, reduction="umap",label=F,pt.size=1, group.by = "SingleR.label")
dev.off()

tiff("Macrophages_umap_split.jpeg", unit="in", width=8, height=3.8, res=500)
DimPlot(Macrophages, reduction="umap",label=F,pt.size=1, group.by = "SingleR.label", split.by = "treatment")
dev.off()

dot_features <- c("Cd80", "Cd86", "Sod2", "Nos2", "Stat1","Il1b","Il18", "Cxcl9", "Cxcl10",
                  "Mrc1", "Msr1","Chil3", "Retnla", "Arg1", "Stat3",
                  "Il10ra","Cd163", "Tgfb1","Mertk", "Vegfa")

tiff("M1 M2 Macrophages.jpeg", unit="in", width=12, height=3.5, res=500)
DotPlot(Macrophages, features = dot_features, group.by="treatment") +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(1,10)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

FeaturePlot(Macrophages, features=c("Spp1"))

################################################################################################
hallmark_pathway <- gmtPathways("h.all.v7.5.1.symbols.gmt")
temp <- Macrophages
temp[["orig.ident"]] <- replace(temp[["orig.ident"]], temp[["orig.ident"]]=="C1", "control")
temp[["orig.ident"]] <- replace(temp[["orig.ident"]], temp[["orig.ident"]]=="C2", "control")
temp[["orig.ident"]] <- replace(temp[["orig.ident"]], temp[["orig.ident"]]=="T1", "TAK981")
temp[["orig.ident"]] <- replace(temp[["orig.ident"]], temp[["orig.ident"]]=="T2", "TAK981")
Idents(temp) <- "orig.ident" #set all idents to orig.ident i.e. control, TAK981, etc

Macrophages_ranked_list2 <- FindMarkers(temp, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)
# order list, pull out gene name and log2fc, and convert genes to uppercase
Macrophages_ranked_list2 <- Macrophages_ranked_list2[order(Macrophages_ranked_list2$avg_log2FC, decreasing = T),]
Macrophages_ranked_list2$Gene.name <- str_to_upper(rownames(Macrophages_ranked_list2))
Macrophages_ranked_list2 <- Macrophages_ranked_list2[,c("Gene.name", "avg_log2FC")]
rownames(Macrophages_ranked_list2) <- NULL
Macrophages_ranked_list2 <- prepare_ranked_list(Macrophages_ranked_list2)
head(Macrophages_ranked_list2)

Macrophages_fgsea_results <- fgseaMultilevel(pathways = hallmark_pathway,
                                     stats = Macrophages_ranked_list2,
                                     minSize = 15,
                                     maxSize = Inf)

Macrophages_fgsea_df<- Macrophages_fgsea_results %>% arrange (desc(NES)) %>% select (pathway, pval, padj, NES, size)
write.csv(Macrophages_fgsea_df, "Macrophages_fgsea.csv",row.names = FALSE)

topPathwaysUp <- Macrophages_fgsea_results[ES > 0][head(order(pval), n=6), pathway]
topPathwaysDown <- Macrophages_fgsea_results[ES < 0][head(order(pval), n=6), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(hallmark_pathway[topPathways],
              stats=Macrophages_ranked_list2,
              fgseaRes=Macrophages_fgsea_results, 
              gseaParam = 0.5)

tiff("Macrophages_gseatable.jpeg", unit="in", width=10, height=3.5, res=500)
plotGseaTable(hallmark_pathway[topPathways],
              stats=Macrophages_ranked_list2,
              fgseaRes=Macrophages_fgsea_results, 
              gseaParam = 0.5)
dev.off()

################################################################################################
#gene correlation 
library(ggpubr)
DefaultAssay(Macrophages) <- "RNA"

sumo_features <- list(c('UBE2I'))
Macrophages <- AddModuleScore(
  Macrophages,
  features = sumo_features,
  name='SUMO_module'
)

M1_features <- list(c("Cd86", "Nos2",'Stat1', "Cxcl9", "Cxcl10"))
M2_features <- list(c("Cd163", "Arg1", "Mrc1","Stat3", "Msr1", "Vegfa"))

Macrophages <- AddModuleScore(
  Macrophages,
  features = M1_features,
  name='M1_module'
)

Macrophages <- AddModuleScore(
  Macrophages,
  features = M2_features,
  name='M2_module'
)

Split_FeatureScatter(seurat_object = Macrophages, feature1 = "SUMO1", feature2 = "M21",
                     split.by = "treatment", group.by = "SingleR.label", num_columns = 2,
                     pt.size = 1)

SUMOvsM1 <- FetchData(Macrophages, c("M1_module1", "M2_module1", "treatment") , slot = "data")

my_comparisons <- list(c("control", "TAK981"))
ggboxplot(SUMOvsM1, x = "treatment", y = "M2_module1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)
ggboxplot(SUMOvsM1, x = "treatment", y = "M1_module1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)
ggboxplot(SUMOvsM1, x = "treatment", y = "Cd40", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)

tiff("M1_Macrophages.jpeg", unit="in", width=2.8, height=5, res=500)
ggboxplot(SUMOvsM1, x = "treatment", y = "M1_module1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)
dev.off()

tiff("M2_Macrophages.jpeg", unit="in", width=2.8, height=5, res=500)
ggboxplot(SUMOvsM1, x = "treatment", y = "M2_module1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)
dev.off()



TAM <- FetchData(Macrophages, c('Cd40', "treatment", "SingleR.label") , slot = "data")
my_comparisons <- list(c("control", "TAK981"))
ggboxplot(TAM, x = "treatment", y = "Cd40", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)

dot_features <- c("Cd40", "Cd47")
VlnPlot(Macrophages, features = c("Cd40"), group.by = "SingleR.label", split.by = "treatment", pt.size=0.01)
dittoDotPlot(Macrophages, vars = dot_features, group.by= "treatment", split.by = "SingleR.label", split.nrow = 8,
             size=10, max.color="red", min.color = "yellow") +
  theme(axis.title.y=element_blank(), axis.text.x = element_text(angle=30))

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
setwd("D:/control_vs_TAK981")
library("scProportionTest")
source("scFunctions.R")

mref <- ImmGenData() #mouse immune cells
sceM <- BaronPancreasData('mouse')
sceM <- sceM[,!is.na(sceM$label)]
library(scuttle)
sceM <- logNormCounts(sceM)
## import all necessary files
## use the following command to unzip tar.gz: tar -xvzf
hallmark_pathway <- gmtPathways("c2.cp.kegg.v7.5.1.symbols")
load('nonimmune_new.RData')
memory.limit(56000)

#perform scTransform
nonimmune <- DietSeurat(nonimmune, assay="RNA")
nonimmune.list <- SplitObject(nonimmune, split.by = "orig.ident")

for (i in 1:length(nonimmune.list)) {
  nonimmune.list[[i]] <- SCTransform(nonimmune.list[[i]], vars.to.regress=c("percent.mito", "percent.ribo"))
}
features <- SelectIntegrationFeatures(object.list = nonimmune.list, nfeatures = 2000)
nonimmune.list <- PrepSCTIntegration(object.list = nonimmune.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = nonimmune.list, normalization.method = "SCT",
                                  anchor.features = features)
nonimmune <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
DefaultAssay(nonimmune) <- "integrated"
nonimmune <- ScaleData(object = nonimmune, features = features,verbose = FALSE)
nonimmune <- RunPCA(nonimmune, verbose = FALSE)
nonimmune <- RunTSNE(nonimmune, dims = 1:40, verbose = FALSE)
nonimmune <- RunUMAP(nonimmune, dims = 1:40, verbose = FALSE)
nonimmune <- FindNeighbors(nonimmune, reduction="umap",dims = 1:2, verbose = FALSE)
nonimmune <- FindClusters(nonimmune, resolution = 0.8,verbose = FALSE)
DimPlot(nonimmune, reduction="umap",label=TRUE,pt.size=1)

#Must normalize data with Lognormalize for further DE analyses
DefaultAssay(nonimmune) <- "RNA"
nonimmune <- NormalizeData(nonimmune,verbose = FALSE)
all.genes <- rownames(nonimmune)
nonimmune <- ScaleData(object = nonimmune, features = all.genes)
nonimmune.markers <- FindAllMarkers(nonimmune, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
nonimmune.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10
write.csv(nonimmune.markers, "nonimmune_markers2.csv", row.names = T)

tiff("nonimmune_cluster_heatmap2.jpeg", unit="in", width=30, height=18, res=500)
DoHeatmap(nonimmune, features = top10$gene) + NoLegend()
dev.off()

tiff("nonimmune_tsne_unlabeled.jpeg", unit="in", width=6, height=5, res=300)
DimPlot(nonimmune, reduction="tsne",label=TRUE,pt.size=1, repel = T)
dev.off()

DimPlot(nonimmune, reduction="umap",label=TRUE,pt.size=1)

temp <- nonimmune
SRKPC <- as.SingleCellExperiment(temp)
SRKPC <- SingleR(test=SRKPC, ref=sceM, assay.type.test = 1, labels = sceM$label)
temp[["SingleR.label"]] <- SRKPC$labels
rm(SRKPC)
DimPlot(nonimmune, reduction="umap",label=TRUE,pt.size=1,group.by = "SingleR.label")
# CAF <- subset(nonimmune, idents=c(0:16)[c(8,13,17)])
VlnPlot(nonimmune, features = c("Ly6c2"))


VlnPlot(nonimmune, features=c("Cd74"),group.by = "SingleR.label", split.by = "treatment")
FeaturePlot(nonimmune, features=c("Cd74"))



# temp <- nonimmune
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==6, "myCAF")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==17, "apCAF") # Saa3, Slpi 
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==22, "iCAF") #Clec3b, Has1, Col14a1
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==21, "PSC") #Acta2 (a-SMA) increase  
# nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==11, "Fibrocytes") #
# nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==25, "Fibrocytes") 

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==13, "Endothelial") #Pecam1

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==19, "ADM") #Acinar-to-ductal metaplasia area

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==23, "Ductal 3") # EpCAM+

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==0, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==1, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==5, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==8, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==9, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==11, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==12, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==20, "Ductal 1")

# nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==3, "mt")

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==3, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==4, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==7, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==10, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==14, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==15, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==16, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==18, "Ductal 2")

VlnPlot(kpc, features = c("Cd74"), group.by="SingleR.label", split.by='treatment')


dot_features <- c("Cd40")

DotPlot(kpc, features = dot_features, split.by = "treatment")

DotPlot(nonimmune, features = dot_features, group.by="SingleR.label") +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(1,10)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

kpc <- subset(nonimmune, subset=(SingleR.label=="ADM" | SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2" | 
                                   SingleR.label=="Ductal 3" | SingleR.label=="Endothelial" |
                                   SingleR.label=="iCAF"| SingleR.label=="PSC" | SingleR.label=="myCAF" |
                                   SingleR.label=="apCAF"))

DimPlot(kpc, reduction="umap",label=TRUE,pt.size=1.5, group.by = "SingleR.label")

Idents(kpc) <- "seurat_clusters"
Idents(kpc) <- "SingleR.label"
levels(kpc) <- c("ADM", "Ductal 1", "Ductal 2", "Ductal 3", "Endothelial", "apCAF", "iCAF", "myCAF", "PSC")
kpc$SingleR.label <- factor(x=kpc$SingleR.label)
levels(kpc$SingleR.label) <- c("ADM", "Ductal 1", "Ductal 2", "Ductal 3", "Endothelial", "apCAF", "iCAF", "myCAF", "PSC")

tiff("nonimmune_cluster_heatmap_SingleR.jpeg", unit="in", width=35, height=18, res=500)
DoHeatmap(kpc, features = top10$gene, group.by = "SingleR.label") + NoLegend()
dev.off()

tiff("nonimmune_umap.jpeg", unit="in", width=6, height=5, res=300)
DimPlot(kpc, reduction="umap", label=TRUE,pt.size=1, group.by = "SingleR.label",
        repel = TRUE)
dev.off()

tiff("nonimmune_umap_individual.jpeg", unit="in", width=15, height=5, res=300)
DimPlot(kpc, reduction="umap", label=TRUE,pt.size=1, split.by = "treatment", 
        repel = TRUE)
dev.off()


###########################################################################################################################
# Group-wise cell proportion plot (https://erilu.github.io/single-cell-rnaseq-analysis/#group-wise_analysis)
source("scFunctions.R")

tiff("nonimmune_celltype_proportion.jpeg", unit="in", width=9, height=6, res=300)
plot_group_proportions(kpc, graph.type = "dodge")
dev.off()

tiff("nonimmune_celltype_proportion_stacked.jpeg", unit="in", width=3, height=6, res=300)
plot_group_proportions(kpc, graph.type = "stacked")
dev.off()

#Proportion comparison - montecarlo
kpc_prop <- sc_utils(kpc)

prop_test <- permutation_test(
  kpc_prop, cluster_identity = "SingleR.label",
  sample_1 = "control", sample_2 = "TAK981",
  sample_identity = "treatment"
)
permutation_plot(prop_test)+ scale_colour_discrete(labels = function(x) str_wrap(x, width = 7)) 

tiff("kpc_celltype_proportion_test.jpeg", unit="in", width=7, height=3, res=500)
permutation_plot(prop_test)+ scale_colour_discrete(labels = function(x) str_wrap(x, width = 7))
dev.off()

####################################################################################################################
ductal <- subset(kpc, subset=(SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2"| SingleR.label=="Ductal 3" | SingleR.label=="ADM"))
CAF <- subset(kpc, subset=(SingleR.label=="myCAF" | SingleR.label=="iCAF"| SingleR.label=="apCAF"| SingleR.label=="PSC"))
endo <- subset(kpc, subset=(SingleR.label=="Endothelial"))

####################################################################################################################### 
#common fibroblast marker
tiff("kpc_Col1a1.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(kpc, features = c("Col1a1"),group.by = "SingleR.label")
dev.off()
#common CAF marker https://cancerdiscovery.aacrjournals.org/content/9/8/1102
tiff("kpc_Dcn.jpeg", unit="in", width=6, height=4, res=300)
VlnPlot(kpc, features = c("Onecut2"),group.by = "SingleR.label") #DCN is higher in iCAF
dev.off()

Prss2 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Prss2")) +
  ylab("Prss2") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Ctrb1 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Ctrb1")) +
  ylab("Ctrb1") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Reg1 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Reg1")) +
  ylab("Reg1") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Krt19 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Krt19")) +
  ylab("Krt19") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Onecut2 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Onecut2")) +
  ylab("Onecut2") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
S100a4 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("S100a4")) +
  ylab("S100a4") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Mki67 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Mki67")) +
  ylab("Mki67") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Hmgb2 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Hmgb2")) +
  ylab("Hmgb2") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Clu <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Clu")) +
  ylab("Clu") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Krt18 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Krt18")) +
  ylab("Krt18") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Pecam1 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Pecam1")) +
  ylab("Pecam1") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Cd34 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Cd34")) +
  ylab("Cd34") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
H2Ab1 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("H2-Ab1")) +
  ylab("H2-Ab1") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Saa3 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Saa3")) +
  ylab("Saa3") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Slpi <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Slpi")) +
  ylab("Slpi") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Col1a1 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Col1a1")) +
  ylab("Col1a1") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Col4a1 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Col4a1")) +
  ylab("Col4a1") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Has1 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Has1")) +
  ylab("Has1") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Postn <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Postn")) +
  ylab("Postn") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Thy1 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Thy1")) +
  ylab("Thy1") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Pdgfrb <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Pdgfrb")) +
  ylab("Pdgfrb") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Rgs5 <- VlnPlot(kpc,group.by = "SingleR.label", features = c("Rgs5")) +
  ylab("Rgs5") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_text(size = 25),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

tiff("kpc_vlnplot.jpeg", unit="in", width=24, height=28, res=500)
Prss2/Ctrb1/Krt19/Onecut2/S100a4/Mki67/Hmgb2/Clu/Krt18/Pecam1/Cd34/Col1a1/H2Ab1/Saa3/Has1/Postn/Thy1/Rgs5
dev.off()

dot_features <- c("Ifnb1")

tiff("kpc_marker_dotplot.jpeg", unit="in", width=12, height=7, res=500)
DotPlot(kpc, features = dot_features) +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(0,10)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


DotPlot(kpc, features = dot_features, group.by = "treatment")

genes <- c("Acta2", "Cald1", "Clu","Col1a1","Vim", "Snai1",
           "Cdh1", "Cdh3", "Epcam", "Egfr", "Klf4", "Krt18", "Krt19", "Met", "Tgfa")

# Annotating and ordering cells by some meaningful feature(s):
dittoHeatmap(ductal, genes,
             annot.by = c("SingleR.label", "treatment"))


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
ductal <- subset(nonimmune, subset=(SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2"))
ductal1 <- subset(nonimmune, subset=(SingleR.label=="Ductal 1"))
ductal2 <- subset(nonimmune, subset=(SingleR.label=="Ductal 2"))
CAF <- subset(nonimmune, subset=(SingleR.label=="iCAF"|SingleR.label=="apCAF"|SingleR.label=="myCAF"))
####################################################################################################################

#GSEA
hallmark_pathway <- gmtPathways("h.all.v7.4.symbols.gmt.txt")

kpc.subset <- ductal2
kpc.subset[["orig.ident"]] <- replace(kpc.subset[["orig.ident"]], kpc.subset[["orig.ident"]]=="C1", "control")
kpc.subset[["orig.ident"]] <- replace(kpc.subset[["orig.ident"]], kpc.subset[["orig.ident"]]=="C2", "control")
kpc.subset[["orig.ident"]] <- replace(kpc.subset[["orig.ident"]], kpc.subset[["orig.ident"]]=="T1", "TAK981")
kpc.subset[["orig.ident"]] <- replace(kpc.subset[["orig.ident"]], kpc.subset[["orig.ident"]]=="T2", "TAK981")
Idents(kpc.subset) <- "orig.ident" #set all idents to orig.ident i.e. control, TAK981, etc

kpc_ranked_list <- FindMarkers(kpc.subset, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)
temp_list <- kpc_ranked_list
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
write.csv(fgsea_df, "ductal2_fgsea.csv",row.names = FALSE)

topPathwaysUp <- fgsea_results[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgsea_results[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

kpc_gseatable <- plotGseaTable(hallmark_pathway[topPathways],
                                      stats=kpc_ranked_list,
                                      fgseaRes=fgsea_results, 
                                      gseaParam = 0.5)
tiff("kpc_gseatable.jpeg", unit="in", width=10, height=10, res=500)
plotGseaTable(hallmark_pathway[topPathways],
              stats=kpc_ranked_list,
              fgseaRes=fgsea_results, 
              gseaParam = 0.5)
dev.off()

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
ductal <- subset(kpc, subset=(SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2" ))

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

macrophage_gseatable <- plotGseaTable(hallmark_pathway[topPathways],
                                      stats=Macrophages_ranked_list,
                                      fgseaRes=Macrophages_fgsea_results, 
                                      gseaParam = 0.5)

tiff("macrophage_gseatable.jpeg", unit="in", width=10, height=3.5, res=500)
plotGseaTable(hallmark_pathway[topPathways],
              stats=Macrophages_ranked_list,
              fgseaRes=Macrophages_fgsea_results, 
              gseaParam = 0.5)
dev.off()


waterfall_plot(ductal_fgsea_results, "Pathways enriched in TAK981-treated vs control ductal cells")
# example of pathway highly enriched in treated
plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_ALPHA_RESPONSE" , ductal_ranked_list)
plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_GAMMA_RESPONSE" , ductal_ranked_list)

tiff("ductal_gsea.jpeg", unit="in", width=10, height=7, res=300)
waterfall_plot(ductal_fgsea_results, "Pathways enriched in TAK981-treated vs control ductal cells")
dev.off()



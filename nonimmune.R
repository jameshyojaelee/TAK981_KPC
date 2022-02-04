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
source("scFunctions.R")

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

#perform scTransform
nonimmune <- DietSeurat(nonimmune, assay="RNA")
nonimmune.list <- SplitObject(nonimmune, split.by = "orig.ident")

for (i in 1:length(nonimmune.list)) {
  nonimmune.list[[i]] <- SCTransform(nonimmune.list[[i]], vars.to.regress=c("nCount_RNA","percent.mito", "percent.ribo"))
}
features <- SelectIntegrationFeatures(object.list = nonimmune.list, nfeatures = 2000)
nonimmune.list <- PrepSCTIntegration(object.list = nonimmune.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = nonimmune.list, normalization.method = "SCT",
                                  anchor.features = features)
memory.limit(56000)
nonimmune <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
DefaultAssay(nonimmune) <- "integrated"
nonimmune <- ScaleData(object = nonimmune, features = features, vars.to.regress=c("nCount_RNA","percent.mito", "percent.ribo"),verbose = FALSE)
nonimmune <- RunPCA(nonimmune, verbose = FALSE)
nonimmune <- RunTSNE(nonimmune, dims = 1:40, verbose = FALSE)
nonimmune <- RunUMAP(nonimmune, dims = 1:40, verbose = FALSE)
nonimmune <- FindNeighbors(nonimmune, dims = 1:40, verbose = FALSE)
nonimmune <- FindClusters(nonimmune, resolution = 1.2,verbose = FALSE)
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
write.csv(nonimmune.markers, "nonimmune_markers.csv", row.names = T)

tiff("nonimmune_cluster_heatmap.jpeg", unit="in", width=30, height=18, res=500)
DoHeatmap(nonimmune, features = top10$gene) + NoLegend()
dev.off()

tiff("nonimmune_tsne_unlabeled.jpeg", unit="in", width=6, height=5, res=300)
DimPlot(nonimmune, reduction="tsne",label=TRUE,pt.size=1, repel = T)
dev.off()

DimPlot(nonimmune, reduction="tsne",label=TRUE,pt.size=1)

temp <- nonimmune
SRKPC <- as.SingleCellExperiment(temp)
SRKPC <- SingleR(test=SRKPC, ref=sceM, assay.type.test = 1, labels = sceM$label)
temp[["SingleR.label"]] <- SRKPC$labels
rm(SRKPC)
DimPlot(nonimmune, reduction="umap",label=TRUE,pt.size=1,group.by = "SingleR.label")
# CAF <- subset(nonimmune, idents=c(0:16)[c(8,13,17)])
VlnPlot(nonimmune, features = c("Cd74"))

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==7, "myCAF")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==12, "iCAF") #Clec3b, Has1, Col14a1
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==16, "PSC") #Acta2 (a-SMA) increase  
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==11, "apCAF")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==15, "apCAF") 

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==10, "Endothelial") #Pecam1

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==9, "ADM") #Acinar-to-ductal metaplasia area

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==14, "Ductal 3") # EpCAM+

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==0, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==1, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==4, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==5, "Ductal 1")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==6, "Ductal 1")

nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==2, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==3, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==4, "mt")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==8, "Ductal 2")
nonimmune[["SingleR.label"]]<- replace(nonimmune[["SingleR.label"]], nonimmune[["seurat_clusters"]]==13, "Ductal 2")

VlnPlot(nonimmune, features = c("Mki67"), group.by="SingleR.label")

dot_features <- c("Gata6", "Cldn18", "Tff1", "S100a2", "Krt17")
DotPlot(nonimmune, features = dot_features, group.by="SingleR.label") +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(1,10)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


kpc <- subset(nonimmune, subset=(SingleR.label=="ADM" | SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2" | 
                                   SingleR.label=="Ductal 3" | SingleR.label=="Endothelial" |
                                   SingleR.label=="iCAF"| SingleR.label=="PSC" | SingleR.label=="myCAF" |
                                   SingleR.label=="Fibrocytes" | SingleR.label=="apCAF"))

DimPlot(kpc, reduction="umap",label=TRUE,pt.size=2, group.by = "SingleR.label")

tiff("nonimmune_cluster_heatmap_SingleR.jpeg", unit="in", width=35, height=18, res=500)
DoHeatmap(kpc, features = top10$gene, group.by = "SingleR.label") + NoLegend()
dev.off()

tiff("nonimmune_umap.jpeg", unit="in", width=6, height=5, res=300)
DimPlot(kpc, reduction="umap", label=TRUE,pt.size=1, group.by = "SingleR.label",
        repel = TRUE)
dev.off()

tiff("nonimmune_umap_individual.jpeg", unit="in", width=15, height=5, res=300)
DimPlot(kpc, reduction="umap", label=TRUE,pt.size=1, group.by = "SingleR.label", split.by = "treatment", 
        repel = TRUE)
dev.off()

###########################################################################################################################
ductal <- subset(kpc, subset=(SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2" | SingleR.label=="Ductal 3" | SingleR.label=="ADM"))
CAF <- subset(kpc, subset=(SingleR.label=="myCAF" | SingleR.label=="iCAF"| SingleR.label=="apCAF"))
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

Prss3 <- VlnPlot(kpc, features = c("Amy2a2"), group.by = "SingleR.label") +
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
Mki67 <- VlnPlot(kpc, features = c("Mki67"), group.by = "SingleR.label") +
  ylab("Mki67") +
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
Prss3/Krt19/Mki67/Epcam/Pecam1/Col1a1/Dcn/Acta2/Postn/Rgs5
dev.off()


Idents(kpc) <- "SingleR.label"
levels(kpc) <- c("ADM", "Ductal 1", "Ductal 2", "Ductal 3", "Endothelial", "apCAF", "iCAF", "myCAF", "PSC")

dot_features <- c("Prss2","Ctrb1","Reg1", "Krt19", "Onecut2", "S100a4", "Cldn18",
                  "Mki67","Pcna","Hmgb2","Aldh1a1","Clu","Krt18","Isg15","Pecam1","Cd34",
                  "Cd74","H2-Ab1","Col1a1","Col14a1","Has1","Postn","Thy1","Rgs5", "Pdgfrb")

tiff("kpc_marker_dotplot.jpeg", unit="in", width=12, height=7, res=500)
DotPlot(kpc, features = dot_features) +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(0,10)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
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

####################################################################################################################
# Group-wise cell proportion plot (https://erilu.github.io/single-cell-rnaseq-analysis/#group-wise_analysis)
source("scFunctions.R")

tiff("nonimmune_celltype_proportion.jpeg", unit="in", width=9, height=6, res=300)
plot_group_proportions(kpc, graph.type = "dodge")
dev.off()

tiff("nonimmune_celltype_proportion_stacked.jpeg", unit="in", width=3, height=6, res=300)
plot_group_proportions(kpc, graph.type = "stacked")
dev.off()

tiff("nonimmune_proportion_heatmap.jpeg", unit="in", width=8, height=4, res=300)
plot_heatmap_proportions(kpc, graph.type = "by.cell")
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

topPathwaysUp <- fgsea_results[ES > 0][head(order(pval), n=3), pathway]
topPathwaysDown <- fgsea_results[ES < 0][head(order(pval), n=3), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

kpc_gseatable <- plotGseaTable(hallmark_pathway[topPathways],
                                      stats=kpc_ranked_list,
                                      fgseaRes=fgsea_results, 
                                      gseaParam = 0.5)
tiff("kpc_gseatable.jpeg", unit="in", width=10, height=3.5, res=500)
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



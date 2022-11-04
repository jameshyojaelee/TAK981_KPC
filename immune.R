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
load("new_immune.RData")
use_python("C:/ProgramData/Miniconda3/python.exe")
py_config()

mref <- ImmGenData() #mouse immune cells
sceM <- BaronPancreasData('mouse')
sceM <- sceM[,!is.na(sceM$label)]
library(scuttle)
sceM <- logNormCounts(sceM)
memory.limit(130000)
immune <- DietSeurat(immune, assay="RNA")
immune.list <- SplitObject(immune, split.by = "orig.ident")

for (i in 1:length(immune.list)) {
  immune.list[[i]] <- SCTransform(immune.list[[i]], vars.to.regress=c("percent.mito", "percent.ribo"))
}
features <- SelectIntegrationFeatures(object.list = immune.list, nfeatures = 2000)
immune.list <- PrepSCTIntegration(object.list = immune.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = immune.list, normalization.method = "SCT",
                                  anchor.features = features)
immune <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
DefaultAssay(immune) <- "integrated"
immune <- ScaleData(object = immune, features=features)
immune <- RunPCA(immune, verbose = FALSE)
immune <- RunTSNE(immune, dims = 1:40, verbose = FALSE)
immune <- RunUMAP(immune, dims = 1:40,verbose = FALSE)
immune <- FindNeighbors(immune, dims = 1:40, verbose = FALSE)
immune <- FindClusters(immune, resolution = 0.8)
DimPlot(immune, reduction="umap",label=TRUE,pt.size=1)

# lognormalize for marker genes
DefaultAssay(immune) <- "RNA"
immune <- NormalizeData(immune) #LogNormalize
all.genes <- rownames(immune)
immune <- ScaleData(object = immune, features = all.genes)
immune.markers <- FindAllMarkers(immune, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
immune.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> immune.top10
write.csv(immune.markers, "immune_markers2.csv", row.names = T)
jpeg("immuune_cluster_heatmap2.jpeg", width=3000, height=2000)
DoHeatmap(immune, features = immune.top10$gene) + NoLegend()
dev.off()

tiff("immune_unlabeled.jpeg", unit="in", width=8, height=5, res=500)
DimPlot(immune, reduction="umap",label=TRUE,pt.size=0.6)
dev.off()

VlnPlot(immune, features = c("Cd19"), group.by = "SingleR.label")
VlnPlot(immune, features = c("Ccr2"))

VlnPlot(immune, features=c("Cd74"),group.by = "SingleR.label", split.by = "treatment")
FeaturePlot(immune, features=c("Cd40"))



DimPlot(immune, reduction="tsne",label=TRUE,pt.size=1, group.by = "SingleR.label", repel=T)
#cluster 13 is not immune cell. some fibroblasts and fibrocytes express cd45. 
temp <- immune
# immune <- temp
SRimmune <- as.SingleCellExperiment(temp)
SRimmune <- SingleR(test=SRimmune, ref=mref, assay.type.test = 1, labels = mref$label.main)
temp[["SingleR.label"]] <- SRimmune$labels
rm(SRimmune)
DimPlot(temp, reduction="tsne", label=TRUE,pt.size=1, group.by = "SingleR.label", repel = TRUE)
subset(temp, subset=(seurat_clusters==22))

DimPlot(immune, reduction="umap",label=TRUE,pt.size=1, repel=T)
VlnPlot(immune, features = c("Spp1"))

dot_features <- c("Cd40", "Cd47", "Cd80", "Cd86", "Stat1")
VlnPlot(immune, features = c("Spp1"), group.by = "SingleR.label", split.by = "treatment", pt.size=0.01)
dittoDotPlot(immune, vars = dot_features, group.by= "treatment", split.by = "SingleR.label", split.nrow = 8,
             size=10, max.color="red", min.color = "yellow") +
  theme(axis.title.y=element_blank(), axis.text.x = element_text(angle=30))

temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==7, "B cells")
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==12, "T/NK cells")
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==14, "T/NK cells")
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==8, "Granulocytes") #Il1b+ Arg2+  https://www.science.org/doi/10.1126/sciimmunol.aay6017?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==20, "PMN-MDSC") # https://rupress.org/jem/article/218/4/e20201803/211778/Analysis-of-classical-neutrophils-and
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==15, "doublet")

# https://www.cell.com/pb-assets/products/nucleus/nucleus-phagocytes/rnd-systems-dendritic-cells-br.pdf
# https://www.science.org/doi/10.1126/scitranslmed.abf5058
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==11, "cDC2-Itgax") #expressing Itgax and Mgl2, Ear2, 
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==16, "cDC1-Ccl22") #Relb+ Itgae(Cd103)- CD8a-  Xcr1- Cd209a- Clec9a- Ly75(Cd205)+ 
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==17, "cDC1-Clec9a") #Itgae(Cd103)+ CD8a-  Xcr1+  Clec9a+ Batf3+
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==19, "pDC") #Cd11b- Siglech+ Bst2+ Tlr7+ Clec9a+ Ly6c2+ 

# M1 high tumor also proliferative
# https://www.nature.com/articles/s41598-020-73624-w
# https://www.nature.com/articles/s41467-017-01711-0
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==0, "Macro-MHCII") #MHC Class II high macropahge (cd74 is invariant MHCii) https://www.frontiersin.org/articles/10.3389/fimmu.2018.01132/full
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==2, "Macro-C1q")

temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==1, "Macro-MHCII") 
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==3, "Macro-MHCII")

temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==4, "Mono/MDSC") #high Il1b, Ly6c2, Ccr2. pro-tumorogensis genes such as Lyz2, Ifitm3, Vim, S100a6, https://www.jimmunol.org/content/206/1_Supplement/101.01
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==9, "Mono-Ly6c") #Ly6C+CCR2+ indicates BMDM. highly expressed genes are related to IFN response (IFIT, ISG, OAS). 

temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==6, "Macro-Proliferating") #Ki67 high

temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==5, "Macro-Spp1") #high in Arg1. may be MDSC. determine after InferCNV
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==10, "Macro-Spp1") 
temp[["SingleR.label"]] <- replace(temp[["SingleR.label"]], temp[["seurat_clusters"]]==13, "Macro-Spp1") 

#remove doublet cluster
immune <- subset(temp, subset=(SingleR.label=="Macro-Spp1" |SingleR.label=="Mono/MDSC" | SingleR.label=="Macro-MHCII" | SingleR.label=="Macro-C1q"  |
                                 SingleR.label=="Mono-Ly6c" |SingleR.label=="Macro-Proliferating" | SingleR.label=="pDC" | 
                                 SingleR.label=="cDC1-Clec9a" | SingleR.label=="cDC1-Ccl22" |SingleR.label== "cDC2-Itgax"|
                                 SingleR.label=="T/NK cells" | SingleR.label=="PMN-MDSC" | SingleR.label=="Granulocytes" |
                                 SingleR.label=="B cells"))
DimPlot(immune, reduction="umap", label=TRUE,pt.size=1, group.by = "SingleR.label", repel = TRUE) + theme(legend.position = "none")
# 
# Idents(immune) <- "SingleR.label"
# levels(immune) <- c("Macro-Spp1", "Mono/MDSC", "Macro-C1q", "Macro-Isg15", "Macro-Proliferating", "pDC", "cDC1-Clec9a",
#                     "cDC1-Ccl22", "cDC2-Itgax", "T/NK cells", "PMN-MDSC", "Granulocytes","B cells", "Basophils")

tiff("immune_umap_SingleR.jpeg", unit="in", width=8, height=6, res=500)
DimPlot(immune, reduction="umap", label=TRUE,pt.size=1, group.by = "SingleR.label", repel = TRUE) + theme(legend.position = "none")
dev.off()
tiff("immune_umap_SingleR_treatment.jpeg", unit="in", width=18, height=6, res=500)
DimPlot(immune, reduction="umap", label=TRUE,pt.size=1, group.by = "SingleR.label", split.by = "treatment", repel = TRUE) + theme(legend.position = "none")
dev.off()

tiff("immune_celltype_proportion_treatment.jpeg", unit="in", width=9, height=6, res=500)
plot_group_proportions(immune, graph.type = "dodge")
dev.off()

tiff("immune_celltype_proportion_treatment_stacked.jpeg", unit="in", width=3, height=6, res=500)
plot_group_proportions(immune, graph.type = "stacked")
dev.off()

G0s2 <- VlnPlot(immune, features = c("G0s2"), group.by = "SingleR.label") +
  ylab("G0s2") +
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

H2Ab1 <- VlnPlot(immune, features = c("H2-Ab1"), group.by = "SingleR.label") +
  ylab("H2-Ab1") +
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
Chil3 <- VlnPlot(immune, features = c("Chil3"), group.by = "SingleR.label") +
  ylab("Chil3") +
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
Ccr2 <- VlnPlot(immune, features = c("Ccr2"), group.by = "SingleR.label") +
  ylab("Ccr2") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Siglech <- VlnPlot(immune, features = c("Siglech"), group.by = "SingleR.label") +
  ylab("Siglech") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
Ngp <- VlnPlot(immune, features = c("Ngp"), group.by = "SingleR.label") +
  ylab("Ngp") +
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

VlnPlot(immune, features = c("Cd40"), group.by = "SingleR.label", split.by = "treatment") +
  ylab("Cd40") +
  scale_y_continuous(position="left")+
  theme(axis.title.x =element_blank(),
        axis.text.x = element_text(size = 25),
        axis.title.y.left = element_text(size= 25, margin=margin(r = 30)),
        legend.position = 'none',
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

tiff("immune_vlnplot.jpeg", unit="in", width=24, height=20, res=500)
Cd79a/Ccl22/Clec9a/Itgax/G0s2/C1qa/H2Ab1/Mki67/Spp1/Ly6c2/Ccr2/Siglech/Ngp/Cd3e
dev.off()


dot_features <- c("Cd79a","Ccl22","Clec9a","Itgax","G0s2","C1qc","","Mki67","Spp1","Ly6c2","Isg15","Ccr2","Chil3","Siglech","Ngp","Cd3e", "Cd27")

tiff("immune_marker_dotplot.jpeg", unit="in", width=10, height=8, res=500)
DotPlot(immune, features = dot_features, group.by="SingleR.label") +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(0,10)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

VlnPlot(immune, features = c("Ccl5"), group.by = "SingleR.label", split.by = "treatment")
temp <- immune
Idents(temp) <- "seurat_clusters"
Idents(temp) <- "SingleR.label"

tiff("immune_CD40.jpeg", unit="in", width=8, height=6, res=500)
FeaturePlot(immune, feature=c("Cd40"),pt.size=0.7)
dev.off()

dot_features <- c("Ccr2","Chil3","Vcan","Spp1","Arg1","Vegfa", "Mki67","Top2a","Ube2c","Ly6c2","Isg15", "Plac8" ,"C1qb","Mafb", "Maf") 
dittoDotPlot(immune, vars = dot_features, group.by = "SingleR.label",size=10, max.color="red", min.color = "yellow") +
  theme(axis.title.y=element_blank(),
        axis.text.x = element_text(angle=30)) 

DotPlot(TNK, features = c("Tnfaip3"), group.by = "treatment") + scale_size(range = c(1,10))

gene_sig <- c("Batf3")
comparisons <- list(c("control", "TAK981"))
vp_case1(gene_signature = gene_sig, file_name = "test", test_sign = comparisons, y_max = 7)


#Proportion comparison - montecarlo
library("scProportionTest")
immune_prop <- sc_utils(immune)

prop_test <- permutation_test(
  immune_prop, cluster_identity = "SingleR.label",
  sample_1 = "control", sample_2 = "TAK981",
  sample_identity = "treatment"
)
permutation_plot(prop_test)+ scale_colour_discrete(labels = function(x) str_wrap(x, width = 7))

tiff("immune_celltype_proportion_test.jpeg", unit="in", width=7, height=3, res=500)
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
Macrophages <- subset(immune, subset=(SingleR.label=="Macro-Spp1" |SingleR.label=="Mono/MDSC" |
                                    SingleR.label=="Macro-Ly6c" |SingleR.label=="Macro-Proliferating" | SingleR.label=="Macro-C1q"))
DimPlot(Macrophages, reduction="umap", label=TRUE,pt.size=1, group.by = "SingleR.label", repel = TRUE)
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
hallmark_pathway <- gmtPathways("h.all.v7.4.symbols.gmt.txt")

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

topPathwaysUp <- Macrophages_fgsea_results[ES > 0][head(order(pval), n=3), pathway]
topPathwaysDown <- Macrophages_fgsea_results[ES < 0][head(order(pval), n=3), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

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

waterfall_plot(Macrophages_fgsea_results, "Pathways enriched in TAK981 treated vs control Macrophages")

tiff("waterfall_control_vs_TAK981_macrophages.jpeg", unit="in", width=9, height=6, res=300)
waterfall_plot(Macrophages_fgsea_results, "Pathways enriched in TAK981-treated vs control macrophages")
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
#Boxplot module
DefaultAssay(Macrophages) <- "RNA"
FeatureScatter(Macrophages, "SUMO2", "CD163", slot="counts")

sumo_features <- list(c('Sumo1','Sumo2','Sumo3','Ube2i', 'Uba2'))
Macrophages <- AddModuleScore(
  Macrophages,
  features = sumo_features,
  name='SUMO_score'
)

M1_features <- list(c("Cd80", 'Cd86','Stat1', 'Ccl2', 'Ccl3', 'Ccl4','Cxcl10', "Nos2"))
M2_features <- list(c("Mrc1", "Chil3", "Arg1","Cd163", "Stat3","Stat6", "Msr1", "Tgfb1", "Vegfa"))

Macrophages <- AddModuleScore(
  Macrophages,
  features = M1_features,
  name='M1_score'
)

Macrophages <- AddModuleScore(
  Macrophages,
  features = M2_features,
  name='M2_score'
)


FeatureScatter(Macrophages, 'SUMO_score1', 'M2_score1', group.by = "treatment")

Split_FeatureScatter(seurat_object = Macrophages, feature1 = "SUMO1", feature2 = "M21",
                     split.by = "treatment", group.by = "SingleR.label", num_columns = 2,
                     pt.size = 1)

SUMOvsM1 <- FetchData(Macrophages, c('M1_score1', "SUMO_score1", "treatment") , slot = "data")
SUMOvsM2 <- FetchData(Macrophages, c('M2_score1', "SUMO_score1", "treatment") , slot = "data")

VlnPlot(Macrophages, features = c('SUMO_score1', 'M2_score1'), group.by = "treatment")

my_comparisons <- list(c("control", "TAK981"))
ggboxplot(SUMOvsM1, x = "treatment", y = "SUMO_score1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)
ggboxplot(SUMOvsM1, x = "treatment", y = "M1_score1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)

tiff("SUMO_Macrophages.jpeg", unit="in", width=4, height=6, res=500)
ggboxplot(SUMOvsM2, x = "treatment", y = "SUMO_score1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)
dev.off()

tiff("M1_Macrophages.jpeg", unit="in", width=4, height=6, res=500)
ggboxplot(SUMOvsM1, x = "treatment", y = "M1_score1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)
dev.off()
tiff("M2_Macrophages.jpeg", unit="in", width=4, height=6, res=500)
ggboxplot(SUMOvsM2, x = "treatment", y = "M2_score1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)
dev.off()

#3D featurescatter
data <- FetchData(Macrophages, vars = c('MRC1', 'CD86', 'UBE2I'))
ggplot(data = data) +
  geom_point(mapping = aes(x = MRC1, y = CD86, color = UBE2I), position = 'jitter') +
  scale_color_gradientn(colors = c('lightgrey', 'blue')) +
  cowplot::theme_cowplot()

#########################

MDSC_score <- list(c("Il1b", "Vegfa", "Stat3"))

Macrophages <- AddModuleScore(
  Macrophages,
  features = M1_features,
  name='MDSC_score'
)

FeatureScatter(Macrophages, 'SUMO_score1', 'MDSC_score1', group.by = "treatment")
SUMOvsMDSC <- FetchData(Macrophages, c('MDSC_score1', "SUMO_score1", "treatment") , slot = "data")

tiff("MDSC_score.jpeg", unit="in", width=4, height=6, res=500)
ggboxplot(SUMOvsMDSC, x = "treatment", y = "MDSC_score1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)
dev.off()





############################################################################################################
load("DC.RData")
DC <- subset(Myeloid, subset=SingleR.label=="DC")

dot_features <- c("Itgam","Itgax","H2-Ab1","Cd40", "Tlr7", "Tlr9")

DotPlot(DC, features = dot_features, group.by="orig.ident") + scale_size(range = c(1,15)) + RotatedAxis() + 
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")

TAN <- subset(immune, subset=SingleR.label=="Granulocytes")




#################################################################################################################
treatment_split <- SplitObject(immune, split.by = "treatment")
control <- treatment_split$control
TAK981 <- treatment_split$TAK981

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
control_cellchatDB <- CellChatDB.mouse 
control_cellchat@DB <- control_cellchatDB
control_cellchat <- subsetData(control_cellchat)
control_cellchat <- identifyOverExpressedGenes(control_cellchat)
control_cellchat <- identifyOverExpressedInteractions(control_cellchat)
control_cellchat <- projectData(control_cellchat, PPI.mouse)
control_cellchat <- computeCommunProb(control_cellchat)
# control_cellchat <- filterCommunication(control_cellchat, min.cells = 40)
control_cellchat <- computeCommunProbPathway(control_cellchat)
control_cellchat <- aggregateNet(control_cellchat)
groupSize <- as.numeric(table(control_cellchat@idents))
par(mfrow = c(2,2), xpd=TRUE)
netVisual_circle(control_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(control_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


TAK981_cellchat <- GetAssayData(TAK981, assay = "RNA", slot = "data") # normalized data matrix
Idents(TAK981) <- "SingleR.label"
labels <- Idents(TAK981)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
TAK981_cellchat <- createCellChat(object = TAK981_cellchat, meta = meta, group.by = "group")
TAK981_cellchat <- addMeta(TAK981_cellchat, meta = meta, meta.name = "labels")
TAK981_cellchat <- setIdent(TAK981_cellchat, ident.use = "labels")
levels(TAK981_cellchat@idents) 
groupSize <- as.numeric(table(TAK981_cellchat@idents)) 
TAK981_cellchatDB <- CellChatDB.mouse 
TAK981_cellchat@DB <- TAK981_cellchatDB
TAK981_cellchat <- subsetData(TAK981_cellchat)
TAK981_cellchat <- identifyOverExpressedGenes(TAK981_cellchat)
TAK981_cellchat <- identifyOverExpressedInteractions(TAK981_cellchat)
TAK981_cellchat <- projectData(TAK981_cellchat, PPI.mouse)
TAK981_cellchat <- computeCommunProb(TAK981_cellchat)
# TAK981_cellchat <- filterCommunication(TAK981_cellchat, min.cells = 10)
TAK981_cellchat <- computeCommunProbPathway(TAK981_cellchat)
TAK981_cellchat <- aggregateNet(TAK981_cellchat)
groupSize <- as.numeric(table(TAK981_cellchat@idents))
par(mfrow = c(2,2), xpd=TRUE)
netVisual_circle(TAK981_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(TAK981_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



# Comparison Analysis
group.new = levels(control_cellchat@idents)
TAK981_cellchat <- liftCellChat(TAK981_cellchat, group.new)
object.list <- list(control = control_cellchat, TAK981 = TAK981_cellchat)
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}
combined_cellchat <- mergeCellChat(object.list, add.names = names(object.list))
combined_cellchat
# Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(combined_cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(combined_cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(2,2), xpd=TRUE)
# netVisual_diffInteraction(combined_cellchat, weight.scale = T)
netVisual_diffInteraction(combined_cellchat, weight.scale = T, measure = "weight")

tiff("immune_interaction_web.jpeg", unit="in", width=6, height=6, res=500)
netVisual_diffInteraction(combined_cellchat, weight.scale = T, measure = "weight")
dev.off()

gg1 <- netVisual_heatmap(combined_cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(combined_cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg2

# Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

tiff("immune_interaction_weight_comparison.jpeg", unit="in", width=8, height=4, res=500)
patchwork::wrap_plots(plots = gg)
dev.off()

# Hierarchy plot
pathways.show <- c("PD-L1") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
vertex.receiver = c(5) # Left portion of hierarchy plot the shows signaling to dermal cells and right portion shows signaling to epidermal cells
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

netAnalysis_signalingChanges_scatter(combined_cellchat, idents.use = "T/NK cells")
gg1 <- netAnalysis_signalingChanges_scatter(combined_cellchat, idents.use = "T/NK cells")
gg2 <- netAnalysis_signalingChanges_scatter(combined_cellchat, idents.use = "myCAF")
gg3 <- netAnalysis_signalingChanges_scatter(combined_cellchat, idents.use = "apCAF")
patchwork::wrap_plots(plots = list(gg1,gg2,gg3))



#Part II: Identify the conserved and context-specific signaling pathways
#functional
combined_cellchat <- computeNetSimilarityPairwise(combined_cellchat, type = "functional")
combined_cellchat <- netEmbedding(combined_cellchat, type = "functional", umap.method = "uwot")
combined_cellchat <- netClustering(combined_cellchat, type = "functional")
netVisual_embeddingPairwise(combined_cellchat, type = "functional", label.size = 3.5)

#structural
combined_cellchat <- computeNetSimilarityPairwise(combined_cellchat, type = "structural")
combined_cellchat <- netEmbedding(combined_cellchat, type = "structural",umap.method = "uwot")
combined_cellchat <- netClustering(combined_cellchat, type = "structural")
netVisual_embeddingPairwise(combined_cellchat, type = "structural", label.size = 2)

# Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(combined_cellchat, type = "functional")

#Identify and visualize the conserved and context-specific signaling pathways
gg1 <- rankNet(combined_cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(combined_cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

unique(combined_cellchat@idents$joint)

#outgoing signal
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#incoming signal
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#Overal
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
#Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(combined_cellchat, thresh=0.05, sources.use = 2, targets.use = c(1,3,5,6),  comparison = c(1,2), angle.x = 45)

gg1 <- netVisual_bubble(combined_cellchat, sources.use = 5, targets.use = c(2,15,17,18,19,20,21,22,23),  comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(combined_cellchat, sources.use = 5, targets.use = c(2,15,17,18,19,20,21,22,23),  comparison = c(1,2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2



#Identify dysfunctional signaling by using differential expression analysis
pos.dataset = "TAK981"
features.name = pos.dataset
combined_cellchat <- identifyOverExpressedGenes(combined_cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(combined_cellchat, features.name = features.name)
net.up <- subsetCommunication(combined_cellchat, net = net, datasets = "TAK981",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(combined_cellchat, net = net, datasets = "control",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, combined_cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, combined_cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(combined_cellchat, pairLR.use = pairLR.use.up, sources.use = 1, targets.use = c(15,16,17,18,19,20,21,22,23), comparison = c(1, 2), 
                        angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(combined_cellchat, pairLR.use = pairLR.use.down, sources.use = 1, targets.use = c(15,16,17,18,19,20,21,22,23), comparison = c(1, 2),
                        angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 1, targets.use = c(5), slot.name = 'net', net = net.up, 
                     lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 1, targets.use = c(5), slot.name = 'net', net = net.down, 
                     lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))



# Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("CD80") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("CD80") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))


# Chord diagram
pathways.show <- c("CD86") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

# NET Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(object.list[[1]]@idents)
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}

#Net chord diagram - gene
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 1, targets.use = c(15:23), lab.cex = 0.5, title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
}

# compare all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2, 3, 4), targets.use = c(8,10),  title.name = paste0("Signaling received by Inflam.DC and .TC - ", names(object.list)[i]), legend.pos.x = 10)
}

# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1), targets.use = c(15:23),slot.name = "netP", title.name = paste0("Signaling pathways sending from fibroblast - ", names(object.list)[i]), legend.pos.x = 10)
}
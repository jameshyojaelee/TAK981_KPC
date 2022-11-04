#CAF cell analysis
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
setwd("X:/control_vs_TAK981")
library("scProportionTest")
source("scFunctions.R")
# 
# mref <- ImmGenData() #mouse immune cells
# sceM <- BaronPancreasData('mouse')
# sceM <- sceM[,!is.na(sceM$label)]
# library(scuttle)
# sceM <- logNormCounts(sceM)
# ## import all necessary files
## use the following command to unzip tar.gz: tar -xvzf
hallmark_pathway <- gmtPathways("h.all.v7.4.symbols.gmt.txt")
load('CAF.RData')

save(CAF, file="CAF.RData")

temp_CAF <- CAF
#perform scTransform
CAF <- DietSeurat(CAF, assay="RNA")
DefaultAssay(CAF) <- "RNA"
CAF <- SCTransform(CAF, method = "glmGamPoi")
CAF <- RunPCA(CAF, verbose = FALSE)
CAF <- RunTSNE(CAF, dims = 1:40, verbose = FALSE)
CAF <- RunUMAP(CAF, dims = 1:40, verbose = FALSE)
CAF <- FindNeighbors(CAF,dims=1:40, verbose = FALSE)
CAF <- FindClusters(CAF, resolution = 1,verbose = FALSE)
DimPlot(CAF, reduction="umap",label=TRUE,pt.size=2)

DefaultAssay(CAF) <- "RNA"
CAF <- NormalizeData(CAF) 
all.genes <- rownames(CAF)
CAF <- ScaleData(object = CAF, features = all.genes)
CAF.markers <- FindAllMarkers(CAF, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CAF.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

write.csv(CAF.markers, "CAF_markers.csv", row.names = T)

tiff("CAF_cluster_heatmap.jpeg", unit="in", width=30, height=18, res=500)
DoHeatmap(CAF, features = top10$gene) + NoLegend()
dev.off()

DoHeatmap(CAF, features = top10$gene, group.by = "SingleR.label") + NoLegend()

FeaturePlot(CAF, features=c("Cd44", "Col1a1", "Col1a1","Col4a1", "Col6a1"), split.by = "treatment")
VlnPlot(CAF, features=c("Cd44", "Col1a1", "Col1a1","Col4a1", "Col6a1"), group.by = "SingleR.label",split.by = "treatment")

VlnPlot(CAF, "Saa3")
VlnPlot(CAF, "Krt19", group.by = "SingleR.label")

DimPlot(CAF, reduction="umap",label=TRUE,pt.size=2)
CAF[["SingleR.label"]]<- replace(CAF[["SingleR.label"]], CAF[["seurat_clusters"]]==0, "apCAF") # Krt19
CAF[["SingleR.label"]]<- replace(CAF[["SingleR.label"]], CAF[["seurat_clusters"]]==1, "myCAF") # Pdgfrb, Postn
CAF[["SingleR.label"]]<- replace(CAF[["SingleR.label"]], CAF[["seurat_clusters"]]==2, "apCAF") # Krt19 Saa3, Krt8, H2-Ab1
CAF[["SingleR.label"]]<- replace(CAF[["SingleR.label"]], CAF[["seurat_clusters"]]==3, "myCAF") 
CAF[["SingleR.label"]]<- replace(CAF[["SingleR.label"]], CAF[["seurat_clusters"]]==4, "PSC") 
CAF[["SingleR.label"]]<- replace(CAF[["SingleR.label"]], CAF[["seurat_clusters"]]==5, "myCAF") # 
CAF[["SingleR.label"]]<- replace(CAF[["SingleR.label"]], CAF[["seurat_clusters"]]==6, "iCAF") # Igf1, Ly6a, Il6, Pdgfra, Clec3b


tiff("CAF_umap_split.jpeg", unit="in", width=6.5, height=3, res=500)
DimPlot(CAF, reduction="umap",label=TRUE,pt.size=2, group.by = "SingleR.label", split.by = "treatment")
dev.off()


tiff("CAF_umap.jpeg", unit="in", width=4, height=3, res=500)
DimPlot(CAF, reduction="umap",label=F,pt.size=2, group.by = "SingleR.label")
dev.off()

tiff("CAF_celltype_proportion.jpeg", unit="in", width=9, height=6, res=300)
plot_group_proportions(CAF, graph.type = "dodge")
dev.off()

tiff("CAF_celltype_proportion_stacked.jpeg", unit="in", width=3, height=6, res=300)
plot_group_proportions(CAF, graph.type = "stacked")
dev.off()

library("scProportionTest")
CAF_prop <- sc_utils(CAF)
prop_test <- permutation_test(
  CAF_prop, cluster_identity = "SingleR.label",
  sample_1 = "control", sample_2 = "TAK981",
  sample_identity = "treatment"
)

tiff("CAF_celltype_permutation.jpeg", unit="in", width=5, height=1.5, res=300)
permutation_plot(prop_test)+ scale_colour_discrete(labels = function(x) str_wrap(x, width = 7))
dev.off()

pt <- table(CAF$SingleR.label, CAF$treatment)
pt <- as.data.frame(pt)
library(RColorBrewer)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) + 
  coord_polar("y", start=0) +
  scale_fill_manual(values = rev(brewer.pal(3, "RdYlBu"))) +
  theme(legend.title = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

tiff("CAF_pie.jpeg", unit="in", width=3, height=1.5, res=800)
pie
dev.off()

df <- as.data.frame(pt)
df$Var1 <- factor(df$Var1)
df$Var2 <- factor(df$Var2) 
ggplot(data=pt, aes(x=" ", y=Freq, group=Var1, colour=Var1, fill=Var1)) +
  geom_bar(width = 1, stat = "identity",position = position_fill()) +
  coord_polar("y", start=0) + 
  facet_grid(.~ Var2) +theme_void() + 
  theme(legend.title = element_blank())-> pie

# https://www.cell.com/cancer-cell/pdf/S1535-6108(21)00339-1.pdf
VlnPlot(CAF, features = c("H2-Ab1"), group.by = 'SingleR.label') #apCAF
VlnPlot(CAF, features = c("Pdgfra"), group.by = 'SingleR.label') #iCAF    induced by IL1/JAK-STAT3
VlnPlot(CAF, features = c("Acta2"), group.by = 'SingleR.label') #myCAF    induced by TGFb/SMAD2/3

VlnPlot(CAF, features = c("S100a4"), group.by = 'SingleR.label', split.by = "treatment")

dot_features <- c("Krt19","H2-Ab1","Clu","Slpi",
                  "Ly6c1", "Pdgfra","Igf1","Il6", "Cxcl1", "Ccl2", "Cxcl12",
                  "Pdgfrb", "Postn","Acta2", "Vim", "Fap",
                  "Rgs5","Col18a1")

tiff("CAF_dotplot.jpeg", unit="in", width=9, height=4, res=300)
dittoDotPlot(CAF, vars = dot_features, group.by = "SingleR.label",split.nrow = 4,
             size=12, max.color="red", min.color = "yellow") +
  theme(axis.title.y=element_blank(),
        axis.text.x = element_text(angle=30))
dev.off()

tiff("kpc_collagen_total.jpeg", unit="in", width=10, height=8, res=500)
FeaturePlot(CAF, features=c("Col1a1", 'Col1a2', 'Col4a1', 'Col6a1'))
dev.off()

########################################################################################################
#GSVA
hallmark_pathway <- gmtPathways("h.all.v7.4.symbols.gmt.txt")

kpc.subset <- CAF
Idents(kpc.subset) <- "SingleR.label" #set all idents to orig.ident i.e. control, TAK981, etc

kpc_ranked_list <- FindMarkers(kpc.subset, ident.1 = "myCAF", ident.2 =  "apCAF", min.pct = 0.1, logfc.threshold = 0)
# temp_list <- kpc_ranked_list
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
write.csv(fgsea_df, "myCAFvsapCAF_fgsea.csv",row.names = FALSE)

topPathwaysUp <- fgsea_results[ES > 0][head(order(pval), n=3), pathway]
topPathwaysDown <- fgsea_results[ES < 0][head(order(pval), n=7), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

kpc_gseatable <- plotGseaTable(hallmark_pathway[topPathways],
                               stats=kpc_ranked_list,
                               fgseaRes=fgsea_results, 
                               gseaParam = 0.5)

tiff("myCAFvsapCAF_gseatable.jpeg", unit="in", width=11, height=5, res=500)
plotGseaTable(hallmark_pathway[topPathways],
              colwidths = c(5, 1, 0.5, 0.5, 0),
              stats=kpc_ranked_list,
              fgseaRes=fgsea_results, 
              gseaParam = 0.5)
dev.off()

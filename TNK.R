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
library(EnhancedVolcano)
library(scCustomize)
library(ggpubr)
library(RColorBrewer)
library("scProportionTest")
setwd("X:/control_vs_TAK981")
source("scFunctions.R")
mref <- ImmGenData()
hallmark_pathway <- gmtPathways("h.all.v7.4.symbols.gmt.txt")

load("TNK.RData")

#TNK analysis
TNK <- DietSeurat(TNK, assay="RNA")
TNK <- SCTransform(TNK, method = "glmGamPoi")
TNK <- RunPCA(TNK, verbose = FALSE)
TNK <- RunTSNE(TNK, dims = 1:40, verbose = FALSE)
TNK <- RunUMAP(TNK, dims = 1:40,verbose = FALSE)
TNK <- FindNeighbors(TNK, dims = 1:40, verbose = FALSE)
TNK <- FindClusters(TNK, resolution = 1.4)
DimPlot(TNK, reduction="umap",label=TRUE,pt.size=2,split.by = "treatment")

DefaultAssay(TNK) <- "RNA"
TNK <- NormalizeData(TNK) 
all.genes <- rownames(TNK)
TNK <- ScaleData(object = TNK, features = all.genes)
TNK.markers <- FindAllMarkers(TNK, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TNK.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> TNK.top10
write.csv(TNK.markers, "TNK_markers.csv", row.names = T)

tiff("TNK_cluster_heatmap.jpeg", unit="in", width=12, height=12, res=500)
DoHeatmap(TNK, features = TNK.top10$gene, group.by = "SingleR.label") + NoLegend()
dev.off()

FeaturePlot(TNK,reduction = "umap", features = "Cd8")

SRTNK <- as.SingleCellExperiment(TNK)
SRTNK <- SingleR(test=SRTNK, ref=mref, assay.type.test = 1, labels = mref$label.main)
TNK[["SingleR.label"]] <- SRTNK$labels
rm(SRTNK)

TNK[["SingleR.label"]] <- replace(TNK[["SingleR.label"]], TNK[["seurat_clusters"]]==0, "CD4 TEM")
TNK[["SingleR.label"]] <- replace(TNK[["SingleR.label"]], TNK[["seurat_clusters"]]==1, "Teff1")
TNK[["SingleR.label"]] <- replace(TNK[["SingleR.label"]], TNK[["seurat_clusters"]]==2, "Naive-like")
TNK[["SingleR.label"]] <- replace(TNK[["SingleR.label"]], TNK[["seurat_clusters"]]==3, "Teff2")
TNK[["SingleR.label"]] <- replace(TNK[["SingleR.label"]], TNK[["seurat_clusters"]]==4, "Treg")
TNK[["SingleR.label"]] <- replace(TNK[["SingleR.label"]], TNK[["seurat_clusters"]]==5, "NK")

# TNK[["SingleR.label"]] <- replace(TNK[["SingleR.label"]], TNK[["SingleR.label"]]=="Tgd", "T cells")
# TNK[["SingleR.label"]] <- replace(TNK[["SingleR.label"]], TNK[["SingleR.label"]]=="ILC", "NK cells")
#TNK <- subset(TNK, subset=(SingleR.label=="T cells" |SingleR.label == "NK cells"))
Idents(TNK) <- "seurat_clusters"
DimPlot(TNK, reduction="umap", label=TRUE,pt.size=3)
DimPlot(TNK, reduction="umap", label=TRUE,pt.size=3, group.by = "SingleR.label")

tiff("TNK_UMAP_split.jpeg", unit="in", width=8, height=5, res=500)
DimPlot(TNK, reduction="umap",label=TRUE,pt.size=2, repel = T, group.by = "SingleR.label",split.by = "treatment")
dev.off()

tiff("TNK_UMAP_labeled.jpeg", unit="in", width=5, height=5, res=500)
DimPlot(TNK, reduction="umap", label=TRUE,pt.size=3, repel = T,group.by = "SingleR.label")
dev.off()


VlnPlot(TNK, features = c("Cd4"))

dot_features <- c("Cd3d", "Cd4","Sell", "Cd44","Il7r", "Cxcr4",
                  "Ccr7","Tcf7", "Lef1","Ccl5","Gzma","Gzmb", "Klrd1",
                  "Cd8a", "Cd8b1", "Nkg7", "Gzmk", "Havcr2", "Lag3", "Pdcd1", "Ctla4",
                  "Foxp3", "Il2ra")

tiff("TNK_dotplot.jpeg", unit="in", width=10, height=4, res=500)
DotPlot(TNK, features = dot_features, group.by="SingleR.label") +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(1,10)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

dot_features <- c("Cd69","Ctla4","Pdcd1", "Lag3", "Havcr2", "Tigit")

tiff("TNK_split_featureplot.jpeg", unit="in", width=8, height=16, res=500)
FeaturePlot(TNK, features = dot_features, split.by = "treatment")
dev.off()

#Differential expression Teff1 vs Teff2
Teff_markers <- FindMarkers(TNK, assay="RNA", ident.1=1, ident.2=3, min.pct = 0, logfc.threshold = 0)
EnhancedVolcano(Teff_markers, lab = rownames(Teff_markers),
                x = "avg_log2FC", y = "p_val",
                pCutoff = 1e-02, FCcutoff = 0.5, labSize = 5,
                legendPosition = "none", drawConnectors = F,
                title = "Teff1 vs Teff2")

dittoDotPlot(TNK, vars = dot_features, group.by= "treatment", split.by = "SingleR.label", split.nrow = 8,
             size=10, max.color="red", min.color = "yellow") +
  theme(axis.title.y=element_blank(), axis.text.x = element_text(angle=30))

tiff("TNK_celltype_proportion.jpeg", unit="in", width=6, height=4, res=300)
plot_group_proportions(TNK, graph.type = "dodge")
dev.off()

tiff("TNK_celltype_proportion_stacked.jpeg", unit="in", width=4.5, height=6, res=300)
plot_group_proportions(TNK, graph.type = "stacked")
dev.off()

library("scProportionTest")
TNK_prop <- sc_utils(TNK)
prop_test <- permutation_test(
  TNK_prop, cluster_identity = "SingleR.label",
  sample_1 = "control", sample_2 = "TAK981",
  sample_identity = "treatment"
)

tiff("TNK_celltype_permutation.jpeg", unit="in", width=5, height=1.5, res=300)
permutation_plot(prop_test)+ scale_colour_discrete(labels = function(x) str_wrap(x, width = 7))
dev.off()



TNK.temp <- subset(TNK)
Idents(TNK.temp) <- "treatment" #set all idents to orig.ident i.e. control, TAK981, etc

TNK_ranked_list <- FindMarkers(TNK.temp, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.5, logfc.threshold = 0,
                               test.use = "MAST")

write.csv(TNK_ranked_list, "TNK_DE_50percent.csv",row.names = T)


EnhancedVolcano(TNK_ranked_list,
                lab = rownames(TNK_ranked_list),
                x = 'avg_log2FC',
                y = 'p_val')

tiff("control_vs_TAK981_TNK.jpeg", unit="in", width=10, height=10, res=500)
EnhancedVolcano(TNK_ranked_list, lab = rownames(TNK_ranked_list),
                x = "avg_log2FC", y = "p_val",
                pCutoff = 0.05, FCcutoff = 0, labSize = 5, pointSize = 2,
                selectLab = c("Sell", "Cd44","Il7r",
                              "Ccr7", "Ccl5","Gzmb", "Klrd1", "Tnfrsf4", "Tnfrsf9", "Tnfrsf22", "Tnfrsf23",
                              "Cd8a", "Nkg7", "Havcr2", "Lag3", "Pdcd1", "Ctla4",
                              "Foxp3", "Il2ra", "Tigit", "Tox"),
                labCol = 'black',
                labFace = 'bold',
                parseLabels = TRUE,
                legendPosition = "none", 
                title = "Control vs TAK981 T/NK",
                drawConnectors = TRUE,
                widthConnectors = 0.6,
                arrowheads = FALSE,
                colConnectors = 'black',
                max.overlaps = 20)
 dev.off()

####################################################################################################################
#CD8 T cell analysis
NK <- subset(TNK, subset= Ncr1 > 0)
CTL <- subset(TNK, subset=(Cd3d > 0 | Cd3e > 0 | Cd3g >0))

DimPlot(CTL, reduction="tsne", label=TRUE,pt.size=3, group.by = "SingleR.label")

#Il2ra(Cd25), Tnfrsf4(Ox40), Tnfrsf9(41bb)
dot_features <- c("Cd28", "Cd69","Cd44","Tnfrsf4", "Tnfrsf9", "Nkg7", "Gzmb", 
                  "Ctla4","Pdcd1", "Tigit", "Lag3", "Havcr2")

tiff("Cd3+ Vlnplot.jpeg", unit="in", width=8, height=7, res=300)
VlnPlot(CTL, dot_features, group.by = "treatment")
dev.off()

tiff("Cd3+ actexh.jpeg", unit="in", width=10, height=3, res=300)
DotPlot(CTL, features = dot_features, group.by='treatment') +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(1,15)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()



dot_features <- c("Atf4", "Batf",
                  "Jun", "Junb", "Jund",
                  "Fos", "Fosb")
tiff("CTL Ap-1.jpeg", unit="in", width=6, height=2, res=300)
DotPlot(CTL, features = dot_features, group.by="treatment") +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(1,10)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

#NFAT pathway (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3694398/)
#Ppp3cc (Calcineurin)
#NFAT1 = NFATc2
#NFAT2= NFATc1
#NFAT4= NFATc3
#NFAT5

dot_features <- c("Cracr2a","Nfatc1", "Nfatc2", "Ppp3ca","Ppp3cc")
tiff("Cd8 NFAT.jpeg", unit="in", width=8, height=5, res=300)
DotPlot(cd8t, features = dot_features, group.by="treatment") + scale_size(range = c(1,20)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()

#NF-kb pathway
#IKK-alpha = CHUK
#IKK-beta = Ikbkb
#NEMO = Ikbkg

dot_features <- c("Nfkb1", "Nfkb2",
                  "Chuk","Ikbkb", "Ikbkg", "Rel")

tiff("Cd8 NFKB.jpeg", unit="in", width=8, height=5, res=300)
DotPlot(CTL, features = dot_features, group.by="treatment") + scale_size(range = c(1,20)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()

# Trac - TCR receptor alpha constant "Trac","Trbc1","Trbc2",
# Cd28 PI3k signaling
# Pik4ca - PI3k
# PRkcq PKC-theta
dot_features <- c("Cd28","Icos","Pik3ca", "Prkcq", "Mapk1",
                  "Plcg1","Grb2", "Vav1", "Card9")

tiff("PI3K CD28 signaling NfKB.jpeg", unit="in", width=8, height=5, res=300)
DotPlot(cd8t, features = dot_features, group.by="treatment") + scale_size(range = c(1,20)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
dev.off()

#NKT-like CD8t cells
#Cd49b = Itga2
#NKG2D = Klrk1
#CD94 = Klrd1
#NK1.1 = Klrb1
#iNKT1 genes - T-bet(Tbx21), Ifng, 
#iNKT2 genes: Il4, PLZF(Zbtb16)
#iNKT17 genes: Il17, ROR??t(Rorc),
#CD69 for activated NK cells
dot_features <- c("Ncam1", "Itga2","Klrk1","Klrd1", "Klrb1", "Klrb1a","Klrb1b", "Klrb1c", "Klrb1d",
                  "Trbc1", "Trac",
                  "Tbx21", "Ifng",
                  "Il4", "Zbtb16",
                  "Il17a","Rorc",
                  "Cd69", "Gzmb")

#tiff("PI3K CD28 signaling NfKB.jpeg", unit="in", width=8, height=5, res=300)
DotPlot(cd8t, features = dot_features, group.by="treatment") + scale_size(range = c(1,20)) + 
  RotatedAxis() + scale_colour_gradient2(low = "#FBF2AE", mid = "#FC8A09", high = "#D10F0F")
#dev.off()

#################################################################################################################
hallmark_pathway <- gmtPathways("h.all.v7.5.1.symbols.gmt")

CTL[["orig.ident"]] <- replace(CTL[["orig.ident"]], CTL[["orig.ident"]]=="C1", "control")
CTL[["orig.ident"]] <- replace(CTL[["orig.ident"]], CTL[["orig.ident"]]=="C2", "control")
CTL[["orig.ident"]] <- replace(CTL[["orig.ident"]], CTL[["orig.ident"]]=="T1", "TAK981")
CTL[["orig.ident"]] <- replace(CTL[["orig.ident"]], CTL[["orig.ident"]]=="T2", "TAK981")
Idents(CTL) <- "orig.ident" #set all idents to orig.ident i.e. control, TAK981, etc

CTL_ranked_list <- FindMarkers(CTL, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)
# order list, pull out gene name and log2fc, and convert genes to uppercase
CTL_ranked_list <- CTL_ranked_list[order(CTL_ranked_list$avg_log2FC, decreasing = T),]
CTL_ranked_list$Gene.name <- str_to_upper(rownames(CTL_ranked_list))
CTL_ranked_list <- CTL_ranked_list[,c("Gene.name", "avg_log2FC")]
rownames(CTL_ranked_list) <- NULL
CTL_ranked_list <- prepare_ranked_list(CTL_ranked_list)
head(CTL_ranked_list)

CTL_fgsea_results <- fgseaMultilevel(pathways = hallmark_pathway,
                                     stats = CTL_ranked_list,
                                     minSize = 15,
                                     maxSize = Inf)

CTL_fgsea_df<- CTL_fgsea_results %>% arrange (desc(NES)) %>% select (pathway, pval, padj, NES, size)
write.csv(CTL_fgsea_df, "CTL_fgsea.csv",row.names = FALSE)

topPathwaysUp <- CTL_fgsea_results[ES > 0][head(order(pval), n=3), pathway]
topPathwaysDown <- CTL_fgsea_results[ES < 0][head(order(pval), n=3), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(hallmark_pathway[topPathways],
             stats=CTL_ranked_list,
             fgseaRes=CTL_fgsea_results, 
             gseaParam = 0.5)

tiff("CTL_gseatable.jpeg", unit="in", width=10, height=2, res=500)
plotGseaTable(hallmark_pathway[topPathways],
              stats=CTL_ranked_list,
              fgseaRes=CTL_fgsea_results, 
              gseaParam = 0.5)
dev.off()


#################################################################################################################

#ProjecTILs
library(ProjecTILs)
ref <- load.reference.map("ref_TILAtlas_mouse_v1.rds")
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
DimPlot(ref,label = T, cols = refCols)

Tonly <- subset(TNK, subset=(Cd3d > 0 | Cd3e > 0))
DimPlot(Tonly, reduction="umap",label=TRUE,pt.size=1)

query.projected <- make.projection(Tonly, ref=ref)
plot.projection(ref, query.projected)

#Predict cell states
query.projected <- cellstate.predict(ref=ref, query=query.projected)
table(query.projected$functional.cluster)

plot.statepred.composition(ref, query.projected,metric = "Percent")

plot.states.radar(ref, query=query.projected)


# Compare states across treatments
TNK <- subset(immune, subset=(SingleR.label=="T/NK cells"))

query.control <- subset(query.projected, subset=(treatment=="control"))
query.perturb <- subset(query.projected, subset=(treatment=="TAK981"))


plot.statepred.composition(ref, query.control,metric = "Percent")
plot.statepred.composition(ref, query.perturb,metric = "Percent")

plot.states.radar(ref, query=query.control)
plot.states.radar(ref, query=query.perturb)


plot.states.radar(ref, query=list("Control" = query.control, "Query" = query.perturb))

discriminantGenes <- find.discriminant.genes(ref=ref, query=query.perturb,
                                             query.control=query.control, state="Treg")
head(discriminantGenes,n=10)


library(EnhancedVolcano)
EnhancedVolcano(discriminantGenes, lab = rownames(discriminantGenes),
                x = "avg_log2FC", y = "p_val",
                pCutoff = 1e-09, FCcutoff = 0.5, labSize = 5, 
                legendPosition = "none", drawConnectors = F, 
                title = "control vs. TAK981 (Treg)")

rand.list <- ProjecTILs:::randomSplit(query.projected, n=2, seed=1)
discriminantGenes <- find.discriminant.genes(ref=ref, query=rand.list[[1]], 
                                             query.control=rand.list[[2]], state="Treg")

EnhancedVolcano(discriminantGenes, lab = rownames(discriminantGenes),
                x = "avg_log2FC", y = "p_val",
                pCutoff = 1e-09, FCcutoff = 0.5, labSize = 5,
                legendPosition = "none", drawConnectors = F,
                title = "Random split (Treg)")

###########################################################################################################################
#Differential Expression

library(EnhancedVolcano)

DefaultAssay(CTL) <- "RNA"

CTL.temp <- CTL
Idents(CTL.temp) <- "treatment" #set all idents to orig.ident i.e. control, TAK981, etc
CTL_ranked_list <- FindMarkers(CTL.temp, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)
EnhancedVolcano(CTL_ranked_list,
                lab = rownames(CTL_ranked_list),
                x = 'avg_log2FC',
                y = 'p_val')




################################################################################################
#gene correlation 
DefaultAssay(TNK) <- "RNA"

ACT_features <- list(c("Cd69", "Cd28", "Il2ra", "Gzmb", "Prf1", "Ifng", "Cd44", "Tnfrsf4", "Tnfrsf9"))
TNK <- AddModuleScore(TNK, features = ACT_features, name='Act_score')

Tex_features <- list(c('Pdcd1', "TNKa4", "Havcr2", "Lag3", "Tigit"))
TNK <- AddModuleScore(TNK, features = Tex_features,name='Tex_score')

Tex_score <- FetchData(TNK, c('Tex_score1', "Act_score1", "treatment") , slot = "data")

my_comparisons <- list(c("control", "TAK981"))
ggplot(Tex_score, aes(x=treatment, y=Tex_score1, fill= treatment)) + 
  geom_violin(width=0.5) +
  geom_boxplot(width=0.1, color="black") +
  stat_compare_means(comparisons = my_comparisons) + 
  theme(axis.title.x =element_blank(),
        panel.background = element_rect(fill="white", colour='grey'),
        panel.grid.major = element_line(color="grey")) +
  scale_fill_manual(values=c("#56B4E9","#E69F00", "#8b0000")) -> p
p
tiff("Tex_score.jpeg", unit="in", width=5, height=6, res=500)
p
dev.off()

####################################################################################################################
DimPlot(CTL, split.by = "treatment")

df <- as.data.frame(CTL@assays$RNA@data)
write.csv(df, file="TNK_matrix.csv")

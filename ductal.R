#ductal cell analysis
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
load("ductal.RData")
#perform scTransform
ductal <- DietSeurat(ductal, assay="RNA")
ductal.list <- SplitObject(ductal, split.by = "orig.ident")

for (i in 1:length(ductal.list)) {
  ductal.list[[i]] <- SCTransform(ductal.list[[i]], vars.to.regress=c("percent.mito", "percent.ribo"))
}
features <- SelectIntegrationFeatures(object.list = ductal.list, nfeatures = 3000)
ductal.list <- PrepSCTIntegration(object.list = ductal.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = ductal.list, normalization.method = "SCT",
                                  anchor.features = features)
ductal <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
DefaultAssay(ductal) <- "integrated"
ductal <- ScaleData(object = ductal, features = features,verbose = FALSE)
ductal <- RunPCA(ductal, verbose = FALSE)
ductal <- RunTSNE(ductal, dims = 1:40, verbose = FALSE)
ductal <- RunUMAP(ductal, dims = 1:40, verbose = FALSE)
ductal <- FindNeighbors(ductal, reduction='umap',dims=1:2, verbose = FALSE)
ductal <- FindClusters(ductal, resolution = 0.2,verbose = FALSE)
DimPlot(ductal, reduction="umap",label=TRUE,pt.size=1)
DimPlot(ductal, reduction="umap",label=TRUE,pt.size=1, group.by = "SingleR.label")

DefaultAssay(ductal) <- "RNA"
ductal <- NormalizeData(ductal,verbose = FALSE)
all.genes <- rownames(ductal)
ductal <- ScaleData(object = ductal, features = all.genes)
ductal.markers <- FindAllMarkers(ductal, assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ductal.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> ductal.top10

tiff("ductal_umap.jpeg", unit="in", width=5, height=3.5, res=500)
DimPlot(ductal, reduction="umap",pt.size=1,group.by = "SingleR.label")
dev.off()

tiff("ductal_cluster_heatmap_SingleR.jpeg", unit="in", width=35, height=18, res=500)
DoHeatmap(ductal, features = ductal.top10$gene, group.by = "SingleR.label") + NoLegend()
dev.off()

ductal1 <- subset(ductal, subset=(SingleR.label=="Ductal 1"))
ductal2 <- subset(ductal, subset=(SingleR.label=="Ductal 2"))

VlnPlot(ductal, "Epcam")
ductal[["SingleR.label"]] <- replace(ductal[["SingleR.label"]], ductal[["seurat_clusters"]]==0, "Ductal 1") #
ductal[["SingleR.label"]] <- replace(ductal[["SingleR.label"]], ductal[["seurat_clusters"]]==1, "Ductal 2")
ductal[["SingleR.label"]] <- replace(ductal[["SingleR.label"]], ductal[["seurat_clusters"]]==2, "Ductal 1") 
ductal[["SingleR.label"]] <- replace(ductal[["SingleR.label"]], ductal[["seurat_clusters"]]==3, "Ductal 2") #Ki67 high
ductal[["SingleR.label"]] <- replace(ductal[["SingleR.label"]], ductal[["seurat_clusters"]]==4, "Ductal 1") #
ductal[["SingleR.label"]] <- replace(ductal[["SingleR.label"]], ductal[["seurat_clusters"]]==6, "Ductal 2") #
ductal[["SingleR.label"]] <- replace(ductal[["SingleR.label"]], ductal[["seurat_clusters"]]==5, "Ductal 1")
ductal[["SingleR.label"]] <- replace(ductal[["SingleR.label"]], ductal[["seurat_clusters"]]==7, "Ductal 3")
ductal[["SingleR.label"]] <- replace(ductal[["SingleR.label"]], ductal[["seurat_clusters"]]==8, "ADM")
ductal[["SingleR.label"]] <- replace(ductal[["SingleR.label"]], ductal[["seurat_clusters"]]==9, "Ductal 1")

dot_features <- c("Prss2","Ctrb1","Reg1", "Krt19", "Onecut2", "S100a4",
                  "Mki67","Pcna","Hmgb2","Epcam","Krt18","Isg15")

tiff("ductal_marker_dotplot.jpeg", unit="in", width=8, height=2.5, res=500)
DotPlot(ductal, features = dot_features, group.by = "SingleR.label") +
  scale_colour_gradient2(low = "#FBF2AE", mid = "#FFB809", high = "#D10F0F") + 
  scale_size(range = c(1,10)) + RotatedAxis() + theme(axis.title = element_blank())+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


FeaturePlot(ductal, features=c("Isg15"))
VlnPlot(ductal, features=c(""), group.by = "SingleR.label")
VlnPlot(ductal, features=c("Isg15"), group.by = "SingleR.label",split.by = "treatment")


########################################################################################################################
# Group-wise cell proportion plot (https://erilu.github.io/single-cell-rnaseq-analysis/#group-wise_analysis)
source("scFunctions.R")

tiff("ductal_celltype_proportion.jpeg", unit="in", width=9, height=6, res=300)
plot_group_proportions(ductal, graph.type = "dodge")
dev.off()

tiff("ductal_celltype_proportion_stacked.jpeg", unit="in", width=3, height=6, res=300)
plot_group_proportions(ductal, graph.type = "stacked")
dev.off()

#Proportion comparison - montecarlo
ductal_prop <- sc_utils(ductal)

prop_test <- permutation_test(
  ductal_prop, cluster_identity = "SingleR.label",
  sample_1 = "control", sample_2 = "TAK981",
  sample_identity = "treatment"
)
permutation_plot(prop_test)+ scale_colour_discrete(labels = function(x) str_wrap(x, width = 7))

tiff("ductal_celltype_proportion_test.jpeg", unit="in", width=5, height=1.5, res=500)
permutation_plot(prop_test)+ scale_colour_discrete(labels = function(x) str_wrap(x, width = 7))
dev.off()

####################################################################################################################
temp <- subset(ductal, subset=(SingleR.label=="Ductal 1"|SingleR.label=="Ductal 2"))
Idents(temp) <- "SingleR.label"

dot_features <- c("Cdh2", "Vim", "Fn1", "S100a4",
                  "Mki67", "Top2a",
                  "Ccnb1", "Ccnb2",
                  "H2-Ab1", "Cd74",
                  "Foxq1", "Egr1", "Elf3", "Snai2", "Aldh1a1",
                  "Krt7", "S100a2", "Gata6","Tff1", "Cldn18")

dot_features <- c( "H2-D1", "H2-K1", "H2-Q4","H2-Aa", "H2-Ab1", "Cd74")
# DotPlot(temp, features = dot_features, split.by = "treatment") +
#   scale_size(range = c(1,10)) + RotatedAxis() + theme(axis.title = element_blank())+
#   geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
#   guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

tiff("ductal1_2_MHC.jpeg", unit="in", width=4, height=4, res=500)
dittoDotPlot(temp, vars = dot_features, group.by = "treatment", split.by="SingleR.label", split.nrow = 4,
             size=10, max.color="red", min.color = "yellow") + 
  theme(axis.title.y=element_blank(), axis.text.x = element_text(angle=45))
dev.off()


FeaturePlot(ductal, features=c("Foxq1"),pt.size=1.5,reduction="tsne", split.by = "treatment")
VlnPlot(ductal, features = c("Ifnb1"), group.by = "SingleR.label",  split.by = "treatment")
VlnPlot(ductal, features = c("Aldh1a1", "Pou2f1", "Nes"), group.by = "SingleR.label", split.by = "treatment")
# Epithelial vs Mesenchymal markers 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6777805/#sd
# epithelial markers: Ocln, Gjb1, and Tjp1, Cldn4, Cldn3, Epcam
# 
# # cancer stem cell markers: Cldn3, Epcam, Cdh1 (E-cadherin)
# tiff("ductal_Aldh1a1.jpeg", unit="in", width=6, height=4, res=300)
# VlnPlot(ductal, features = c("Aldh1a1"), group.by = "SingleR.label", split.by = "treatment")
# dev.off()
# tiff("ductal_Pou2f1.jpeg", unit="in", width=6, height=4, res=300)
# VlnPlot(ductal, features = c("Pou2f1"), group.by = "SingleR.label", split.by = "treatment")
# dev.off()
# tiff("ductal_Nes.jpeg", unit="in", width=6, height=4, res=300)
# VlnPlot(ductal, features = c("Nes"), group.by = "SingleR.label", split.by = "treatment")
# dev.off()
# 
# # Epithelial markers: Cldn3, Epcam, Cdh1 (E-cadherin)
# tiff("ductal_Cldn3.jpeg", unit="in", width=6, height=4, res=300)
# VlnPlot(ductal, features = c("Cldn3"), group.by = "SingleR.label", split.by = "treatment")
# dev.off()
# tiff("ductal_Epcam.jpeg", unit="in", width=6, height=4, res=300)
# VlnPlot(ductal, features = c("Epcam"), group.by = "SingleR.label", split.by = "treatment")
# dev.off()
# tiff("ductal_Cdh1.jpeg", unit="in", width=6, height=4, res=300)
# VlnPlot(ductal, features = c("Cdh1"), group.by = "SingleR.label", split.by = "treatment")
# dev.off()
# 
# 
# # mesenchymal markers: S100a4, Vim, Cdh2 (N-cadherin)
# tiff("ductal_S100a4.jpeg", unit="in", width=6, height=4, res=300)
# VlnPlot(ductal, features = c("S100a4"), group.by = "SingleR.label", split.by = "treatment")
# dev.off()
# tiff("ductal_Vim.jpeg", unit="in", width=6, height=4, res=300)
# VlnPlot(ductal, features = c("Vim"), group.by = "SingleR.label", split.by = "treatment")
# dev.off()
# tiff("ductal_Cdh2.jpeg", unit="in", width=6, height=4, res=300)
# VlnPlot(ductal, features = c("Cdh2"), group.by = "SingleR.label", split.by = "treatment")
# dev.off()
# 
# 
# #Proliferation markers: Mki67, Top2a
# tiff("ductal_Mki67.jpeg", unit="in", width=6, height=4, res=300)
# VlnPlot(ductal, features = c("Mki67"), group.by = "SingleR.label", split.by = "treatment")
# dev.off()
# tiff("ductal_Top2a.jpeg", unit="in", width=6, height=4, res=300)
# VlnPlot(ductal, features = c("Top2a"), group.by = "SingleR.label", split.by = "treatment")
# dev.off()
# 
# #Cell Cycle marker
# tiff("ductal_Ccnb2.jpeg", unit="in", width=6, height=4, res=300)
# VlnPlot(ductal, features = c("Ccnb2"), group.by = "SingleR.label", split.by = "treatment")
# dev.off()

#GSEA
hallmark_pathway <- gmtPathways("c2.cp.kegg.v7.5.1.symbols.gmt")

temp <- subset(ductal, subset=(SingleR.label=="Ductal 1" | SingleR.label=="Ductal 2"))

temp[["orig.ident"]] <- replace(temp[["orig.ident"]], temp[["orig.ident"]]=="C1", "control")
temp[["orig.ident"]] <- replace(temp[["orig.ident"]], temp[["orig.ident"]]=="C2", "control")
temp[["orig.ident"]] <- replace(temp[["orig.ident"]], temp[["orig.ident"]]=="T1", "TAK981")
temp[["orig.ident"]] <- replace(temp[["orig.ident"]], temp[["orig.ident"]]=="T2", "TAK981")
Idents(temp) <- "orig.ident" #set all idents to orig.ident i.e. control, TAK981, etc

temp_ranked_list <- FindMarkers(temp, ident.1 = "TAK981", ident.2 =  "control", min.pct = 0.1, logfc.threshold = 0)
# order list, pull out gene name and log2fc, and convert genes to uppercase
temp_ranked_list <- temp_ranked_list[order(temp_ranked_list$avg_log2FC, decreasing = T),]
temp_ranked_list$Gene.name <- str_to_upper(rownames(temp_ranked_list))
temp_ranked_list <- temp_ranked_list[,c("Gene.name", "avg_log2FC")]
rownames(temp_ranked_list) <- NULL
temp_ranked_list <- prepare_ranked_list(temp_ranked_list)
head(temp_ranked_list)

temp_fgsea_results <- fgseaMultilevel(pathways = hallmark_pathway,
                                        stats = temp_ranked_list,
                                        minSize = 15,
                                        maxSize = Inf)

temp_fgsea_df<- temp_fgsea_results %>% arrange (desc(NES)) %>% select (pathway, pval, padj, NES, size)
write.csv(temp_fgsea_df, "ductal1AND2_fgsea.csv",row.names = FALSE)

topPathwaysUp <- temp_fgsea_results[ES > 0][head(order(pval), n=13), pathway]
topPathwaysDown <- temp_fgsea_results[ES < 0][head(order(pval), n=6), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(hallmark_pathway[topPathways],
              colwidths = c(5, 1, 0.8, 1.2, 1.2),
              stats=temp_ranked_list,
              fgseaRes=temp_fgsea_results, 
              gseaParam = 0.5)

tiff("ductal3_gseatable.jpeg", unit="in", width=10, height=3.5, res=500)
plotGseaTable(hallmark_pathway[topPathways],
              colwidths = c(5, 1, 0.8, 1.2, 1.2),
              stats=temp_ranked_list,
              fgseaRes=temp_fgsea_results, 
              gseaParam = 0.5)
dev.off()


########################################################################################################
#GSVA
hallmark_pathway <- gmtPathways("h.all.v7.4.symbols.gmt.txt")

kpc.subset <- ductal
Idents(kpc.subset) <- "SingleR.label" #set all idents to orig.ident i.e. control, TAK981, etc

kpc_ranked_list <- FindMarkers(kpc.subset, ident.1 = "Ductal 1", ident.2 =  "Ductal 2", min.pct = 0.1, logfc.threshold = 0)
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
write.csv(fgsea_df, "ductal1vs2_fgsea.csv",row.names = FALSE)

topPathwaysUp <- fgsea_results[ES > 0][head(order(pval), n=13), pathway]
topPathwaysDown <- fgsea_results[ES < 0][head(order(pval), n=6), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

kpc_gseatable <- plotGseaTable(hallmark_pathway[topPathways],
                               stats=kpc_ranked_list,
                               fgseaRes=fgsea_results, 
                               gseaParam = 0.5)

tiff("ductal1vs2_gseatable.jpeg", unit="in", width=11, height=5, res=500)
plotGseaTable(hallmark_pathway[topPathways],
              colwidths = c(5, 1, 0.5, 0.5, 0),
              stats=kpc_ranked_list,
              fgseaRes=fgsea_results, 
              gseaParam = 0.5)
dev.off()

###############################################################################
pt <- table(ductal$SingleR.label, ductal$treatment)
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
        axis.text.y = element_blank()) -> pie

tiff("ductal_pie.jpeg", unit="in", width=5, height=4, res=800)
pie
dev.off()


df <- as.data.frame(pt)
df$Var1 <- factor(df$Var1)
df$Var2 <- factor(df$Var2) 
ggplot(data=pt, aes(x=" ", y=Freq, group=Var1, colour=Var1, fill=Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  facet_grid(.~ Var2) +theme_void()



################################################################################################
#gene correlation 

library(ggpubr)
DefaultAssay(ductal) <- "RNA"

sumo_features <- list(c('Ube2i'))
ductal <- AddModuleScore(
  ductal,
  features = sumo_features,
  name='SUMO_module'
)

classical <- list(c("Krt20","Tff1", "Tff2",'Tspan8', "Ctse","Vsig2", "Lgals4")) #Lyz, Ceacam6 not there
basal <- list(c("Spp1", "Clu", "Mmp7","Col18a1", "Slc2a1", "Krt7", "Krt17", "S100a2", "Dhrs9")) #Ctgf not there

ductal <- AddModuleScore(
  ductal,
  features = classical,
  name='classical'
)

ductal <- AddModuleScore(
  ductal,
  features = basal,
  name='basal'
)


classic_basal <- FetchData(ductal, c('classical1', 'basal1', "treatment", "SingleR.label") , slot = "data")

VlnPlot(ductal, features = c('classical1', 'basal1'), group.by = "SingleR.label")

my_comparisons <- list(c("control", "TAK981"))
ggboxplot(classic_basal, x = "treatment", y = "classical1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)
ggboxplot(classic_basal, x = "treatment", y = "basal1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)

tiff("M1_ductal.jpeg", unit="in", width=4, height=6, res=500)
ggboxplot(SUMOvsM1, x = "treatment", y = "M1_module1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)
dev.off()

tiff("M2_ductal.jpeg", unit="in", width=4, height=6, res=500)
ggboxplot(SUMOvsM2, x = "treatment", y = "M2_module1", color='treatment', palette = "jco") + stat_compare_means(comparisons = my_comparisons)
dev.off()

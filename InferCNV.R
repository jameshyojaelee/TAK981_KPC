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
setwd("D:/TAK981_KPC/control_vs_TAK981")
load('cellchat_new.RData')

#InferCNV
kpc_matrix <- GetAssayData(scKPC, slot="counts")
annotations <- as.data.frame(scKPC$SingleR.label)
write.table(annotations, file="cellAnnotations.txt", sep="\t", col.names = FALSE,quote = FALSE)

kpc_cnv <- CreateInfercnvObject(raw_counts_matrix=kpc_matrix,
                                annotations_file="cellAnnotations.txt",
                                delim="\t",
                                gene_order_file="vM28.annotation.txt",
                                ref_group_names=c("Endothelial"))
# memory.limit(56000)
kpc_cnv = infercnv::run(kpc_cnv,
                        cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                        out_dir="InferCNV",  # dir is auto-created for storing outputs
                        cluster_by_groups=T,   # cluster
                        denoise=T,
                        num_threads = 64,
                        HMM=T)

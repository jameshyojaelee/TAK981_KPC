# scRNA-seq functions

get_group_proportions <- function (seuratobj, group.by = "SingleR.label") {
  if (group.by == "SingleR.label") {
    seuratobj[["SingleR.label"]] <- seuratobj@meta.data$SingleR.label
  }
  # gets the total number of cells within each group
  total_populations <- seuratobj@meta.data %>% group_by(treatment) %>% summarize (total.pop = n())
  # gets the proportion of cells for each cell type within a group by dividing by the total
  count_populations <- seuratobj@meta.data %>% group_by_at(vars(group.by, "treatment")) %>% summarize (n = n())
  count_populations <- left_join(count_populations, total_populations, by = "treatment")
  count_populations <- count_populations %>% mutate (proportion = n/total.pop)
  count_populations
}

ggtheme <- function (base_size = 11, base_family = "") {
  theme_bw() %+replace% 
    theme(
      panel.grid.major  = element_line(color = "#808080"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "#808080", fill = NA),
      axis.title.x = element_text(colour = "#808080", margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(colour = "#808080", angle=90, margin = margin(t = 0, r = 20, b = 0, l = 0)),
      axis.line = element_line(color = "#808080"),
      axis.ticks = element_line(color = "#808080"),
      axis.text = element_text(color = "black")
    )
}


plot_group_proportions <- function (seuratobj, graph.type) {
  # get the proportions using get_group_proportions()
  count_populations <- get_group_proportions(seuratobj)
  # plot the proportions
  if (graph.type == "dodge"){
    ggplot(count_populations, aes (x = SingleR.label, y = proportion))+
      geom_bar (aes(fill = treatment), stat = "identity", position = "dodge") +ggtheme()+
      theme(axis.text.x = element_text(angle=45, hjust = 1),
            axis.title.x = element_blank()) + 
      scale_fill_manual("legend", values = c("control" = "#818181",  "TAK981" = "#D30000"))
  } else if (graph.type == "stacked") {
    ggplot(count_populations, aes (x = treatment, y = proportion))+
      geom_bar (aes(fill = SingleR.label), stat = "identity", position = "fill") + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
      scale_fill_discrete(name = "Cell type") + theme_minimal() +
      theme(legend.text = element_text(size=12),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1))
  }
  else
    print("invalid graph type")
}

#Plot heamap proportions
plot_heatmap_proportions <- function (seuratobj, graph.type = "by.cell") {
  count_populations <- get_group_proportions(seuratobj)
  reformatted <- count_populations %>% reshape2::dcast(formula = treatment~SingleR.label, value.var = "proportion")
  heatmap_matrix <- as.matrix(reformatted[,-1])
  rownames(heatmap_matrix) <- reformatted[,1]
  if (graph.type == "by.cell") {
    pheatmap::pheatmap(heatmap_matrix, scale = "column",cluster_rows = F)
  } else if (graph.type == "by.sample") {
    pheatmap::pheatmap(cor(t(heatmap_matrix)))
  }
  else
    print("invalid graph type")
}





#fGSEA
# formats the ranked list for the fgsea() function
prepare_ranked_list <- function(ranked_list) { 
  # if duplicate gene names present, average the values
  if( sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(.~Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$avg_log2FC, decreasing = T),]
  }
  # omit rows with NA values
  ranked_list <- na.omit(ranked_list)
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}

plot_enrichment <- function (geneset, pathway, ranked_list) {
  plotEnrichment(geneset[[pathway]], ranked_list)+labs (title = pathway)
}


waterfall_plot <- function (fsgea_results, graph_title) {
  fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot(aes(reorder(short_name,NES), NES)) +
    geom_bar(stat= "identity", aes(fill = padj<0.05))+
    coord_flip()+
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
    theme(axis.text.y = element_text(size = 7), 
          plot.title = element_text(hjust = 1))
}




# Violinplot p-value
vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(pbmc, features = signature,
            pt.size = 0.1, 
            group.by = "Response", 
            y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  purrr::map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 14, height = 8)
}


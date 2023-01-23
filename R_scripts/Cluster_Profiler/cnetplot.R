# cnetplot example
?cnetplot

cnetplot(
  enrich,
  showCategory = 5,
  foldChange = gene_list,
  layout = "kk",
  colorEdge = T,
  circular = T,
  # node_label = "category",
  cex_category = 0.5,
  cex_gene = 0.5,
  cex_label_category = 0.5,
  cex_label_gene = 0.25,
  color_category = "palegreen4",
  color_gene = "#B3B3B3",
  shadowtext = "category"
  # color.params = list(foldChange = gene_list, edge = T, category = "#E5C494", gene = "#B3B3B3"),
  # cex.params = list(category_node = 1, gene_node = 1, category_label = 1, gene_label = 0.25),
  #hilight.params = list(category = ", alpha_hilight = 1, alpha_no_hilight = 0.3)
)

cnetplot(
  enrich,
  showCategory = 5,
  foldChange = gene_list,
  layout = "circle",
  colorEdge = T,
  circular = T,
  # node_label = "category",
  cex_category = 0.5,
  cex_gene = 0.5,
  cex_label_category = 0.5,
  cex_label_gene = 0.25,
  color_category = "palegreen4",
  color_gene = "#B3B3B3",
  shadowtext = "category"
  # color.params = list(foldChange = gene_list, edge = T, category = "#E5C494", gene = "#B3B3B3"),
  # cex.params = list(category_node = 1, gene_node = 1, category_label = 1, gene_label = 0.25),
  #hilight.params = list(category = ", alpha_hilight = 1, alpha_no_hilight = 0.3)
)




#' Heatmap of average expression values per group
#'
#' @param object Seurat object
#' @param features Which features to plot
#' @param assay Which assay to get expression values from
#' @param slot Which slot to pull the assay data from
#' @param group.by Which metadata slot to group cells by
#' @param max_zed Cutoff for maximum z-scores (absolute value)
#' @param y_text_size Size of text on y axis (features)
#' @param scale_rows Whether to scale the expression per row
#' @param plot_dendro Whether to plot accompanying x-axis and y-axis dendrograms
#' @param dot.scale Value by which to scale the dots in DotPlot
#' @param cluster_cols Whether to cluster the columns
#' @details Plots the average expression value of the specified genes/features
#' for each group in the Seurat object, optionally using hierarchical
#' clustering for prettier visualization
#'
#' @return Returns a patchwork object of a heatmap and dendrogram if specified
#'
#' @importFrom rlang %||%
#'
DoGenesetHeatmap <- function(seurat_object,
                             assay = "AUC",
                             slot = "counts",
                             group.by = "seurat_clusters",
                             diff_exp_results,
                             max_zed = 3,
                             plot_rownames = FALSE,
                             label_features = NULL,
                             y_text_size = 10,
                             scale_rows = TRUE,
                             p_val_choice = 0.01,
                             logFC = 0.4,
                             plot_all = FALSE,
                             plot_dendro = FALSE,
                             auc_choice = 0.6){

  if(all(c("logFC", "auc", "padj") %in% colnames(diff_exp_results))){
    message("diff_exp_results in Presto format")
    diff_exp_results <- diff_exp_results %>%
      dplyr::mutate(is_significant = ifelse(auc > auc_choice, "*", "")) %>%
      dplyr::mutate(cluster = group)
    sig_genes <- dplyr::filter(diff_exp_results, is_significant == "*") %>% pull(feature) %>% unique()
  } else if (all(c("p_val_adj", "avg_log2FC")) %in% colnames(diff_exp_results)){
    message("diff_exp_results in Seurat format")
    diff_exp_results <- diff_exp_results %>%
      dplyr::mutate(is_significant = ifelse((p_val_adj < p_val_choice) & (avg_log2FC > avg_log2FC_choice), "*", "")) %>%
      dplyr::mutate(feature = gene)
    sig_genes <- dplyr::filter(diff_exp_results, is_significant == "*") %>% pull(feature) %>% unique()
  } else {
    stop("diff_exp_results must be in Seurat or Presto format")
  }

  my_data <- TrueAverageExpression(seurat_object, slot = slot, assay = assay, group.by = group.by)
  if(scale_rows){
    my_data <- as.data.frame(t(scale(t(my_data), center = TRUE, scale = TRUE)))
  }

  if(!plot_all){
    my_data <- my_data[sig_genes,]
  }

  my_data[my_data > max_zed] <- max_zed
  my_data[my_data < -max_zed] <- -max_zed
  hclust_results <- hclust(dist(my_data))
  feature_order <- hclust_results$labels[hclust_results$order]
  my_data_long <- my_data %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "feature") %>%
    tidyr::pivot_longer(-feature, values_to = "Z_sc", names_to = "cluster") %>%
    dplyr::mutate(feature = factor(feature, levels = feature_order))


  my_data_long <- left_join(my_data_long, diff_exp_results, by = c("cluster", "feature")) %>%
    dplyr::mutate(feature = factor(feature, levels = feature_order)) %>%
    dplyr::mutate(cluster = factor(cluster, levels = colnames(my_data)))

  g_heatmap <- ggplot2::ggplot(my_data_long, aes(x = cluster, y = feature, fill = Z_sc))+
    ggplot2::geom_tile()+
    scale_fill_viridis_c()+
    ylab(NULL)+
    ggplot2::geom_text(aes(label = is_significant), vjust = 0.75)+
    scale_y_discrete(position = "left", expand = c(0.01, 0.01))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = y_text_size),
          plot.margin = margin(l = 0.1),
          panel.background = element_rect(fill = "white"),
          axis.text.y = element_text(size = y_text_size),
          legend.text = element_text(size = y_text_size * 0.75),
          legend.title = element_text(size = y_text_size * 0.75),
          legend.key.width = unit(5, "pt"),
          legend.key.height = unit(10, "pt"))+
    ylab(NULL)+
    xlab(NULL)
  g_dendro <- ggdendro::ggdendrogram(data = hclust_results)+
    scale_y_reverse(expand = c(0.01,0.01))+
    coord_flip()+
    theme(axis.text.y = element_blank(), axis.text.x = element_blank(), plot.margin = margin(r = 0))+
    scale_x_discrete(expand = c(0.01,0.01))
  if(plot_dendro){
    g_dendro + g_heatmap + patchwork::plot_layout(widths = c(1,2), guides = "collect")
  } else {
    g_heatmap
  }
}

#' Heatmap of average expression values per group with stars
#'
#' @param object Seurat object
#' @details TBD!!!
#'
#' @return Returns a patchwork object of a heatmap and dendrogram if specified
#'
#' @importFrom rlang %||%
#'
#' @export
DoStarHeatmap <- function(object,
                          diff_exp_results = NULL,
                          assay = "ADT",
                          slot = "data",
                          group.by = "seurat_clusters",
                          max_zed = 3,
                          plot_rownames = FALSE,
                          label_features = TRUE,
                          y_text_size = 10,
                          scale_rows = TRUE,
                          p_val_choice = 0.01,
                          logFC = 0.4,
                          plot_all = FALSE,
                          plot_dendro = FALSE,
                          subset_features = NULL,
                          auc_choice = 0.6,
                          cluster_cols = TRUE){

  if(!is.null(diff_exp_results)){
    if(all(c("logFC", "auc", "padj") %in% colnames(diff_exp_results))){
      message("diff_exp_results in Presto format")
      diff_exp_results <- diff_exp_results %>%
        dplyr::mutate(is_significant = ifelse(auc > auc_choice, "*", "")) %>%
        dplyr::mutate(cluster = group)
      sig_genes <- dplyr::filter(diff_exp_results, is_significant == "*") %>% pull(feature) %>% unique()
    } else if (all(c("p_val_adj", "avg_log2FC") %in% colnames(diff_exp_results))){
      message("diff_exp_results in Seurat format")
      diff_exp_results <- diff_exp_results %>%
        dplyr::mutate(is_significant = ifelse((p_val_adj < p_val_choice) & (avg_log2FC > logFC), "*", "")) %>%
        dplyr::mutate(feature = gene)
      sig_genes <- dplyr::filter(diff_exp_results, is_significant == "*") %>% pull(feature) %>% unique()
    } else {
      stop("diff_exp_results must be in Seurat or Presto format")
    }
  }

  my_data <- AverageExpression(object, slot = slot, assay = assay, group.by = group.by)[[assay]]
  if(!is.null(subset_features)){
    my_data <- my_data[subset_features,]
  }
  if(scale_rows){
    my_data <- as.data.frame(t(scale(t(my_data), center = TRUE, scale = TRUE)))
  }

  if(!plot_all){
    my_data <- my_data[rownames(my_data) %in% sig_genes,]
  }

  my_data[my_data > max_zed] <- max_zed
  my_data[my_data < -max_zed] <- -max_zed
  hclust_results <- hclust(dist(my_data))
  feature_order <- hclust_results$labels[hclust_results$order]
  hclust_results_2 <- hclust(dist(t(my_data)))
  cluster_order <- hclust_results_2$labels[hclust_results_2$order]


  my_data_long <- my_data %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "feature") %>%
    pivot_longer(-feature, values_to = "Z_sc", names_to = "cluster") %>%
    dplyr::mutate(feature = factor(feature, levels = feature_order))

  my_data_long <- left_join(my_data_long, diff_exp_results, by = c("cluster", "feature")) %>%
    dplyr::mutate(feature = factor(feature, levels = feature_order)) %>%
    dplyr::mutate(cluster = factor(cluster, levels = colnames(my_data)))

  if(cluster_cols){
    my_data_long <- my_data_long %>%
      dplyr::mutate(cluster = factor(cluster, levels = cluster_order))
  }

  g_heatmap <- ggplot(my_data_long, aes(x = cluster, y = feature, fill = Z_sc))+
    geom_tile(color = "black")+
    scale_fill_viridis_c()+
    ylab(NULL)+
    geom_text(aes(label = is_significant), vjust = 0.75)+
    scale_x_discrete(position = "bottom", expand = c(0,0))+
    scale_y_discrete(position = "right", expand = c(0.01, 0.01))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = y_text_size),
          plot.margin = margin(t = 0),
          panel.background = element_rect(fill = "white"),
          axis.text.y = element_text(size = y_text_size),
          legend.text = element_text(size = y_text_size * 0.75),
          legend.title = element_text(size = y_text_size * 0.75),
          legend.key.width = unit(5, "pt"),
          legend.key.height = unit(10, "pt"))+
    xlab(NULL)
  g_dendro_rows <- ggdendro::ggdendrogram(data = hclust_results)+
    scale_y_reverse(expand = c(0,0))+
    scale_x_discrete(expand = c(0.01,0.01))+
    coord_flip()+
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = margin(r = 0, t = 0, b = 0),
          panel.background = element_rect("white"))
  g_dendro_cols <- ggdendro::ggdendrogram(data = hclust_results_2)+
    scale_y_discrete(expand = c(0,0))+
    scale_x_discrete(expand = c(0.01,0.01))+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, l = 0, b = 0))

  if(plot_dendro){
    patchwork::plot_spacer() +g_dendro_cols+ g_dendro_rows + g_heatmap +
      patchwork::plot_layout(widths = c(1,8), heights = c(1,8),
                             ncol = 2,
                             nrow = 2)
  } else {
    g_heatmap
  }
}



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
#' @param cluster_rows Whether to cluster the rows
#' @details Plots the average expression value of the specified genes/features
#' for each group in the Seurat object, optionally using hierarchical
#' clustering for prettier visualization
#'
#' @return Returns a patchwork object of a heatmap and dendrogram if specified
#'
#' @importFrom rlang %||%
#'
#' @export
DoClusteredHeatmap <- function(object,
                               features = NULL,
                               assay = "RNA",
                               slot = "scale.data",
                               group.by = "seurat_clusters",
                               max_zed = 3,
                               y_text_size = 10,
                               scale_rows = FALSE,
                               plot_dendro = FALSE,
                               dot.scale = 5,
                               cluster_cols = TRUE,
                               cluster_rows = TRUE){

  if(slot == "scale.data" & !all(features %in% rownames(object[[assay]]@scale.data))){
    object <- ScaleData(object, assay = assay, features = features)
  }
  my_data <- AverageExpression(object, slot = slot, assay = assay, group.by = group.by)[[assay]]
  if(!all(features %in% rownames(my_data))){
    stop("verify that all features are in the appropriate slot and assay (may need to ScaleData manually)")
  }
  my_data <- my_data[features,]

  if(scale_rows){
    my_data <- as.data.frame(t(scale(t(my_data), center = TRUE, scale = TRUE)))
  }

  my_data[my_data > max_zed] <- max_zed
  my_data[my_data < -max_zed] <- -max_zed

  my_data_long <- my_data %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "feature") %>%
    pivot_longer(-feature, values_to = "Z_sc", names_to = "cluster") %>%
    dplyr::mutate(feature = factor(feature, levels = features)) %>%
    dplyr::mutate(cluster = factor(cluster, levels = colnames(my_data)))

  if(cluster_cols){
    hclust_results_2 <- hclust(dist(t(my_data)))
    cluster_order <- hclust_results_2$labels[hclust_results_2$order]
    my_data_long <- my_data_long %>%
      dplyr::mutate(cluster = factor(cluster, levels = cluster_order))
  }
  if(cluster_rows){
    hclust_results <- hclust(dist(my_data))
    feature_order <- hclust_results$labels[hclust_results$order]
    my_data_long <- my_data_long %>%
      dplyr::mutate(feature = factor(feature, levels = feature_order))
  }

  label_size <- y_text_size
  heatmap_theme <- theme(axis.text.x = element_text(size = label_size),
                         axis.text.y = element_text(size = label_size),
                         legend.text = element_text(size = label_size * .75),
                         legend.title = element_text(size = label_size * .75),
                         legend.key.width = unit(5, "pt"),
                         legend.key.height = unit(10, "pt"))
  base_heatmap <- ggplot(my_data_long, aes(x = cluster, y = feature, fill = Z_sc))+
    geom_tile()+
    scale_x_discrete(position = "bottom")+
    RotatedAxis()+
    ylab(NULL)+
    xlab(NULL)+
    scale_fill_viridis_c(option = "B", direction = 1)+
    heatmap_theme
  if(cluster_rows){
    g_dendro_rows <- ggdendro::ggdendrogram(data = hclust_results)+
      scale_y_reverse(expand = c(0,0))+
      scale_x_discrete(expand = c(0.01,0.01))+
      coord_flip()+
      theme(axis.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = margin(r = 0, t = 0, b = 0),
            panel.background = element_rect("white"))
  }
  if(cluster_cols){
    g_dendro_cols <- ggdendro::ggdendrogram(data = hclust_results_2)+
      scale_y_discrete(expand = c(0,0))+
      scale_x_discrete(expand = c(0.01,0.01))+
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = margin(t = 0, r = 0, l = 0, b = 0))
  }

  if(cluster_rows & cluster_cols & plot_dendro){
    patchwork::plot_spacer() +g_dendro_cols+ g_dendro_rows + base_heatmap +
      patchwork::plot_layout(widths = c(1,8), heights = c(1,8),
                             ncol = 2,
                             nrow = 2)
  } else {
    base_heatmap
  }
}




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
#' @param cluster_rows Whether to cluster the rows
#' @details Plots the average expression value of the specified genes/features
#' for each group in the Seurat object, optionally using hierarchical
#' clustering for prettier visualization
#'
#' @return Returns a patchwork object of a heatmap and dendrogram if specified
#'
#' @importFrom rlang %||%
#'
#' @export
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

  my_data <- AverageExpression(seurat_object, slot = slot, assay = assay, group.by = group.by)[[assay]]
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
    pivot_longer(-feature, values_to = "Z_sc", names_to = "cluster") %>%
    dplyr::mutate(feature = factor(feature, levels = feature_order))


  my_data_long <- left_join(my_data_long, diff_exp_results, by = c("cluster", "feature")) %>%
    dplyr::mutate(feature = factor(feature, levels = feature_order)) %>%
    dplyr::mutate(cluster = factor(cluster, levels = colnames(my_data)))

  g_heatmap <- ggplot(my_data_long, aes(x = cluster, y = feature, fill = Z_sc))+
    geom_tile()+
    scale_fill_viridis_c()+
    ylab(NULL)+
    geom_text(aes(label = is_significant), vjust = 0.75)+
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



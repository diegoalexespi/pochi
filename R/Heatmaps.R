#' Heatmap of average expression values per group with stars
#'
#' @param object Seurat object
#' @param viridis_option The viridis color scale to use
#' @param viridis_direction The direction of the viridis color scale
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
                          logFC_choice = 0.4,
                          auc_choice = NULL,
                          plot_all = FALSE,
                          plot_dendro = FALSE,
                          subset_features = NULL,
                          cluster_cols = TRUE,
                          cluster_rows = TRUE,
                          viridis_option = "C",
                          viridis_direction = 1){

  Idents(object) <- group.by
  og.levels <- levels(Idents(object))

  if(!is.null(diff_exp_results)){
    if(all(c("logFC", "auc", "padj") %in% colnames(diff_exp_results))){
      message("diff_exp_results in Presto format")
      if(!is.null(auc_choice)){
        message("using auc_choice and ignoring logFC_choice and p_val_choice")
        diff_exp_results <- diff_exp_results %>%
          dplyr::mutate(is_significant = ifelse(auc > auc_choice, "*", "")) %>%
          dplyr::mutate(cluster = group)
        sig_genes <- dplyr::filter(diff_exp_results, is_significant == "*") %>% dplyr::pull(feature) %>% unique()
      } else {
        diff_exp_results <- diff_exp_results %>%
          dplyr::mutate(is_significant = ifelse((padj < p_val_choice) & (logFC > logFC_choice), "*", "")) %>%
          dplyr::mutate(cluster = group)
        sig_genes <- dplyr::filter(diff_exp_results, is_significant == "*") %>% dplyr::pull(feature) %>% unique()
      }

    } else if (all(c("p_val_adj", "avg_log2FC") %in% colnames(diff_exp_results))){
      message("diff_exp_results in Seurat format")
      diff_exp_results <- diff_exp_results %>%
        dplyr::mutate(is_significant = ifelse((p_val_adj < p_val_choice) & (avg_log2FC > logFC_choice), "*", "")) %>%
        dplyr::mutate(feature = gene)
      sig_genes <- dplyr::filter(diff_exp_results, is_significant == "*") %>% dplyr::pull(feature) %>% unique()
    } else {
      stop("diff_exp_results must be in Seurat or Presto format")
    }
  } else {
    stop("must input a diff_exp_results (in Seurat or Presto format)")
  }

  my_data <- TrueAverageExpression(object, slot = slot, assay = assay, group.by = group.by)
  if(!is.null(subset_features)){
    my_data <- my_data[subset_features,]
  }
  if(!plot_all){
    my_data <- my_data[rownames(my_data) %in% sig_genes,]
  }

  if(scale_rows){
    my_data <- as.data.frame(t(scale(t(my_data), center = TRUE, scale = TRUE)))
  }

  my_data[my_data > max_zed] <- max_zed
  my_data[my_data < -max_zed] <- -max_zed

  my_data_long <- my_data %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "feature") %>%
    tidyr::pivot_longer(-feature, values_to = "Z_sc", names_to = "cluster") %>%
    dplyr::mutate(cluster = factor(cluster, levels = colnames(my_data))) %>%
    dplyr::left_join(diff_exp_results, by = c("cluster", "feature"))

  if(cluster_cols){
    hclust_results_2 <- hclust(dist(t(my_data)))
    cluster_order <- hclust_results_2$labels[hclust_results_2$order]
    my_data_long <- my_data_long %>%
      dplyr::mutate(cluster = factor(cluster, levels = cluster_order))
  } else {
    message("i")
    my_data_long <- my_data_long %>%
      dplyr::mutate(cluster = factor(cluster, levels = og.levels))
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
  base_heatmap <- ggplot2::ggplot(my_data_long, aes(x = cluster, y = feature, fill = Z_sc))+
    ggplot2::geom_tile()+
    ggplot2::geom_text(aes(label = is_significant), vjust = 0.75)+
    scale_x_discrete(position = "bottom")+
    RotatedAxis()+
    ylab(NULL)+
    xlab(NULL)+
    scale_fill_viridis_c(option = viridis_option, direction = viridis_direction)+
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
#' @param viridis_option The viridis color scale to use
#' @param viridis_direction The direction of the viridis color scale
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
                               slot = "data",
                               group.by = "seurat_clusters",
                               max_zed = 3,
                               y_text_size = 10,
                               scale_rows = TRUE,
                               plot_dendro = FALSE,
                               dot.scale = 5,
                               cluster_cols = TRUE,
                               cluster_rows = TRUE,
                               viridis_option = "B",
                               viridis_direction = 1){

  my_data <- TrueAverageExpression(object, slot = slot, assay = assay, group.by = group.by)
  my_data <- my_data[features,]

  if(scale_rows){
    my_data <- as.data.frame(t(scale(t(my_data), center = TRUE, scale = TRUE)))
  }

  my_data[my_data > max_zed] <- max_zed
  my_data[my_data < -max_zed] <- -max_zed

  my_data_long <- my_data %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "feature") %>%
    tidyr::pivot_longer(-feature, values_to = "Z_sc", names_to = "cluster") %>%
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
  base_heatmap <- ggplot2::ggplot(my_data_long, aes(x = cluster, y = feature, fill = Z_sc))+
    ggplot2::geom_tile()+
    scale_x_discrete(position = "bottom")+
    RotatedAxis()+
    ylab(NULL)+
    xlab(NULL)+
    scale_fill_viridis_c(option = viridis_option, direction = viridis_direction)+
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




#' Average feature expression across clustered samples in a Seurat object using fast sparse matrix methods
#'
#' @param object Seurat object
#' @param group.by Ident with sample clustering information (default is
#' seurat_clusters). Can also input a vector with multiple idents specified and
#' TrueAverageExpression will return specified idents concatenated into new
#' a new output ident with group.by.delim as the separator.
#' @param group.by.delim The delimiter used if group.by is a vector of
#' length > 1.
#' @param assay Assay to average (default is the active assay)
#' @param slot Slot to average (default is counts)
#' @param verbose Boolean or integer, show progress bar (default is TRUE)
#'
#' @references credit to GitHub user @zdebruine on Seurat issues!
#'
#' @return dense matrix of cluster averages
#'
#' @export
TrueAverageExpression <- function(object,
                                  group.by = "seurat_clusters",
                                  group.by.delim = "_",
                                  assay = "RNA",
                                  slot = "data",
                                  verbose = TRUE) {

  my_data <- GetAssayData(object, assay = assay, slot = slot)

  idents <- object@meta.data[,group.by,drop=FALSE] %>%
    tidyr::unite("ident_vector", group.by) %>%
    dplyr::pull(ident_vector)

  # loop through all idents, averaging them in data
  ident.names <- unique(idents)
  if (verbose > 0) pb <- txtProgressBar(char = "=", style = 1, max = length(ident.names), width = 50)
  m <- list()
  for (i in 1:length(ident.names)) {
    m[[i]] <- Matrix::rowMeans(my_data[, which(idents == ident.names[i])])
    if (verbose > 0) setTxtProgressBar(pb = pb, value = i)
  }
  result <- do.call(cbind, m)
  colnames(result) <- ident.names
  return(result)
}



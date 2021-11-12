#' Seurat FeaturePlot enhancement
#'
#' @param object Seurat object
#' @param feature Which feature to plot (singular)
#' @param min.value Which threshold value to set as not-expressed/grey in the
#' FeaturePlot
#' @param assay Which assay the feature is located in
#' @param viridis_choice Which viridis palette to choose
#' @param ... Other arguments to Seurat::FeaturePlot
#' @details Identical to FeaturePlot with the exception that 0-values are set to
#' grey, which replicates the Monocle-type UMAP plots.
#'
#' @return Returns a FeaturePlot
#'
#' @importFrom rlang %||%
#'
#' @export
FeaturePlotPlus <- function(object,
                            feature,
                            min.value = 0,
                            assay = "RNA",
                            viridis_choice = "D",
                            ...){
  DefaultAssay(object) <- assay
  temp_plot <- FeaturePlot(object, features = feature, ...)
  temp_plot$data[,feature] <- ifelse(temp_plot$data[,feature] <= min.value, NA, temp_plot$data[,feature])
  temp_plot +
    viridis::scale_color_viridis(option = viridis_choice, na.value = "grey")
}


#' Seurat DimPlot enhancement
#'
#' @param object Seurat object
#' @param graph.choice Which Graph to utilize for plotting with the DimPlot
#' @param reduction Which reduction to use in the DimPlot
#' @param line.thickness Thickness of lines connecting individual points in the
#' DimPlot based on the graph.choice Graph
#' @param line.color Color of the lines connecting individual points
#' @param raster Raseterization boolean
#' @param ... Other arguments to Seurat::DimPlot
#' @details Identical to DimPlot with the exception that edges can be drawn on
#' the DimPlot between cells connected in the specified Graph.
#'
#' @return Returns a DimPlot
#'
#' @importFrom rlang %||%
#'
#' @export
DimPlotEdges <- function(object,
                         graph.choice = "RNA_nn",
                         reduction = "umap",
                         line.thickness = 0.05,
                         line.color = "grey",
                         raster = FALSE,
                         ...){

  my_graph <- as(Graphs(object, graph.choice), "dgCMatrix")
  my_graph_coordinates <- data.frame(Matrix::summary(my_graph))
  my_graph_coordinates$first_cell <- my_graph@Dimnames[[1]][my_graph_coordinates$i]
  my_graph_coordinates$second_cell <- my_graph@Dimnames[[1]][my_graph_coordinates$j]

  my_dimreduc_coordinates <- Embeddings(object, reduction) %>%
    as.data.frame() %>%
    set_colnames(value = c("coord_x", "coord_y")) %>%
    rownames_to_column(var = "cellID")

  my_graph_coordinates <- my_graph_coordinates %>%
    dplyr::left_join(my_dimreduc_coordinates, by = c(first_cell = "cellID")) %>%
    dplyr::rename("coord_x.1" = coord_x, "coord_y.1" = coord_y) %>%
    dplyr::left_join(my_dimreduc_coordinates, by = c(second_cell = "cellID")) %>%
    dplyr::rename("coord_x.2" = coord_x, "coord_y.2" = coord_y) %>%
    dplyr::mutate(keep_edge = ifelse(first_cell == second_cell, "no", "yes")) %>%
    dplyr::filter(keep_edge == "yes") %>%
    dplyr::filter(!duplicated(paste0(pmax(first_cell, second_cell), pmin(first_cell, second_cell))))

  temp_dimplot <- DimPlot(object, reduction = reduction, raster = raster, ...)+
    ggrastr::rasterize(geom_segment(data = my_graph_coordinates, aes(x = coord_x.1, y = coord_y.1, xend = coord_x.2, yend = coord_y.2), inherit.aes = FALSE, size = line.thickness, color = line.color))

  num_layers <- length(temp_dimplot$layers)
  temp_dimplot$layers <- temp_dimplot$layers[c(num_layers, 1:(num_layers-1))]
  temp_dimplot
}

#' Seurat DimPlot enhancement II
#'
#' @param object Seurat object
#' @param reduction Which reduction to use in the DimPlot
#' @param shape.choice Which shape to use for the points int he DimPlot
#' @param pt.size Size of DimPlot points
#' @param st.size Size of DimPlot stencil around points
#' @param ... Other arguments to Seurat::DimPlot
#' @details Identical to DimPlot with the exception that points are drawn as
#' fill-able circles with a black border.
#'
#' @return Returns a DimPlot
#'
#' @importFrom rlang %||%
#'
#' @export
DimPlotPlus <- function(object,
                        reduction = "umap",
                        shape.choice = 21,
                        pt.size = 1.75,
                        st.size = 0.25,
                        ...){
  temp_dimplot <- DimPlot(object, reduction = reduction, ...)
  if(shape.choice == 21){
    temp_dimplot$layers[[1]]$mapping$fill <- temp_dimplot$layers[[1]]$mapping$colour
    temp_dimplot$layers[[1]]$mapping$colour <- NULL
    temp_dimplot$layers[[1]]$aes_params$shape <- 21
    temp_dimplot$layers[[1]]$aes_params$colour <- "black"
    temp_dimplot$layers[[1]]$aes_params$size <- pt.size
    temp_dimplot$layers[[1]]$aes_params$stroke <- st.size
  }
  temp_dimplot
}


#' Seurat DotPlot enhancement with hierarchical clustering
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
#' @param cluster_cols Whether to cluster the columns as well as the rows
#' @details Identical to DotPlot with the exception that the genes/features are
#' clustered using hierarchical clustering for prettier visualization
#'
#' @return Returns a patchwork object of a DotPlot and dendrogram if specified
#'
#' @importFrom rlang %||%
#'
#' @export
DoClusteredDotPlot <- function(object,
                               features = NULL,
                               assay = "RNA",
                               slot = "data",
                               group.by = "seurat_clusters",
                               max_zed = 3,
                               y_text_size = 10,
                               scale_rows = TRUE,
                               plot_dendro = FALSE,
                               auc_choice = 0.6,
                               dot.scale = 5,
                               cluster_cols = TRUE){



  # if(!is.null(diff_exp_results)){
  #   if(all(c("logFC", "auc", "padj") %in% colnames(diff_exp_results))){
  #     message("diff_exp_results in Presto format")
  #     diff_exp_results <- diff_exp_results %>%
  #       dplyr::mutate(is_significant = ifelse(auc > auc_choice, "*", "")) %>%
  #       dplyr::mutate(cluster = group)
  #     sig_genes <- dplyr::filter(diff_exp_results, is_significant == "*") %>% pull(feature) %>% unique()
  #   } else if (all(c("p_val_adj", "avg_log2FC") %in% colnames(diff_exp_results))){
  #     message("diff_exp_results in Seurat format")
  #     diff_exp_results <- diff_exp_results %>%
  #       dplyr::mutate(is_significant = ifelse((p_val_adj < p_val_choice) & (avg_log2FC > logFC), "*", "")) %>%
  #       dplyr::mutate(feature = gene)
  #     sig_genes <- dplyr::filter(diff_exp_results, is_significant == "*") %>% pull(feature) %>% unique()
  #   } else {
  #     stop("diff_exp_results must be in Seurat or Presto format")
  #   }
  #   if(!is.null(subset_features)){
  #     my_data <- my_data[subset_features,]
  #   }
  #   if(!plot_all){
  #     my_data <- my_data[rownames(my_data) %in% sig_genes,]
  #   }
  my_data <- AverageExpression(object, slot = slot, assay = assay, group.by = group.by)[[assay]]
  my_data <- my_data[features,]
  if(scale_rows){
    my_data <- as.data.frame(t(scale(t(my_data), center = TRUE, scale = TRUE)))
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
    dplyr::mutate(feature = factor(feature, levels = feature_order)) %>%
    dplyr::mutate(cluster = factor(cluster, levels = colnames(my_data)))

  # my_data_long <- left_join(my_data_long, diff_exp_results, by = c("cluster", "feature")) %>%
  #   dplyr::mutate(feature = factor(feature, levels = feature_order))

  if(cluster_cols){
    my_data_long <- my_data_long %>%
      dplyr::mutate(cluster = factor(cluster, levels = cluster_order))
  }


  label_size <- y_text_size
  dotplot_theme <- theme(axis.text.x = element_text(size = label_size),
                         axis.text.y = element_text(size = label_size),
                         panel.background = element_rect(color = "black", size = 1.5),
                         legend.text = element_text(size = label_size * .75),
                         legend.title = element_text(size = label_size * .75),
                         legend.key.width = unit(5, "pt"),
                         legend.key.height = unit(10, "pt"))
  base_dotplot <- DotPlot(object,
                          features = feature_order,
                          group.by = group.by,
                          assay = assay,
                          dot.scale = dot.scale)+
    coord_flip()+
    scale_x_discrete(position = "top")+
    scale_y_discrete(limits = cluster_order)+
    RotatedAxis()+
    ylab(NULL)+
    xlab(NULL)+
    scale_color_viridis_c(option = "D", direction = -1)+
    dotplot_theme


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
    patchwork::plot_spacer() +g_dendro_cols+ g_dendro_rows + base_dotplot +
      patchwork::plot_layout(widths = c(1,8), heights = c(1,8),
                             ncol = 2,
                             nrow = 2)
  } else {
    base_dotplot
  }
}

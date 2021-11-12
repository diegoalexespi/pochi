
# just a wrapper around uwot to return the kNN (and SNN) that uwot constructs,
# which isn't done in Seurat
BuildAnnoyUMAP <- function(
  seurat_object,
  dims,
  use_raw = FALSE,
  genes_exclude = "",
  reduction = "pca",
  graph_key = "annoy",
  assay = DefaultAssay(seurat_object),
  slot = "scale.data",
  umap_method = 'uwot',
  n_neighbors = 30L,
  n_components = 2L,
  metric = 'cosine',
  n_epochs = NULL,
  prune_snn = 1/15,
  learning_rate = 1.0,
  min_dist = 0.3,
  spread = 1.0,
  set_op_mix_ratio = 1.0,
  local_connectivity = 1L,
  repulsion_strength = 1,
  negative_sample_rate = 5,
  a = NULL,
  b = NULL,
  uwot_sgd = FALSE,
  seed_use = 42,
  metric_kwds = NULL,
  angular_rp_forest = FALSE,
  reduction_key = 'UMAP_',
  verbose = TRUE,
  reduction_name = "umap",
  return_model = TRUE,
  save_graphs = TRUE,
  init_pos = "spectral",
  scale = FALSE) {
  if (!is.null(x = seed_use)) {
    set.seed(seed = seed_use)
  }

  if(use_raw){
    data_use <- t(as.matrix(Seurat::GetAssayData(seurat_object, slot = slot, assay = assay)))
    data_use <- data_use[,!(colnames(data_use) %in% genes_exclude)]
  } else {
    data_use <- Seurat::Embeddings(seurat_object[[reduction]])[, dims]
  }

  if(save_graphs) graph_vector <- c("nn", "fgraph") else graph_vector <- c()

  umap_results <- uwot::umap(
    X = data_use,
    init = init_pos,
    n_neighbors = as.integer(x = n_neighbors),
    n_components = as.integer(x = n_components),
    metric = metric,
    n_epochs = n_epochs,
    learning_rate = learning_rate,
    min_dist = min_dist,
    spread = spread,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    repulsion_strength = repulsion_strength,
    negative_sample_rate = negative_sample_rate,
    a = a,
    b = b,
    fast_sgd = uwot_sgd,
    verbose = verbose,
    ret_extra = graph_vector,
    scale = scale,
    ret_model = TRUE
  )

  if(save_graphs){
    nn_metric <- umap_results$nn[[metric]]
    nn_metric_vertices <- nrow(nn_metric$idx)
    nn_metric_sparse_dist <- Matrix::sparseMatrix(i = rep(1:nn_metric_vertices, each = n_neighbors),
                                                  j = as.vector(t(nn_metric$idx)),
                                                  x = as.vector(t(nn_metric$dist)))
    nn_metric_sparse_adj <- Matrix::sparseMatrix(i = rep(1:nn_metric_vertices, each = n_neighbors),
                                                 j = as.vector(t(nn_metric$idx)),
                                                 x = 1)
    nn_metric_sparse_snn <- Matrix::tcrossprod(nn_metric_sparse_adj)
    nn_metric_sparse_snn@x <- nn_metric_sparse_snn@x / (2 * n_neighbors - nn_metric_sparse_snn@x)
    nn_metric_sparse_snn@x[nn_metric_sparse_snn@x < prune_snn] <- 0
    nn_metric_sparse_snn <- Matrix::drop0(nn_metric_sparse_snn)
    nn_fuzzy <- umap_results$fgraph

    nn_list <- list(nn_metric_sparse_dist, nn_metric_sparse_adj, nn_metric_sparse_snn, nn_fuzzy)

    for(i in 1:length(nn_list)){
      colnames(nn_list[[i]]) <- colnames(seurat_object)
      rownames(nn_list[[i]]) <- colnames(seurat_object)
      nn_list[[i]] <- as.Graph(nn_list[[i]])
      DefaultAssay(nn_list[[i]]) <- assay
    }

    seurat_object[[paste0(graph_key, paste0("_nn_", metric))]] <- nn_list[[1]]
    seurat_object[[paste0(graph_key, "_nn_nn")]] <- nn_list[[2]]
    seurat_object[[paste0(graph_key, "_nn_snn")]] <- nn_list[[3]]
    seurat_object[[paste0(graph_key, "_fuzzy")]] <- nn_list[[4]]
  }

  umap_output <- umap_results$embedding
  colnames(x = umap_output) <- paste0(reduction_key, 1:ncol(x = umap_output))
  rownames(x = umap_output) <- rownames(data_use)
  umap_reduction <- Seurat::CreateDimReducObject(
    embeddings = umap_output,
    key = reduction_key,
    assay = assay,
    global = TRUE
  )
  #added for compatibility with @immunogenomics/symphony
  if(return_model){
    Misc(umap_reduction, slot = "model") <- umap_results
  }

  seurat_object[[reduction_name]] <- umap_reduction
  return(seurat_object)
}


#code adapted from @immunogenomics/singlecellmethods
BuildSNNfromKNN <- function(seurat_object,
                            knn_choice = "RNA_nn",
                            snn_name = paste0(knn_choice, "_snn"),
                            prune_snn = 1/15,
                            nn_k = NULL) {
  my_knn <- as(Graphs(seurat_object, knn_choice), "dgCMatrix")
  if(is.null(nn_k)){
    nn_k <- unique(Matrix::rowSums(my_knn))
    if(length(nn_k) != 1){
      stop("your kNN graph doesn't have equal NN for each node")
    }
  }
  message("building SNN from ", knn_choice)
  my_snn <- Matrix::tcrossprod(as(my_knn, "dgCMatrix"))
  my_snn@x <- my_snn@x / (2 * nn_k - my_snn@x)
  my_snn@x[my_snn@x < (prune_snn)] <- 0
  my_snn <- Matrix::drop0(my_snn)
  seurat_object[[snn_name]] <- as.Graph(my_snn)
  return(seurat_object)
}



#very naive version - may crash in some user cases
ToMonocle3 <- function(seurat_object,
                       scale_all = FALSE,
                       assay_to_keep = "RNA",
                       reduction_for_projection = "PCA",
                       UMAP_cluster_slot = NULL){


  # We prep the seurat object by creating gene loadings for ALL genes in the Seurat scale.data slot.
  # This is done to allow downstream monocle3 functions on gene_modules to work appropriately.
  message(paste0("Creating gene loadings using ", reduction_for_projection))
  projection_assay <- seurat_object[[reduction_for_projection]]@assay.used
  if(scale_all){
    message("Scaling all genes in assay used for computing reduction_for_projection")
    seurat_genes <- rownames(seurat_object[[projection_assay]])
    remaining_genes <- setdiff(seurat_genes, rownames(seurat_object[[projection_assay]]@scale.data))
    if(projection_assay == "SCT"){
      seurat_object <- Seurat::GetResidual(seurat_object, features = remaining_genes, assay = projection_assay, umi.assay = "RNA")
    } else {
      seurat_object <- Seurat::ScaleData(seurat_object, features = rownames(seurat_object[[projection_assay]]))
    }
  }
  seurat_object <- Seurat::ProjectDim(seurat_object, reduction = reduction_for_projection, assay = projection_assay)

  ##################

  message("Initializing CDS object")

  #Extract Seurat's log-transformed values
  expression_matrix <- Seurat::GetAssayData(seurat_object, assay = assay_to_keep, slot = "counts")
  #Extract Seurat meta_data
  meta_data <- seurat_object@meta.data
  #Extract gene names from Seurat object SCT slot to make CDS
  seurat_genes <- data.frame(gene_short_name = rownames(seurat_object[[assay_to_keep]]),
                             row.names = rownames(seurat_object[[assay_to_keep]]))
  new_cds <- monocle3::new_cell_data_set(expression_data = expression_matrix, cell_metadata = meta_data, gene_metadata = seurat_genes)

  ##################

  message("Making an SCE object from the Seurat object to facilitate transfer of information from SCE to CDS")
  sce <- as.SingleCellExperiment(seurat_object, assay = assay_to_keep)
  message("Loading in all Seurat reductions (PCA, HARMONY, UMAP, etc.) into CDS")
  SingleCellExperiment::reducedDims(new_cds) <- SingleCellExperiment::reducedDims(sce)
  message("Loading in specified Seurat assay into CDS")
  SummarizedExperiment::assays(new_cds) <- SummarizedExperiment::assays(sce)
  message("Loading in Seurat gene names into CDS")
  SummarizedExperiment::rowData(new_cds) <- SummarizedExperiment::rowData(sce)
  SummarizedExperiment::rowData(new_cds)$gene_short_name <-  row.names(new_cds)
  message("Loading in Seurat gene loadings into CDS")
  new_cds@preprocess_aux$gene_loadings <- seurat_object@reductions[[reduction_for_projection]]@feature.loadings.projected

  ##################

  message("Get user specified selected clusters (or active idents) from Seurat and load into CDS")
  if(is.null(UMAP_cluster_slot)){
    list_cluster <- Idents(seurat_object)
  } else {
    Idents(seurat_object) <- UMAP_cluster_slot
    list_cluster <- Idents(seurat_object)
  }
  new_cds@clusters[["UMAP"]]$clusters <- list_cluster
  #The next two commands are run in order to allow "order_cells" to be run in monocle3
  rownames(new_cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
  colnames(SingleCellExperiment::reducedDims(new_cds)[["UMAP"]]) <- NULL

  ##################

  message("Setting all cells as belonging to one partition (multiple partitions not supported yet)")
  recreate_partition <- c(rep(1, length(new_cds@colData@rownames)))
  names(recreate_partition) <- new_cds@colData@rownames
  recreate_partition <- as.factor(recreate_partition)
  new_cds@clusters[["UMAP"]]$partitions <- recreate_partition

  ##################
  message("Done")
  new_cds
}

MetaDataPlot <- function(seurat_object,
                         group.by,
                         split.by,
                         as_freq = TRUE,
                         label_count = FALSE,
                         text_size = 3,
                         label_angle = 45){

  seurat_metadata <- seurat_object@meta.data
  grouped_count <- seurat_metadata %>%
    dplyr::group_by(!!sym(split.by), !!sym(group.by)) %>%
    dplyr::tally() %>%
    dplyr::mutate(y = n)
  if(as_freq){
    grouped_count <- dplyr::mutate(grouped_count, y = y/sum(y))
  }
  g <- ggplot(grouped_count, aes(x = !!sym(split.by), y = y, fill = !!sym(group.by)))+
    geom_col(color = "black")+
    theme_cowplot()
  if(label_angle == 90){
    g <- g + theme(axis.text.x = element_text(angle  = label_angle, hjust = 1, vjust = 0.5))
  } else if (label_angle == 45){
    g <- g + theme(axis.text.x = element_text(angle  = label_angle, hjust = 1, vjust = 1))
  }

  if(as_freq){
    g <- g + ylab("Frequency")
  } else {
    g <- g + ylab("Count")
    if(label_count){
      g <- g + geom_text(aes(label = stat(y), group = !!sym(split.by)),
                         stat = 'summary', fun = sum, vjust = -1, size = text_size)
    }
  }
  g
}

# AbundancePlot <- function(seurat_object,
#                           group.by,
#                           split.by,
#                           replicate.by,
#                           perform.t.test = TRUE,
#                           ncols = 3,
#                           paired = FALSE,
#                           set_minimum_to_zero = FALSE){
#
#   seurat_metadata <- seurat_object@meta.data
#   seurat_metadata_filtered <- seurat_metadata %>%
#     dplyr::group_by(!!sym(replicate.by), !!sym(split.by), !!sym(group.by)) %>%
#     dplyr::tally() %>%
#     dplyr::mutate(Frequency = n/sum(n)) %>%
#     dplyr::ungroup() %>%
#     droplevels() %>%
#     tidyr::complete(nesting(!!sym(replicate.by), !!sym(split.by)), !!sym(group.by), fill = list(Frequency = 0, n = 0))
#
#   split.by.values <- unique(seurat_metadata_filtered[[split.by]])
#   if(length(split.by.values) != 2){
#     stop("split.by must denote a column of only 2 groups")
#   }
#
#   frequencies_1 <- seurat_metadata_filtered %>%
#     dplyr::filter(!!sym(split.by) == split.by.values[1]) %>%
#     dplyr::ungroup() %>%
#     dplyr::select(!!sym(replicate.by), !!sym(group.by), Frequency) %>%
#     tidyr::pivot_wider(names_from = !!sym(group.by), values_from = Frequency, values_fill = 0)
#
#   frequencies_2 <- seurat_metadata_filtered %>%
#     dplyr::filter(!!sym(split.by) == split.by.values[2]) %>%
#     dplyr::ungroup() %>%
#     dplyr::select(!!sym(replicate.by), !!sym(group.by), Frequency) %>%
#     tidyr::pivot_wider(names_from = !!sym(group.by), values_from = Frequency, values_fill = 0)
#
#   group.by.names <- sort(colnames(frequencies_1))
#   group.by.names <- group.by.names[group.by.names != replicate.by]
#
#   sinaplot_list <- lapply(1:length(group.by.names), function(i){
#     target_group <- group.by.names[i]
#     if(perform.t.test){
#       if(paired){
#         t_test_results <- t.test(frequencies_1[[target_group]], frequencies_2[[target_group]], paired = TRUE)
#       } else {
#         t_test_results <- t.test(frequencies_1[[target_group]], frequencies_2[[target_group]])
#       }
#     }
#     frequencies_1_temp <- data.frame(split_by = split.by.values[1],
#                                      replicate_by = frequencies_1[[replicate.by]],
#                                      Frequency = frequencies_1[[target_group]])
#     frequencies_2_temp <- data.frame(split_by = split.by.values[2],
#                                      replicate_by = frequencies_2[[replicate.by]],
#                                      Frequency = frequencies_2[[target_group]])
#
#     temp_df <- rbind(frequencies_1_temp, frequencies_2_temp)
#
#
#     g <- ggplot(temp_df, aes(x = split_by,
#                              y = Frequency,
#                              color = replicate_by,
#                              fill = split_by,
#                              group = replicate_by))+
#       geom_violin(inherit.aes = FALSE, aes(x = split_by,
#                                            y = Frequency,
#                                            group = split_by,
#                                            fill = split_by))+
#       ggbeeswarm::geom_beeswarm(size = 3)+
#       geom_line()+
#       scale_y_continuous(labels = function(x) paste0(round(x*100, 2), "%"))+
#       ggtitle(target_group)+
#       scale_fill_manual(values = c("black", "grey"))+
#       NoLegend()
#     if(perform.t.test){
#       g <- g + annotate("text", x = 1.5, y = 1.1 * max(temp_df$Frequency), label = paste0("p-value: ", round(t_test_results$p.value, 5)))
#     }
#     if(set_minimum_to_zero){
#       g <- g + coord_cartesian(ylim = c(0,NA))
#     }
#     g
#   })
#   cowplot::plot_grid(plotlist = sinaplot_list, ncol = ncols)
# }

AbundancePlotMulti <- function(seurat_object,
                               group.by,
                               split.by,
                               replicate.by,
                               draw_paths = FALSE,
                               perform.stat.test = TRUE,
                               test.choice = "t.test",
                               ncols = 3,
                               paired = FALSE,
                               set_minimum_to_zero = FALSE,
                               fill_colors = NULL,
                               point_colors = NULL,
                               p.adjust.method = "BH",
                               filter.p.above = 1,
                               step.increase = 0.15,
                               title_size = 10,
                               p_val_size = 2.5,
                               spacing_above = 1.3,
                               plot_type = "box",
                               sina_shift = TRUE,
                               x_lab = NULL,
                               y_lab = NULL,
                               point_size = 1,
                               same_y_limit = FALSE){

  x_lab <- if(is.null(x_lab)) split.by else x_lab
  y_lab <- if(is.null(y_lab)) "percentage" else y_lab

  seurat_metadata <- seurat_object@meta.data
  seurat_metadata_filtered <- seurat_metadata %>%
    dplyr::group_by(!!sym(replicate.by), !!sym(split.by), !!sym(group.by)) %>%
    dplyr::tally() %>%
    dplyr::mutate(Frequency = n/sum(n)) %>%
    dplyr::ungroup() %>%
    droplevels() %>%
    tidyr::complete(nesting(!!sym(replicate.by), !!sym(split.by)), !!sym(group.by), fill = list(Frequency = 0, n = 0))

  split.by.values <- unique(seurat_metadata_filtered[[split.by]])

  frequencies <- lapply(seq_along(split.by.values), function(i){
    seurat_metadata_filtered %>%
      dplyr::filter(!!sym(split.by) == split.by.values[i]) %>%
      dplyr::ungroup() %>%
      dplyr::select(!!sym(replicate.by), !!sym(group.by), Frequency) %>%
      dplyr::mutate(split_by = split.by.values[i], replicate_by = !!sym(replicate.by), group_by = !!sym(group.by))
  }) %>% do.call(rbind, .)

  group.by.names <- sort(unique(frequencies[["group_by"]]))

  if(is.null(point_colors)){
    point_colors <- "black"
  }
  if(is.null(fill_colors)){
    fill_colors <- viridis::viridis_pal()(length(split.by.values))
  }

  my_plot_list <- lapply(seq_along(group.by.names), function(i){
    target_group <- group.by.names[i]
    target_group_frequencies <- frequencies %>%
      dplyr::filter(group_by == target_group)
    if(perform.stat.test){
      t_test_results <- ggpubr::compare_means(Frequency~split_by,
                                              target_group_frequencies,
                                              paired = paired,
                                              p.adjust.method = p.adjust.method,
                                              method = test.choice)
    }

    #start ggplot
    g <- ggplot(target_group_frequencies, aes(x = split_by,
                                              y = Frequency,
                                              fill = split_by,
                                              group = split_by))+
      scale_fill_manual(values = fill_colors)+
      theme(title = element_text(size = title_size),
            panel.background = element_rect(fill = "white"),
            axis.line.x.bottom = element_line(color = 'black'),
            axis.line.y.left   = element_line(color = 'black'))+
      ggtitle(target_group)+
      xlab(x_lab)+
      NoLegend()

    #plot either violin or boxplot
    if(plot_type == "violin"){
      g <- g + geom_violin()
    } else if (plot_type == "box") {
      g <- g + geom_boxplot(outlier.shape = NA)
    } else {
      stop("plot_type must be one of violin or box")
    }


    #plot points as sina points or regular ggplot2 points
    if(sina_shift){
      g <- g + ggforce::geom_sina(color = point_colors, size = point_size)
    } else {
      g <- g + ggplot2::geom_point(color = point_colors, size = point_size)
    }

    #draw lines between replicates if applicable
    if(draw_paths){
      g <- g + geom_line(inherit.aes = FALSE, aes(x = split_by, y = Frequency, group = replicate_by))
    }

    #set plot minimum to 0 if needed (sometimes helpful)
    if(set_minimum_to_zero){
      ylim_min <- 0
    } else {
      ylim_min <- NA
    }

    #add brackets to statistical testing
    if(perform.stat.test){

      #get appropriate column and at the same time the prefix used
      p_column <- if(p.adjust.method == "none") "p" else "p.adj"

      #filtering if necessary
      t_test_results <- dplyr::filter(t_test_results, !!sym(p_column) <= filter.p.above)

      if(nrow(t_test_results) > 0){
        num_brackets <- nrow(t_test_results)
        max_frequency <- max(target_group_frequencies$Frequency)
        my_bracket_floor <- 1.05*max_frequency
        my_plot_ceiling <- max_frequency + my_bracket_floor*step.increase*num_brackets
        g <- g + ggpubr::geom_bracket(inherit.aes = FALSE,
                                      label.size = p_val_size,
                                      data = t_test_results,
                                      aes(xmin = group1,
                                          xmax = group2,
                                          label = paste0(p_column, ": ",signif(!!sym(p_column), 2))),
                                      y.position = my_bracket_floor,
                                      step.increase = step.increase)+
          scale_y_continuous(limits = c(ylim_min,my_plot_ceiling), labels = function(x) paste0(x*100, "%"), name = y_lab)
      } else {
        g <- g + scale_y_continuous(limits = c(ylim_min,NA), labels = function(x) paste0(x*100, "%"), name = y_lab)
      }

      #if no plotting of p values, do:
    } else {
      g <- g + scale_y_continuous(limits = c(ylim_min,NA), labels = function(x) paste0(x*100, "%"), name = y_lab)
    }

  })

  #get lowest y and highest y limit
  y_plot_limits <- lapply(my_plot_list, function(ggp) layer_scales(ggp)$y$get_limits())
  y_plot_limit_high <- max(sapply(y_plot_limits, "[[", 2))
  y_plot_limit_low <- min(sapply(y_plot_limits, "[[", 1))


  final_plot <- patchwork::wrap_plots(plotlist = my_plot_list, ncol = ncols)
  if(same_y_limit){
    final_plot <- final_plot &
      scale_y_continuous(labels = function(x) paste0(x*100, "%"),limits = c(y_plot_limit_low, y_plot_limit_high), name = y_lab)
  }
  final_plot
}

FindJaccard <- function(genelist,
                        distance.measure = "jaccard",
                        as.distance = FALSE,
                        return.tidy = TRUE){
  genelist_data <- stack(genelist) %>%
    magrittr::set_colnames(c("gene", "geneset")) %>%
    dplyr::select(geneset, gene) %>%
    dplyr::mutate(occurrence = 1)

  genelist_wide <- genelist_data %>%
    tidyr::pivot_wider(id_cols = geneset, names_from = gene, values_from = occurrence, values_fill = 0) %>%
    as.data.frame() %>%
    tibble::column_to_rownames(var = "geneset")

  genelist_matrix <- proxy::dist(genelist_wide, method = distance.measure, diag = TRUE, convert_similarities = as.distance)

  if(return.tidy){
    genelist_matrix <- as.matrix(genelist_matrix)
    diag(genelist_matrix) <- 1
    genelist_matrix <- as.data.frame(genelist_matrix) %>% tibble::rownames_to_column(var = "genelist_A")
    genelist_matrix <- tidyr::pivot_longer(genelist_matrix, !genelist_A, names_to = "genelist_B")
    return(genelist_matrix)
  }

  return(genelist_matrix)

}


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


FeaturePlotPlus <- function(seurat_object,
                            feature,
                            min.value = 0,
                            assay = "RNA",
                            viridis_choice = "D",
                            ...){
  DefaultAssay(seurat_object) <- assay
  temp_plot <- FeaturePlot(seurat_object, features = feature, ...)
  temp_plot$data[,feature] <- ifelse(temp_plot$data[,feature] <= min.value, NA, temp_plot$data[,feature])
  temp_plot +
    viridis::scale_color_viridis(option = viridis_choice, na.value = "grey")
}

DimPlotEdges <- function(seurat_object,
                         graph.choice = "wknn",
                         reduction = "umap",
                         line.thickness = 0.05,
                         line.color = "grey",
                         raster = FALSE,
                         ...){

  my_graph <- as(Graphs(seurat_object, graph.choice), "dgCMatrix")
  my_graph_coordinates <- data.frame(Matrix::summary(my_graph))
  my_graph_coordinates$first_cell <- my_graph@Dimnames[[1]][my_graph_coordinates$i]
  my_graph_coordinates$second_cell <- my_graph@Dimnames[[1]][my_graph_coordinates$j]

  my_dimreduc_coordinates <- Embeddings(seurat_object, reduction) %>%
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

  temp_dimplot <- DimPlot(seurat_object, reduction = reduction, raster = raster, ...)+
    ggrastr::rasterize(geom_segment(data = my_graph_coordinates, aes(x = coord_x.1, y = coord_y.1, xend = coord_x.2, yend = coord_y.2), inherit.aes = FALSE, size = line.thickness, color = line.color))

  num_layers <- length(temp_dimplot$layers)
  temp_dimplot$layers <- temp_dimplot$layers[c(num_layers, 1:(num_layers-1))]
  temp_dimplot
}

DimPlotPlus <- function(seurat_object,
                        reduction = "umap",
                        shape.choice = 21,
                        pt.size = 1.75,
                        st.size = 0.25,
                        ...){
  temp_dimplot <- DimPlot(seurat_object, reduction = reduction, ...)
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


ReadGMT <- function(file_path = "", file_sep = ","){
  #code adapted from clusterProfiler::read.gmt but fixed for comma-separated GMT file
  x <- readLines(file_path)
  res <- strsplit(x, file_sep)
  names(res) <- vapply(res, function(y) y[1], character(1))
  res <- lapply(res, "[", -c(1:2))
  ont2gene <- stack(res)
  ont2gene <- ont2gene[, c("ind", "values")]
  colnames(ont2gene) <- c("term", "gene")
  return(ont2gene)
}

BuildRealWKNN <- function(seurat_object, neighbors_name = "weighted.nn", new_graph_name = "real_wknn", weighted = FALSE){
  if(is.null(seurat_object@neighbors[[neighbors_name]])){
    stop(paste0("need to have a Neighbors object called ", neighbors_name, " in the Seurat object"))
  }
  my_neighbors_object <- seurat_object@neighbors[[neighbors_name]]
  my_wnn_edgelist <- cbind(1:nrow(my_neighbors_object@nn.idx), c(my_neighbors_object@nn.idx)) %>%
    set_colnames(c("i", "j")) %>%
    as.data.frame()
  if(weighted){
    my_wnn_edgelist$weight <- c(my_neighbors_object@nn.dist)
    my_wnn_adjacency <- igraph::get.adjacency(igraph::graph_from_data_frame(my_wnn_edgelist), attr = "weight")
  } else {
    my_wnn_adjacency <- igraph::get.adjacency(igraph::graph_from_data_frame(my_wnn_edgelist))
  }
  my_wnn_adjacency@Dimnames[[1]] <- colnames(seurat_object)
  my_wnn_adjacency@Dimnames[[2]] <- colnames(seurat_object)
  Matrix::diag(my_wnn_adjacency) <- 1
  seurat_object@graphs[[new_graph_name]] <- as.Graph(my_wnn_adjacency)
  return(seurat_object)
}

DoStarHeatmap <- function(seurat_object,
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

  my_data <- AverageExpression(seurat_object, slot = slot, assay = assay, group.by = group.by)[[assay]]
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


DoClusteredDotPlot <- function(seurat_object,
                               diff_exp_results = NULL,
                               features = NULL,
                               assay = "RNA",
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
                               dot.scale = 5,
                               cluster_cols = TRUE){

  my_data <- AverageExpression(seurat_object, slot = slot, assay = assay, group.by = group.by)[[assay]]

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
    if(!is.null(subset_features)){
      my_data <- my_data[subset_features,]
    }
    if(!plot_all){
      my_data <- my_data[rownames(my_data) %in% sig_genes,]
    }

  } else {
    if(!is.null(features)){
      my_data <- my_data[features,]
    } else {
      stop("most provide diff_exp_results or features")
    }
  }



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
  base_dotplot <- DotPlot(seurat_object,
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


DoClusteredHeatmap <- function(seurat_object,
                               features = NULL,
                               assay = "RNA",
                               slot = "scale.data",
                               group.by = "seurat_clusters",
                               max_zed = 3,
                               plot_rownames = FALSE,
                               label_features = TRUE,
                               y_text_size = 10,
                               scale_rows = FALSE,
                               p_val_choice = 0.01,
                               logFC = 0.4,
                               plot_all = FALSE,
                               plot_dendro = FALSE,
                               subset_features = NULL,
                               auc_choice = 0.6,
                               dot.scale = 5,
                               cluster_cols = TRUE,
                               cluster_rows = TRUE){

  if(slot == "scale.data" & !all(features %in% rownames(seurat_object[[assay]]@scale.data))){
    seurat_object <- ScaleData(seurat_object, assay = assay, features = features)
  }
  my_data <- AverageExpression(seurat_object, slot = slot, assay = assay, group.by = group.by)[[assay]]
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


#some inspo from the awesome sceasy github
SeuratToAnnData <- function(seurat_object,
                            outFile = NULL,
                            main.assay = 'RNA',
                            main.layer = 'counts',
                            transfer.other.main.layers = NULL,
                            transfer.other.assays.layers = NULL,
                            transfer.graphs = TRUE) {

  main.layer <- match.arg(main.layer, c('data', 'counts', 'scale.data'))
  transfer.other.main.layers <- transfer.other.main.layers[
    transfer.other.main.layers %in% c('data', 'counts', 'scale.data')]
  transfer.other.main.layers <- transfer.other.main.layers[transfer.other.main.layers != main.layer]

  message("getting main.layer data from main.assay and assigning to X")
  X <- Seurat::GetAssayData(object = seurat_object,
                            assay = main.assay,
                            slot = main.layer)

  message("getting metadata")
  obs <- seurat_object@meta.data

  message("getting main.assay metadata")
  var <- Seurat::GetAssay(seurat_object, assay = main.assay)@meta.features
  var <- cbind(var, feature_name = rownames(Seurat::GetAssayData(seurat_object, assay = main.assay)))

  message("getting all reductions and appending X_")
  obsm <- NULL
  reductions <- names(seurat_object@reductions)
  if (length(reductions) > 0) {
    obsm <- sapply(
      reductions,
      function(name) as.matrix(Seurat::Embeddings(seurat_object, reduction=name)),
      simplify = FALSE
    )
    names(obsm) <- paste0('X_', tolower(names(seurat_object@reductions)))
  }

  layers <- list()
  for (layer in transfer.other.main.layers) {
    message(paste0("transferring ", main.assay, " assay's ", layer,  " layer"))
    mat <- Seurat::GetAssayData(object = seurat_object, assay = main.assay, slot = layer)
    if (all(dim(mat) == dim(X))) layers[[layer]] <- Matrix::t(mat)
  }

  obsm.assays <- NULL
  if (length(transfer.other.assays.layers) > 0) {
    message("getting remaining assay:slot pairs and adding to obsm")
    obsm.assays <- lapply(
      seq_along(transfer.other.assays.layers),
      function(i){
        message(paste0("transferring ",names(transfer.other.assays.layers)[i], " ", transfer.other.assays.layers[i]))
        as.data.frame(Matrix::t((Seurat::GetAssayData(object = seurat_object,
                                                      assay = names(transfer.other.assays.layers)[i],
                                                      slot = transfer.other.assays.layers[[i]]))))
      }
    )
    names(obsm.assays) <- paste0("assay_", names(transfer.other.assays.layers), "_", transfer.other.assays.layers)
    obsm <- c(obsm, obsm.assays)
  }

  obsp <- NULL
  if(transfer.graphs){
    graph.list <- names(seurat_object@graphs)
    message("getting sparse Seurat Graphs from R to python")
    obsp.graphs <- lapply(seq_along(graph.list), function(i){
      message(paste0("getting ", graph.list[i]))
      reticulate::r_to_py(as(SeuratObject::Graphs(seurat_object, graph.list[i]), "dgCMatrix"))
    })
    names(obsp.graphs) <- graph.list
    obsp <- c(obsp, obsp.graphs)
  }

  message("loading in anndata through reticulate")
  anndata <- reticulate::import('anndata', convert = FALSE)


  adata <- anndata$AnnData(
    X = Matrix::t(X),
    obs = obs,
    var = var,
    obsm = obsm,
    obsp = obsp,
    layers = layers
  )

  if (!is.null(outFile))
    adata$write(outFile, compression = 'gzip')

  adata
}

RunSCA <- function(seurat_object,
                   assay = "RNA",
                   slot = "data",
                   reduction = "pca",
                   n.comps = 50L,
                   iters = 3L,
                   reduction.name = "sca",
                   reduction.key = "sca_"){
  starting.reduction <- Embeddings(seurat_object, reduction)
  X <- Matrix::t(GetAssayData(seurat_object, assay = assay, slot = slot))
  shannonca <- reticulate::import("shannonca", convert = FALSE)
  new.reduction = reticulate::py_to_r(shannonca$reduce(X, n_comps=n.comps, iters=iters))
  colnames(new.reduction) <- paste0("SC_", 1:ncol(new.reduction))
  rownames(new.reduction) <- colnames(seurat_object)
  new.reduction <- Seurat::CreateDimReducObject(
    embeddings = new.reduction,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  seurat_object[[reduction.name]] <- new.reduction
  return(seurat_object)
}

AssignQuantiles <- function(seurat_object,
                            feature,
                            assay = "ADT",
                            slot = "data",
                            split.by = "run_10x",
                            quantile.probs = c(0,.25,.5,.75,1)){
  split_objects <- SplitObject(seurat_object, split.by = split.by)
  split_objects <- lapply(seq_along(split_objects), function(i){
    temp_object <- split_objects[[i]]
    temp_vector <- temp_object[[assay]][feature,]
    temp_quantiles <- quantile(temp_vector, probs = quantile.probs)
    temp_cuts <- cut(temp_vector, breaks = temp_quantiles, include.lowest = TRUE)
    temp_cuts_levels <- levels(temp_cuts)
    temp_cuts <- forcats::fct_recode(temp_cuts, !!!setNames(temp_cuts_levels, paste0(feature, ".q", 1:length(temp_cuts_levels))))
    temp_object[[paste0(feature, ".quantile")]] <- temp_cuts
    return(temp_object)
  })
  merged_object <- merge(split_objects[[1]], split_objects[2:length(split_objects)])
  seurat_object[[paste0(feature, ".quantile")]] <- merged_object[[paste0(feature, ".quantile")]]
  return(seurat_object)
}


HotspotComputeAutocorrelations <- function(seurat_object,
                                           neighbors.graph = "RNA_nn_nn",
                                           weight.graph = "RNA_nn_cosine",
                                           latent.dimred = "pca",
                                           model = "danb",
                                           assay = "RNA",
                                           n.jobs = 4L,
                                           FDR.cutoff = 0.01,
                                           Z.cutoff = 0.1,
                                           FDR.cutoff.2 = 0.01,
                                           n.genes = 500,
                                           min.gene.threshold = 15

){

  #fetch graphs from Seurat
  message("extracting Seurat graphs")
  nn.graph <- SeuratObject::Graphs(seurat_object, slot = neighbors.graph)
  dist.graph <- SeuratObject::Graphs(seurat_object, slot = weight.graph)

  #convert Graph objects to Neighbors objects
  nn.neighbors <- SeuratObject::as.Neighbor(nn.graph)
  dist.neighbors <- SeuratObject::as.Neighbor(dist.graph)

  #infer total k neighbors and then set last index (for later use in python)
  k.neighbors <- ncol(nn.neighbors)
  last.k.index <- k.neighbors - 1

  #convert Neighbors object slots to dataframes
  nn.neighbors.idx <- nn.neighbors@nn.idx
  dist.neighbors.dist <- nn.neighbors@nn.dist

  #scikit in python indexes the nn starting at 0 lol
  nn.neighbors.idx <- nn.neighbors.idx - 1
  #convert nn.neighbors.idx from double (default) to integer matrix
  mode(nn.neighbors.idx) <- "integer"

  #start reply python reticulate session
  message("importing pandas")
  #Sys.setenv(RETICULATE_PYTHON = python.dir)
  #reticulate::use_python(python.dir)
  pd <- reticulate::import("pandas", delay_load = TRUE)

  #make pandas dataframe from nn.neighbors.idx
  #use cell barcodes are index
  #use 0:last.k.index as columns
  message("transferring Seurat graphs to dataframes in pandas")
  nn_idx_py <- pd$DataFrame(reticulate::r_to_py(nn.neighbors.idx),
                            index = colnames(seurat_object),
                            columns = 0:last.k.index)

  #make pandas dataframe from dist.neighbors.dist
  #use cell barcodes are index
  #use 0:last.k.index as columns
  dist_py <- pd$DataFrame(reticulate::r_to_py(dist.neighbors.dist),
                          index = colnames(seurat_object),
                          columns = 0:last.k.index)

  #before we create/run Hotspot we will get the counts from the
  #specified assay and the total number of counts per barcode
  assay_counts <- GetAssayData(seurat_object, assay = assay, slot = "counts")
  counts_py <- reticulate::r_to_py(as.data.frame(assay_counts))
  #we also begrudgingly import a dim red to intialize the hotspot object
  latent_dimred <- Embeddings(seurat_object, reduction = latent.dimred)
  latent_py <-  pd$DataFrame(reticulate::r_to_py(latent_dimred),
                             index = colnames(seurat_object),
                             columns = 0:(ncol(latent_dimred)-1))

  #import hotspot and set up hotspot object
  message("importing hotspot")
  hs <- reticulate::import("hotspot")

  message("initializing hotspot object")
  hs_object = hs$Hotspot(counts = counts_py,
                         model = model,
                         latent = latent_py)
  hs_object$neighbors = nn_idx_py
  hs_object$weights = dist_py

  message("computing hotspot autocorrelation")
  autocorrelation_results = hs_object$compute_autocorrelations(jobs = n.jobs)
  autocorrelation_results_filtered_genes <- autocorrelation_results %>%
    dplyr::filter(FDR < FDR.cutoff) %>%
    dplyr::top_n(n.genes, Z) %>%
    rownames()

  message("computing hotspot local correlation")
  local_correlations = hs_object$compute_local_correlations(genes = autocorrelation_results_filtered_genes, jobs = n.jobs)

  message("computing hotspot modules")
  my_modules = hs_object$create_modules(min_gene_threshold=min.gene.threshold,
                                        core_only=TRUE,
                                        fdr_threshold=FDR.cutoff.2)

  my_modules_df <- data.frame(gene = names(my_modules), module = setNames(my_modules,NULL))
  return(my_modules_df)

}



MiloGetIDs <- function(milo_object, milo_results, p_cutoff = 0.1, p_choice = "FDR"){
  milo_sig_results <- milo_results %>%
    dplyr::filter(!!sym(p_choice) < p_cutoff)
  if(nrow(milo_sig_results) < 1){
    stop("No statisticaly significant results in milo_results")
  }
  sig_nbhds_ix <- milo_sig_results[["Nhood"]]
  sig_nbhds_id <- as.numeric(nhoodIndex(milo_object)[sig_nbhds_ix])
  milo_nhoods <- nhoods(milo_object)
  all_ids <- 1:nrow(milo_nhoods)
  nbhd_list <- lapply(seq_along(sig_nbhds_ix), function(i){
    temp_nbhd_ix <- sig_nbhds_ix[i]
    chosen_indices <- all_ids[milo_nhoods[,temp_nbhd_ix] == 1]
    return(chosen_indices)
  })
  names(nbhd_list) <- sig_nbhds_id
  return(nbhd_list)

}


MiloWilcoxGenes <- function(seurat_object, milo_ids){
  milo_ids_names <- names(milo_ids)
  nbhd_genes <- lapply(seq_along(milo_ids), function(i){
    temp_nbhd_ids <- milo_ids[[i]]
    seurat_object[["tempID"]] <- ifelse(1:ncol(seurat_object) %in% temp_nbhd_ids,
                                        yes = "in_cluster", no = "out_cluster")
    presto_results <- presto::wilcoxauc(seurat_object,
                                        group_by = "tempID",
                                        groups_use = c("in_cluster", "out_cluster")) %>%
      dplyr::filter(group == "in_cluster")
    presto_results$nhd_ix <- i
    presto_results$nhd_id <- milo_ids_names[i]
    return(presto_results)
  }) %>%
    do.call(rbind, .)
  return(nbhd_genes)
}


MiloHeatmap <- function(seurat_object, features, milo_ids, slot = "data", assay = "RNA", scale_rows = TRUE, max_zed = Inf, y_text_size = 10){
  features <- unique(features)
  non_milo_ids <- unique(unlist(milo_ids))
  milo_ids[["outgroup"]] <- non_milo_ids
  DefaultAssay(seurat_object) <- assay
  seurat_object <- seurat_object[features,]
  nbhd_expr <- lapply(seq_along(milo_ids), function(i){
    seurat_object[["tempID"]] <- ifelse(1:ncol(seurat_object) %in% milo_ids[[i]], "in_cluster", "out_cluster")
    my_data <- AverageExpression(seurat_object, group.by = "tempID")[[assay]]
    my_data <- my_data[,"in_cluster", drop = FALSE]
    colnames(my_data) <- names(milo_ids)[[i]]
    return(my_data)
  }) %>% do.call(cbind, .)
  nbhd_expr <- as.data.frame(t(scale(t(nbhd_expr), center = TRUE, scale = TRUE)))
  nbhd_expr[nbhd_expr > max_zed] <- max_zed
  nbhd_expr[nbhd_expr < -max_zed] <- -max_zed
  hclust_results <- hclust(dist(nbhd_expr))
  feature_order <- hclust_results$labels[hclust_results$order]
  hclust_results_2 <- hclust(dist(t(nbhd_expr)))
  cluster_order <- hclust_results_2$labels[hclust_results_2$order]

  my_data_long <- nbhd_expr %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "feature") %>%
    pivot_longer(-feature, values_to = "Z_sc", names_to = "cluster") %>%
    dplyr::mutate(feature = factor(feature, levels = feature_order)) %>%
    dplyr::mutate(cluster = factor(cluster, levels = cluster_order))
  ggplot(my_data_long, aes(x = cluster, y = feature, fill = Z_sc))+
    geom_tile()+
    scale_fill_viridis_c()+
    ylab(NULL)+
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
}

#copied/pasted from @immunogenomics/symphony
RunHarmony2 <- function(
  object,
  group.by.vars,
  reduction = 'pca',
  dims.use = NULL,
  theta = NULL,
  lambda = NULL,
  sigma = 0.1,
  nclust = NULL,
  tau = 0,
  block.size = 0.05,
  max.iter.harmony = 10,
  max.iter.cluster = 20,
  epsilon.cluster = 1e-5,
  epsilon.harmony = 1e-4,
  plot_convergence = FALSE,
  verbose = TRUE,
  reference_values = NULL,
  reduction.save = "harmony",
  assay.use = 'RNA',
  project.dim = TRUE,
  ...
) {
  if (reduction == "pca") {
    tryCatch(
      embedding <- Seurat::Embeddings(object, reduction = "pca"),
      error = function(e) {
        if (verbose) {
          message("Harmony needs PCA. Trying to run PCA now.")
        }
        tryCatch(
          object <- Seurat::RunPCA(
            object,
            assay = assay.use, verbose = verbose
          ),
          error = function(e) {
            stop("Harmony needs PCA. Tried to run PCA and failed.")
          }
        )
      }
    )
  } else {
    available.dimreduc <- names(methods::slot(object = object, name = "reductions"))
    if (!(reduction %in% available.dimreduc)) {
      stop("Requested dimension reduction is not present in the Seurat object")
    }
    embedding <- Seurat::Embeddings(object, reduction = reduction)
  }
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(embedding))
  }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rereun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  metavars_df <- Seurat::FetchData(object, group.by.vars)

  harmonyObject <- HarmonyMatrix(
    embedding,
    metavars_df,
    group.by.vars,
    FALSE,
    0,
    theta,
    lambda,
    sigma,
    nclust,
    tau,
    block.size,
    max.iter.harmony,
    max.iter.cluster,
    epsilon.cluster,
    epsilon.harmony,
    plot_convergence,
    TRUE,
    verbose,
    reference_values
  )

  harmonyEmbed <- t(as.matrix(harmonyObject$Z_corr))
  rownames(harmonyEmbed) <- row.names(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.save, "_", seq_len(ncol(harmonyEmbed)))

  harmonyClusters <- t(harmonyObject$R)
  rownames(harmonyClusters) <- row.names(embedding)
  colnames(harmonyClusters) <- paste0('R', seq_len(ncol(harmonyClusters)))

  suppressWarnings({
    harmonydata <- Seurat::CreateDimReducObject(
      embeddings = harmonyEmbed,
      stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
      assay = assay.use,
      key = reduction.save,
      misc=list(R=harmonyClusters)
    )
  })

  object[[reduction.save]] <- harmonydata
  if (project.dim) {
    object <- Seurat::ProjectDim(
      object,
      reduction = reduction.save,
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  return(object)
}

ExtractClusterProportions <- function(seurat_object,
                                      group.by,
                                      split.by,
                                      replicate.by){

  seurat_metadata <- seurat_object@meta.data
  seurat_metadata_filtered <- seurat_metadata %>%
    dplyr::group_by(!!sym(replicate.by), !!sym(split.by), !!sym(group.by)) %>%
    dplyr::tally() %>%
    dplyr::mutate(percent = n/sum(n)) %>%
    dplyr::ungroup() %>%
    droplevels() %>%
    tidyr::complete(nesting(!!sym(replicate.by), !!sym(split.by)), !!sym(group.by), fill = list(percent = 0, n = 0))

  split.by.values <- unique(seurat_metadata_filtered[[split.by]])

  frequencies <- lapply(seq_along(split.by.values), function(i){
    seurat_metadata_filtered %>%
      dplyr::filter(!!sym(split.by) == split.by.values[i]) %>%
      dplyr::ungroup() %>%
      dplyr::select(!!sym(replicate.by), !!sym(group.by), percent) %>%
      dplyr::mutate(split_by = split.by.values[i], replicate_by = !!sym(replicate.by), group_by = !!sym(group.by))
  }) %>% do.call(rbind, .)
  return(frequencies)
}



ModulePlotMulti <- function(seurat_object,
                            assay = "rhAUC",
                            features = rownames(seurat_object[[assay]]),
                            split.by,
                            replicate.by,
                            draw_paths = FALSE,
                            perform.stat.test = TRUE,
                            test.choice = "t.test",
                            ncols = 3,
                            paired = FALSE,
                            set_minimum_to_zero = FALSE,
                            fill_colors = NULL,
                            point_colors = NULL,
                            p.adjust.method = "BH",
                            filter.p.above = 1,
                            step.increase = 0.15,
                            title_size = 10,
                            p_val_size = 2.5,
                            spacing_above = 1.3,
                            plot_type = "box",
                            sina_shift = TRUE,
                            x_lab = NULL,
                            y_lab = NULL,
                            same_y_limit = FALSE){

  x_lab <- if(is.null(x_lab)) split.by else x_lab
  y_lab <- if(is.null(y_lab)) "avg_expression" else y_lab

  my_splitter <- c(split.by, replicate.by)

  #nb the Seurat AverageExpression function returns colnames separated by a "_"
  avg_module_exp <- AverageExpression(seurat_object, assays = assay, group.by = my_splitter)[[1]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "module") %>%
    tidyr::pivot_longer(-module, values_to = "avg_score") %>%
    tidyr::separate(name, into = my_splitter, sep = "_") %>%
    tidyr::complete(nesting(!!sym(replicate.by), !!sym(split.by)), fill = list(avg_score = 0, n = 0)) %>%
    dplyr::mutate(split_by = !!sym(split.by), replicate_by = !!sym(replicate.by))

  facet.by.values <- features
  split.by.values <- unique(avg_module_exp$split_by)

  if(is.null(point_colors)){
    point_colors <- "black"
  }
  if(is.null(fill_colors)){
    fill_colors <- viridis::viridis_pal()(length(split.by.values))
  }

  my_plot_list <- lapply(seq_along(facet.by.values), function(i){
    target_group <- facet.by.values[i]
    target_group_scores <- avg_module_exp %>%
      dplyr::filter(module == target_group)
    if(perform.stat.test){
      t_test_results <- ggpubr::compare_means(avg_score~split_by,
                                              target_group_scores,
                                              paired = FALSE,
                                              p.adjust.method = p.adjust.method,
                                              method = test.choice)
    }
    #start ggplot
    g <- ggplot(target_group_scores, aes(x = split_by,
                                         y = avg_score,
                                         fill = split_by,
                                         group = split_by))+
      scale_fill_manual(values = fill_colors)+
      theme(title = element_text(size = title_size),
            panel.background = element_rect(fill = "white"),
            axis.line.x.bottom = element_line(color = 'black'),
            axis.line.y.left   = element_line(color = 'black'))+
      ggtitle(target_group)+
      xlab(x_lab)+
      NoLegend()

    #plot either violin or boxplot
    if(plot_type == "violin"){
      g <- g + geom_violin()
    } else if (plot_type == "box") {
      g <- g + geom_boxplot(outlier.shape = NA)
    } else {
      stop("plot_type must be one of violin or box")
    }


    #plot points as sina points or regular ggplot2 points
    if(sina_shift){
      g <- g + ggforce::geom_sina(color = point_colors)
    } else {
      g <- g + ggplot2::geom_point(color = point_colors)
    }

    #draw lines between replicates if applicable
    if(draw_paths){
      g <- g + geom_line(inherit.aes = FALSE, aes(x = split_by, y = avg_score, group = replicate_by))
    }

    #set plot minimum to 0 if needed (sometimes helpful)
    if(set_minimum_to_zero){
      ylim_min <- 0
    } else {
      ylim_min <- NA
    }

    #add brackets to statistical testing
    if(perform.stat.test){

      #get appropriate column and at the same time the prefix used
      p_column <- if(p.adjust.method == "none") "p" else "p.adj"

      #filtering if necessary
      t_test_results <- dplyr::filter(t_test_results, !!sym(p_column) <= filter.p.above)

      if(nrow(t_test_results) > 0){
        num_brackets <- nrow(t_test_results)
        max_score <- max(target_group_scores$avg_score)
        my_bracket_floor <- 1.05*max_score
        my_plot_ceiling <- max_score + my_bracket_floor*step.increase*num_brackets
        g <- g + ggpubr::geom_bracket(inherit.aes = FALSE,
                                      label.size = p_val_size,
                                      data = t_test_results,
                                      aes(xmin = group1,
                                          xmax = group2,
                                          label = paste0(p_column, ": ",signif(!!sym(p_column), 2))),
                                      y.position = my_bracket_floor,
                                      step.increase = step.increase)+
          scale_y_continuous(limits = c(ylim_min,my_plot_ceiling))
      } else {
        g <- g + scale_y_continuous(limits = c(ylim_min,NA))
      }

      #if no plotting of p values, do:
    } else {
      g <- g + scale_y_continuous(limits = c(ylim_min,NA))
    }

  })

  #get lowest y and highest y limit
  y_plot_limits <- lapply(my_plot_list, function(ggp) layer_scales(ggp)$y$get_limits())
  y_plot_limit_high <- max(sapply(y_plot_limits, "[[", 2))
  y_plot_limit_low <- min(sapply(y_plot_limits, "[[", 1))

  final_plot <- patchwork::wrap_plots(plotlist = my_plot_list, ncol = ncols)
  if(same_y_limit){
    final_plot <- final_plot &
      scale_y_continuous(limits = c(y_plot_limit_low, y_plot_limit_high), name = y_lab)
  } else {
    final_plot <- final_plot &
      scale_y_continuous(name = y_lab)
  }
  final_plot

}



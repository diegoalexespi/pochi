#code adapted from @immunogenomics/singlecellmethods
FindClustersLeiden <- function(seurat_object,
                               graph_name = NULL,
                               resolution_list = list(0.8),
                               min_cluster_size = 5,
                               verbose = TRUE,
                               pythondir = "~/miniconda3/bin/python",
                               partition_type = "RBConfigurationVertexPartition") {


  reticulate::use_python(pythondir)
  leiden <- reticulate::import("leidenalg")
  ig <- reticulate::import("igraph")

  tmpfile <- tempfile()

  if (!graph_name %in% names(x = seurat_object)) {
    stop("Provided graph.name not present in Seurat object")
  }
  if (!is(object = seurat_object[[graph_name]], class2 = "Graph")) {
    stop("Provided graph.name does not correspond to a graph object.")
  }

  snn <- seurat_object[[graph_name]]
  snn <- as(snn, "dgTMatrix")
  data.table::data.table(i = snn@i, j = snn@j) %>% data.table::fwrite(tmpfile, sep = " ", col.names = FALSE)
  g <- ig$Graph(directed = FALSE)
  g <- g$Read_Edgelist(tmpfile)

  .res <- Reduce(cbind, lapply(resolution_list, function(resolution) {
    cres <- leiden$find_partition(g, leiden[[partition_type]], resolution_parameter = resolution)$membership
    clusters_keep <- names(table(cres))[which(table(cres) >= min_cluster_size)]
    cres[which(!cres %in% clusters_keep)] <- NA
    n_na <- sum(is.na(cres))
    message(sprintf("Resolution %f yielded %d clusters", resolution, length(clusters_keep)))
    if (verbose & n_na > 0) message(sprintf("WARNING: for resolution %f, %d/%d points removed from small clusters", resolution, n_na, length(cres)))
    return(cres)
  }))
  .res <- as.data.frame(.res)
  rownames(.res) <- colnames(seurat_object)
  colnames(.res) <- paste0(DefaultAssay(seurat_object), "_snn_leiden.", resolution_list)

  seurat_object <- AddMetaData(object = seurat_object, metadata = .res)
  return(seurat_object)
}

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
  return_graphs = TRUE,
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

  if(return_graphs) graph_vector <- c("nn", "fgraph") else graph_vector <- c()

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
    scale = scale
  )


  if(return_graphs){
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

    umap_output <- umap_results$embedding

  } else {
    umap_output <- umap_results
  }

  colnames(x = umap_output) <- paste0(reduction_key, 1:ncol(x = umap_output))
  rownames(x = umap_output) <- rownames(data_use)
  umap_reduction <- Seurat::CreateDimReducObject(
    embeddings = umap_output,
    key = reduction_key,
    assay = assay,
    global = TRUE
  )
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

AbundancePlot <- function(seurat_object,
                          group.by,
                          split.by,
                          replicate.by,
                          perform.t.test = TRUE,
                          ncols = 3,
                          paired = FALSE,
                          set_minimum_to_zero = FALSE){

  seurat_metadata <- seurat_object@meta.data
  seurat_metadata_filtered <- seurat_metadata %>%
    dplyr::group_by(!!sym(replicate.by), !!sym(split.by), !!sym(group.by)) %>%
    dplyr::tally() %>%
    dplyr::mutate(Frequency = n/sum(n)) %>%
    dplyr::ungroup() %>%
    droplevels() %>%
    tidyr::complete(nesting(!!sym(replicate.by), !!sym(split.by)), !!sym(group.by), fill = list(Frequency = 0, n = 0))

  split.by.values <- unique(seurat_metadata_filtered[[split.by]])
  if(length(split.by.values) != 2){
    stop("split.by must denote a column of only 2 groups")
  }

  frequencies_1 <- seurat_metadata_filtered %>%
    dplyr::filter(!!sym(split.by) == split.by.values[1]) %>%
    dplyr::ungroup() %>%
    dplyr::select(!!sym(replicate.by), !!sym(group.by), Frequency) %>%
    tidyr::pivot_wider(names_from = !!sym(group.by), values_from = Frequency, values_fill = 0)

  frequencies_2 <- seurat_metadata_filtered %>%
    dplyr::filter(!!sym(split.by) == split.by.values[2]) %>%
    dplyr::ungroup() %>%
    dplyr::select(!!sym(replicate.by), !!sym(group.by), Frequency) %>%
    tidyr::pivot_wider(names_from = !!sym(group.by), values_from = Frequency, values_fill = 0)

  group.by.names <- sort(colnames(frequencies_1))
  group.by.names <- group.by.names[group.by.names != replicate.by]

  boxplot_list <- lapply(1:length(group.by.names), function(i){
    target_group <- group.by.names[i]
    if(perform.t.test){
      if(paired){
        t_test_results <- t.test(frequencies_1[[target_group]], frequencies_2[[target_group]], paired = TRUE)
      } else {
        t_test_results <- t.test(frequencies_1[[target_group]], frequencies_2[[target_group]])
      }
    }
    frequencies_1_temp <- data.frame(split_by = split.by.values[1],
                                     replicate_by = frequencies_1[[replicate.by]],
                                     Frequency = frequencies_1[[target_group]])
    frequencies_2_temp <- data.frame(split_by = split.by.values[2],
                                     replicate_by = frequencies_2[[replicate.by]],
                                     Frequency = frequencies_2[[target_group]])

    temp_df <- rbind(frequencies_1_temp, frequencies_2_temp)


    g <- ggplot(temp_df, aes(x = split_by,
                             y = Frequency,
                             color = replicate_by,
                             fill = split_by,
                             group = replicate_by))+
      geom_violin(inherit.aes = FALSE, aes(x = split_by,
                                           y = Frequency,
                                           group = split_by,
                                           fill = split_by))+
      geom_point(size = 3)+
      geom_line()+
      scale_y_continuous(labels = function(x) paste0(round(x*100, 2), "%"))+
      ggtitle(target_group)+
      scale_fill_manual(values = c("black", "grey"))+
      NoLegend()
    if(perform.t.test){
      g <- g + annotate("text", x = 1.5, y = 1.1 * max(temp_df$Frequency), label = paste0("p-value: ", round(t_test_results$p.value, 5)))
    }
    if(set_minimum_to_zero){
      g <- g + coord_cartesian(ylim = c(0,NA))
    }
    g
  })
  cowplot::plot_grid(plotlist = boxplot_list, ncol = ncols)
}


AbundancePlotMulti <- function(seurat_object,
                               group.by,
                               split.by,
                               replicate.by,
                               draw_paths = FALSE,
                               perform.t.test = TRUE,
                               ncols = 3,
                               paired = FALSE,
                               set_minimum_to_zero = FALSE){

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

  sinaplot_list <- lapply(seq_along(group.by.names), function(i){
    target_group <- group.by.names[i]
    target_group_frequencies <- frequencies %>%
      dplyr::filter(group_by == target_group)
    if(perform.t.test){
      t_test_results <- ggpubr::compare_means(Frequency~split_by,
                                              target_group_frequencies,
                                              paired = FALSE,
                                              method = "t.test")
    }

    g <- ggplot(target_group_frequencies, aes(x = split_by,
                                              y = Frequency,
                                              group = replicate_by))+
      geom_violin(inherit.aes = FALSE, aes(x = split_by,
                                           y = Frequency,
                                           group = split_by,
                                           fill = split_by))+
      scale_fill_manual(values = RColorBrewer::brewer.pal(length(group.by.names), "Greys"))+
      ggbeeswarm::geom_beeswarm()+
      scale_y_continuous(labels = function(x) paste0(round(x*100, 2), "%"))+
      ggtitle(target_group)+
      NoLegend()
    if(draw_paths){
      g <- g + geom_line()
    }
    if(perform.t.test){
      g <- g + ggpubr::geom_bracket(inherit.aes = FALSE,
                                    data = t_test_results,
                                    aes(xmin = group1, xmax = group2, label = paste0("padj ",signif(p.adj, 2))),
                                    y.position = max(target_group_frequencies$Frequency*1.1),
                                    step.increase = 0.15)
    }
    if(set_minimum_to_zero){
      g <- g + coord_cartesian(ylim = c(0,NA))
    }
    g
  })
  cowplot::plot_grid(plotlist = sinaplot_list, ncol = ncols)
}

FindJaccard <- function(genelist, distance.measure = "jaccard", as.distance = FALSE, return.tidy = TRUE){
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

#FeaturePlot2 w grey 0 and assay choice
#Density plot between conditions
FeaturePlotPlus <- function(seurat_object, feature, min.value = 0, assay = "RNA", viridis_choice = "D",...){
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

BuildRealWKNN <- function(seurat_object, new_graph_name = "real_wknn"){
  if(is.null(seurat_object@neighbors$weighted.nn)){
    stop("need to have a Neighbors object called weighted.nn in the Seurat object")
  }
  my_neighbors_object <- seurat_object@neighbors$weighted.nn
  my_wnn_edgelist <- cbind(1:nrow(my_neighbors_object@nn.idx), c(my_neighbors_object@nn.idx)) %>%
    set_colnames(c("i", "j"))
  my_wnn_adjacency <- igraph::get.adjacency(igraph::graph.edgelist(my_wnn_edgelist))
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



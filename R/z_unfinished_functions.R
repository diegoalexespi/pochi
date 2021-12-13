FindClustersLeiden <- function(seurat_object,
                               graph_name = NULL,
                               resolution_list = list(0.8),
                               min_cluster_size = 5,
                               verbose = TRUE,
                               pythondir = "~/miniconda3/bin/python",
                               partition_type = "RBConfigurationVertexPartition") {


  require(reticulate)
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



MiloGetIDs <- function(milo_object, milo_results, nhood_column = "Nhood", p_cutoff = 0.1, p_choice = "FDR"){

  milo_sig_results <- milo_results %>%
    dplyr::filter(!!sym(p_choice) < p_cutoff)
  if(nrow(milo_sig_results) < 1){
    stop("No statisticaly significant results in milo_results")
  }
  sig_nbhds_id <- as.numeric(as.character(milo_sig_results[[nhood_column]]))
  sig_nbhds_ix <- which(unlist(as.numeric(nhoodIndex(milo_object))) == sig_nbhds_id)
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

SeuratAddNhds <- function(seurat_object, nhd_list, nhd_directions, milo_results_name = "milo_results"){
  nhd_df <- stack(nhd_list) %>%
    dplyr::mutate(nhd = as.numeric(as.character(ind)), cells = as.numeric(values)) %>%
    dplyr::select(nhd, cells)
  nhd_df$milo_results <- plyr::mapvalues(nhd_df$nhd, from = nhd_directions$nhd, to = nhd_directions$nhd_dir)
  nhd_dups <- nhd_df %>%
    dplyr::group_by(ind, directions) %>%
    dplyr::summarize()
  return(nhd_dups)
  #seurat_object[[milo_results_name]] <- ifelse(colnames(seurat_object) %in% nhd_df$)
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
                                           min.gene.threshold = 15,
                                           return_list = FALSE){

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
  if(return_list){
    return(list(local_correlations = as.matrix(local_correlations), modules = my_modules_df))
  } else {
    return(my_modules_df)
  }


}




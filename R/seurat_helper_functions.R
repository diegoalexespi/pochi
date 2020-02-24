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

#just a wrapper around uwot to return the kNN (and SNN) that the package constructs, which isn't done in Seurat
BuildAnnoyUMAP <- function(
  seurat_object,
  dims,
  use_raw = FALSE,
  genes_exclude = "",
  reduction = "pca",
  graph_key = "annoy",
  assay = NULL,
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
  init_pos = "spectral") {
  if (!is.null(x = seed_use)) {
    set.seed(seed = seed_use)
  }
  
  if(use_raw){
    data_use <- t(as.matrix(Seurat::GetAssayData(seurat_object, slot = slot, assay = assay)))
    data_use <- data.use[,!(colnames(data_use) %in% genes_exclude)]
  } else {
    data_use <- Seurat::Embeddings(seurat_object[[reduction]])[, dims]
  }
  
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
    ret_nn = TRUE
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
    
    nn_list <- list(nn_metric_sparse_dist, nn_metric_sparse_adj, nn_metric_sparse_snn)
    
    for(i in 1:length(nn_list)){
      colnames(nn_list[[i]]) <- colnames(seurat_object)
      rownames(nn_list[[i]]) <- colnames(seurat_object)
    }
    
    seurat_object[[paste0(graph_key, paste0("_nn_", metric))]] <- Seurat::as.Graph(nn_list[[1]])
    seurat_object[[paste0(graph_key, "_nn_nn")]] <- Seurat::as.Graph(nn_list[[2]])
    seurat_object[[paste0(graph_key, "_nn_snn")]] <- Seurat::as.Graph(nn_list[[3]])
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
  seurat_object[[reduction_name]] <- umap_reduction
  return(seurat_object)
}


#code adapted from @immunogenomics/singlecellmethods
BuildSNN_fromKNN <- function(nn_idx,
                             prune_snn,
                             nn_k = 30) {
  adj <- Matrix::sparseMatrix(i = rep(1:nrow(nn_idx), each = ncol(nn_idx)), 
                              j = c(t(nn_idx)), 
                              x = rep(1, prod(dim(nn_idx))))
  snn <- Matrix::tcrossprod(adj)
  snn@x <- snn@x / (2 * nn_k - snn@x)
  snn@x[snn@x < (prune_snn)] <- 0
  snn
}



#very naive version - may crash in some user cases
ToMonocle3 <- function(seurat_object,
                       scale_all = FALSE,
                       assay = "SCT",
                       reduction_for_projection = "PCA",
                       UMAP_cluster_slot = NULL){
  
  if(scale_all){
    message("Getting residuals for all Seurat genes in chosen assay slot and placing in scale.data")
    seurat_genes <- rownames(seurat_object[[assay]])
    remaining_genes <- setdiff(seurat_genes, rownames(seurat_object[[assay]]@scale.data))
    if(assay == "SCT"){
      seurat_object <- Seurat::GetResidual(seurat_object, features = remaining_genes, assay = assay, umi.assay = "RNA")
    } else {
      seurat_object <- Seurat::ScaleData(seurat_object, features = rownames(seurat_object[[assay]]))
    }
  }
  
  #We prep the seurat object by creating gene loadings for ALL genes in the Seurat scale.data slot. This is done to allow downstream monocle3 functions on gene_modules to work appropriately.
  message("Projecting gene loadings for all Seurat genes in scale.data slot")
  seurat_object <- Seurat::ProjectDim(seurat_object, reduction = reduction_for_projection, assay = assay)
  
  ##################
  
  message("Initializing CDS object")
  
  #Extract Seurat's log-transformed values
  expression_matrix <- Seurat::GetAssayData(seurat_object, assay = assay, slot = "counts")
  #Extract Seurat meta_data
  meta_data <- seurat_object@meta.data
  #Extract gene names from Seurat object SCT slot to make CDS
  seurat_genes <- data.frame(gene_short_name = rownames(seurat_object[[assay]]),
                             row.names = rownames(seurat_object[[assay]]))
  new_cds <- monocle3::new_cell_data_set(expression_data = expression_matrix, cell_metadata = meta_data, gene_metadata = seurat_genes)
  
  ##################
  
  message("Making an SCE object from the Seurat object to facilitate transfer of information from SCE to CDS")
  sce <- as.SingleCellExperiment(seurat_object, assay = assay)
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




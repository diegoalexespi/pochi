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

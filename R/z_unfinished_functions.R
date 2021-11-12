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

#' Calculate Leiden clustering for a Graph object in Seurat
#'
#' @param object Seurat object
#' @param graph.name Name of Graph slot in object to use for Leiden clustering
#' @param resolution Vector of resolution to input into Leiden
#' @param min.cluster.size Minimum cluster size trimmed after clustering
#' @param n.iterations Number of iterations to run Leiden clustering for
#' @param pythondir Specified director for python binary that has the leidenalg
#' library installed
#' @param partition_type The partition_type to specify for Leiden.
#' @details code is adapted from @immunogenomics/singlecellmethods with some
#' tweaks and flexibility added to match the arguments allowed by the
#' leidenalg package. In essence, we are using reticulate to run
#' leidenalg.find_partition(graph.name, partition_type) as specified here
#' https://leidenalg.readthedocs.io/en/stable/intro.html. Arguments are input
#' into the find_partition function here:
#' https://leidenalg.readthedocs.io/en/stable/reference.html#leidenalg.find_partition.
#' Singletons will be excluded for now, but ideally can be grouped in a step
#' afterwards...
#'
#' @return Returns a Seurat object with the leiden clusterings stored as
#' object@meta.data columns
#' @references code is adapted from @immunogenomics/singlecellmethods
#'
#' @importFrom rlang %||%
#'
#' @export
FindClustersLeiden <- function(object,
                               graph.name = NULL,
                               resolution = c(0.8),
                               min.cluster.size = 5,
                               n_iterations = 2,
                               pythondir = "~/miniconda3/bin/python",
                               partition_type = "RBConfigurationVertexPartition") {

  message("Warning: singleton grouping not implemented yet")
  message("Loading reticulate")
  require(reticulate)
  reticulate::use_python(pythondir)
  message("Importing leidenalg")
  leiden <- reticulate::import("leidenalg")
  message("Importing igraph")
  ig <- reticulate::import("igraph")

  tmpfile <- tempfile()

  if (!graph.name %in% names(x = object)) {
    stop("Provided graph.name not present in Seurat object")
  }
  if (!is(object = object[[graph.name]], class2 = "Graph")) {
    stop("Provided graph.name does not correspond to a graph object.")
  }

  message("Transferring Graph ", graph.name, " to temp file for python")
  snn <- object[[graph.name]]
  snn <- as(snn, "dgTMatrix")
  # make it a named vertex data table
  # remember that the dgTMatrix has 0-based indexing
  temp_table <- data.table::data.table(i = snn@Dimnames[[1]][snn@i+1],
                                       j = snn@Dimnames[[2]][snn@j+1],
                                       weight = snn@x) %>%
    data.table::fwrite(tmpfile, sep = " ", col.names = FALSE)

  message("Reading temp graph file into python")
  #read it in using Read_Ncol
  g <- ig$Graph$Read_Ncol(tmpfile, names = TRUE, weights = "if_present", directed = FALSE)


  .res_list <- lapply(resolution, function(resolution) {
    message("Applying Leiden algorithm for resolution: ", resolution)
    cres <- leiden$find_partition(g,
                                  leiden[[partition_type]],
                                  resolution_parameter = resolution,
                                  n_iterations = n_iterations)
    cres_membership <- data.frame(nodes = cres$graph$vs["name"], membership = cres$membership)
    clusters_keep <- names(table(cres_membership$membership))[which(table(cres_membership$membership) >= min.cluster.size)]
    cres_membership$membership[which(!cres_membership$membership %in% clusters_keep)] <- NA
    n_na <- sum(is.na(cres_membership$membership))
    message(sprintf("Resolution %f yielded %d clusters", resolution, length(clusters_keep)))
    if (n_na > 0) message(sprintf("WARNING: for resolution %f, %d/%d points removed from small clusters", resolution, n_na, length(cres_membership$membership)))
    colnames(cres_membership)[2] <- paste0("leiden_res.", resolution)
    return(cres_membership)
  })
  .res <-  .res_list %>%
    purrr::reduce(left_join, by = "nodes")

  message("Finished leiden clustering")
  rownames(.res) <- .res$nodes
  .res$nodes <- NULL
  .res <- .res[colnames(object),]

  object <- AddMetaData(object = object, metadata = .res)
  return(object)
}

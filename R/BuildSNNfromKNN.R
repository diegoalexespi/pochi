#' Derive an SNN from a kNN Graph object in Seurat
#'
#' @param object Seurat object
#' @param knn_choice Name of Graph slot to use for SNN calculation
#' @param snn_name Name of Graph slot to store SNN
#' @param prune_snn Threshold at which to prune SNN (set to 0 below this)
#' @details code is adapted from @immunogenomics/singlecellmethods. Builds a
#' shared nearest-neighbors graph from an input kNN using matrix multiplication
#' so pretty quick.
#'
#' @return Returns a Seurat object with the SNN Graph stored.
#' @references code is adapted from @immunogenomics/singlecellmethods
#'
#' @importFrom rlang %||%
#'
#' @export
BuildSNNfromKNN <- function(object,
                            knn_choice = "RNA_nn",
                            snn_name = paste0(knn_choice, "_snn"),
                            prune_snn = 1/15) {
  my_knn <- as(Graphs(seurat_object, knn_choice), "dgCMatrix")

  nn_k <- unique(Matrix::rowSums(my_knn))
  if(length(nn_k) != 1){
    stop("your kNN graph doesn't have equal NN for each node")
  }

  message("building SNN from ", knn_choice)
  my_snn <- Matrix::tcrossprod(as(my_knn, "dgCMatrix"))
  my_snn@x <- my_snn@x / (2 * nn_k - my_snn@x)
  my_snn@x[my_snn@x < (prune_snn)] <- 0
  my_snn <- Matrix::drop0(my_snn)
  seurat_object[[snn_name]] <- as.Graph(my_snn)
  return(seurat_object)
}

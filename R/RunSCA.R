#' Perform Shannon Component Analysis on a Seurat object
#'
#' @param object Seurat object
#' @param assay Which assay to use for SCA calculation
#' @param slot Which slot in the assay to use for SCA calculation
#' @param n.comps Number of SCs to calculate
#' @param iters How many iterations to perform for SCA
#' @param reduction.name Name of the Reduction object to store the SCA as
#' @param reduction.key Key preceding the numbers of SCA components in
#' reduction.name
#'
#' @return Returns a Seurat object with an SCA calculated
#' @references Demeo and Berger, _biorxiv_ 2021. shannonca.readthedocs.io
#'
#' @importFrom rlang %||%
#'
#' @export
RunSCA <- function(object,
                   assay = "RNA",
                   slot = "data",
                   n.comps = 50L,
                   iters = 3L,
                   reduction.name = "sca",
                   reduction.key = "sca_"){
  X <- Matrix::t(GetAssayData(object, assay = assay, slot = slot))
  shannonca <- reticulate::import("shannonca", convert = FALSE)
  new.reduction = reticulate::py_to_r(shannonca$reduce(X, n_comps=n.comps, iters=iters))
  colnames(new.reduction) <- paste0("SC_", 1:ncol(new.reduction))
  rownames(new.reduction) <- colnames(object)
  new.reduction <- Seurat::CreateDimReducObject(
    embeddings = new.reduction,
    key = reduction.key,
    assay = assay,
    global = TRUE
  )
  object[[reduction.name]] <- new.reduction
  return(object)
}

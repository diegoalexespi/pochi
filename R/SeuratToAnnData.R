#' Save a Seurat object as an anndata h5ad file
#'
#' @param object Seurat object
#' @param outFile Name of the target file to write the anndata object
#' @param main.assay Assay to transfer to AnnData as X
#' @param main.layer Layer to transfer to AnnData as X
#' @param transfer.other.main.layers Which other layers to transfer into AnnData
#' @param transfer.other.assays.layers Which other assay-layer pairs
#' to transfer into AnnData (since AnnData only takes one assay for now)
#' @param transfer.graphs Whether to transfer Graphs to AnnData .uns
#'
#' @return Writes a Seurat object to an .h5ad file after anndata conversion
#'
#' @importFrom rlang %||%
#' @references Code adapted from the awesome sceasy github!
#' @export
SeuratToAnnData <- function(object,
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
  X <- Seurat::GetAssayData(object = object,
                            assay = main.assay,
                            slot = main.layer)

  message("getting metadata")
  obs <- object@meta.data

  message("getting main.assay metadata")
  var <- Seurat::GetAssay(object, assay = main.assay)@meta.features
  var <- cbind(var, feature_name = rownames(Seurat::GetAssayData(object, assay = main.assay)))

  message("getting all reductions and appending X_")
  obsm <- NULL
  reductions <- names(object@reductions)
  if (length(reductions) > 0) {
    obsm <- sapply(
      reductions,
      function(name) as.matrix(Seurat::Embeddings(object, reduction=name)),
      simplify = FALSE
    )
    names(obsm) <- paste0('X_', tolower(names(object@reductions)))
  }

  layers <- list()
  for (layer in transfer.other.main.layers) {
    message(paste0("transferring ", main.assay, " assay's ", layer,  " layer"))
    mat <- Seurat::GetAssayData(object = object, assay = main.assay, slot = layer)
    if (all(dim(mat) == dim(X))) layers[[layer]] <- Matrix::t(mat)
  }

  obsm.assays <- NULL
  if (length(transfer.other.assays.layers) > 0) {
    message("getting remaining assay:slot pairs and adding to obsm")
    obsm.assays <- lapply(
      seq_along(transfer.other.assays.layers),
      function(i){
        message(paste0("transferring ",names(transfer.other.assays.layers)[i], " ", transfer.other.assays.layers[i]))
        as.data.frame(Matrix::t((Seurat::GetAssayData(object = object,
                                                      assay = names(transfer.other.assays.layers)[i],
                                                      slot = transfer.other.assays.layers[[i]]))))
      }
    )
    names(obsm.assays) <- paste0("assay_", names(transfer.other.assays.layers), "_", transfer.other.assays.layers)
    obsm <- c(obsm, obsm.assays)
  }

  obsp <- NULL
  if(transfer.graphs){
    graph.list <- names(object@graphs)
    message("getting sparse Seurat Graphs from R to python")
    obsp.graphs <- lapply(seq_along(graph.list), function(i){
      message(paste0("getting ", graph.list[i]))
      reticulate::r_to_py(as(SeuratObject::Graphs(object, graph.list[i]), "dgCMatrix"))
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

#' Run Harmony with necessary internal object tinkering for Symphony
#'
#' @return Runs Harmony in the implementation of Symphony such that the Harmony
#' model is stored within the Seurat object for later use with Symphony, see
#' https://github.com/immunogenomics/symphony/blob/main/vignettes/utils_seurat.R
#' @references copy-pasted from @immunogenomics/harmony github repo
#'
#' @importFrom rlang %||%
#'
#' @export
RunHarmonyForSymphony <- function(
    object,
    group.by.vars,
    reduction.use = 'pca',
    dims.use = NULL,
    reduction.save = "harmony",
    project.dim = TRUE,
    ...
) {

  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop("Running Harmony on a Seurat object requires Seurat")
  }
  if (!reduction.use %in% Seurat::Reductions(object = object)) {
    stop(paste(reduction.use, "cell embeddings not found in Seurat object.",
               "For a Seurat preprocessing walkthrough, please refer to the vignette"))
  }
  embedding <- Seurat::Embeddings(object, reduction = reduction.use)
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(embedding))
  }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rerun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  metavars_df <- Seurat::FetchData(
    object,
    group.by.vars,
    cells = Seurat::Cells(x = object[[reduction.use]])
  )

  #changes here
  ####
  data_mat <- embedding[,dims.use] #change #1
  data_mat <- Matrix::t(data_mat) #change #2 need to do this because
  # RunHarmony.default does this anyways and we want data_mat to be ready to be
  # used to add the row and col names
  harmonyObject <- harmony::RunHarmony(
    data_mat = data_mat,
    meta_data = metavars_df,
    vars_use = group.by.vars,
    return_object = TRUE, #this is the big change, #3
    ...
  )

  # #this is needed beacuse we need to extract the embeddings for the object
  result_matrix <- as.matrix(harmonyObject$Z_corr)
  row.names(result_matrix) <- row.names(data_mat)
  colnames(result_matrix) <- colnames(data_mat)
  harmonyEmbed <- t(result_matrix)
  # #end extraction

  ####

  reduction.key <- Seurat::Key(reduction.save, quiet = TRUE)
  rownames(harmonyEmbed) <- rownames(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.key, seq_len(ncol(harmonyEmbed)))

  object[[reduction.save]] <- Seurat::CreateDimReducObject(
    embeddings = harmonyEmbed,
    stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
    assay = Seurat::DefaultAssay(object = object[[reduction.use]]),
    key = reduction.key,
    misc = list(harmony_object = harmonyObject)#another change
  )
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


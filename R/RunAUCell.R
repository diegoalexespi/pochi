#' Run AUCell on a Seurat object
#'
#' @param object Seurat object
#' @param assay Assay to use for building rankings. Will use default assay if NULL.
#' @param slot Slot to use for building rankings
#' @param genesets A list of vectors of features for expression programs; each entry should be a vector of feature names. If the list of vectors is named, the resulting AUCell score for the expression program will be a row in the AUCell assay returned.
#' @param ranking.save If TRUE, will save a new Assay `rankings` into `object` with the rankings of each gene per cell.
#' @param ranking.key If not NULL, will pull Assay `rankings` from `object` rather than re-calculating rankings.
#' @param normAUC Whether to normalize the maximum possible AUC per geneset to 1
#' @param aucMaxRank Threshold to calculate the AUC (see details section below and details section of \code{\link[AUCell]{AUCell_calcAUC}})
#' @param verbose Boolean. TRUE to show progress messages, FALSE to hide progress messages
#'
#' @details (Copied verbatim from AUCell::AUCell_calcAUC) _In a simplified way, the AUC value represents the fraction of genes, within the top X genes in the ranking, that are included in the signature. The parameter 'aucMaxRank' allows to modify the number of genes (maximum ranking) that is used to perform this computation. By default, it is set to 5% of the total number of genes in the rankings. Common values may range from 1 to 20%._
#'
#' @return Returns a Seurat object with the AUCell results stored as a \code{\link{Assay}} object within the Seurat object
#' @seealso \code{\link[AUCell]{AUCell_buildRankings}}
#' @seealso \code{\link[AUCell]{AUCell_calcAUC}}
#' @references Aibar et al. (2017) SCENIC: single-cell regulatory network inference and clustering. Nature Methods. doi: 10.1038/nmeth.4463
#'
#' @importFrom rlang %||%
#'
#' @export
RunAUCell <- function(
  object,
  assay = NULL,
  slot = "counts",
  genesets,
  ranking.save = FALSE,
  ranking.key = NULL,
  normAUC = TRUE,
  aucMaxRank = 0.05,
  verbose = TRUE,
  auc_assay_name = "AUC",
  ...
) {
  #SeuratWrappers:::CheckPackage(package = 'AUCell', repository = "bioconductor")
  assay <- assay %||% Seurat::DefaultAssay(object = object)
  my_data <- Seurat::GetAssayData(object = object, assay = assay, slot = slot)
  aucMaxRank <- ceiling(aucMaxRank*nrow(my_data))

  #filter genesets to exclude genesets with >80% of its genes missing
  available_genes <- rownames(my_data)
  geneset_names <- names(genesets)
  genesets_filtered <- lapply(1:length(genesets), function(i){
    target_geneset <- genesets[[i]]
    filtered_geneset <- target_geneset[target_geneset %in% available_genes]
    missing_genes <- target_geneset[!(target_geneset %in% available_genes)]
    if(length(missing_genes)/length(target_geneset) >= .80){
      warning(paste0("The AUC for the geneset ", geneset_names[i], " was not calculated, as 80% or more of its genes are missing."))
      return(NULL)
    }
    return(filtered_geneset)
  })
  names(genesets_filtered) <- geneset_names
  genesets_filtered <- genesets_filtered[lengths(genesets_filtered) > 0]
  geneset_names <- names(genesets_filtered)
  genesets_length <- lengths(genesets_filtered)

  #find the maxAUC for each geneset
  max_aucs <- lapply(1:length(genesets_filtered), function(i){
    x_th <- 1:genesets_length[i]
    x_th <- x_th[x_th<aucMaxRank]
    y_th <- seq_along(x_th)
    inter_1 <- c(x_th, aucMaxRank)
    inter_2 <- diff(inter_1)
    inter_3 <- sum(inter_2 * y_th)
    return(inter_3)
    })


  #get the rankings for each cell
  if(is.null(ranking.key)){
    message("building ranking matrix using sparseMatrixStats")
    gene_names <- rownames(my_data)
    cell_names <- colnames(my_data)
    my_data <- sparseMatrixStats::colRanks(-my_data, preserveShape = TRUE, ties.method = "min")

    dimnames(my_data) <- list(gene_names, cell_names)
    gene_names <- rownames(my_data)
    cell_names <- colnames(my_data)

    message("breaking ties in ranking matrix")
    my_data <- matrixStats::colRanks(my_data, preserveShape = TRUE, ties.method = "random")
    my_data[my_data > aucMaxRank] <- 0
    message("converting to sparse matrix")
    my_dgc <- as(my_data, "dgCMatrix")
    dimnames(my_dgc) = list(gene_names,cell_names)
    if(ranking.save){
      object[["ranking"]] <- CreateAssayObject(counts = my_dgc)
    }
  } else {
    message("using pre-built ranking assay")
    my_dgc <- GetAssayData(object, slot = "counts", assay = "ranking")
    gene_names <- rownames(my_dgc)
    cell_names <- colnames(my_dgc)
  }


  message("calculating AUC for...")
  #calculate the AUC
  geneset_aucs <- lapply(1:length(genesets_filtered), function(i){
    message(geneset_names[[i]])
    my_dgc_subset <- my_dgc[genesets_filtered[[i]],]
    my_dgc_sorted <- apply(my_dgc_subset,2,sort)
    my_dgc_sorted_temp <- rbind(my_dgc_sorted, aucMaxRank)
    my_dgc_sorted_temp <- my_dgc_sorted_temp[-1,]
    my_dgc_sorted_diffs <- my_dgc_sorted_temp - my_dgc_sorted
    column_zeros <- colSums(my_dgc_sorted == 0)
    mds_ranks <- t(matrixStats::colRanks(my_dgc_sorted,
                                         ties.method = "max",
                                         preserveShape=TRUE))
    colrank_matrix <- t(mds_ranks - column_zeros)
    cell_aucs <- colSums(my_dgc_sorted_diffs * colrank_matrix)
    return(cell_aucs/max_aucs[[i]])
    }) %>% do.call(rbind, .)
  colnames(geneset_aucs) <- cell_names
  rownames(geneset_aucs) <- geneset_names
  object[[auc_assay_name]] <- CreateAssayObject(counts = geneset_aucs)
  return(object)
}








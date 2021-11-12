#' Extract cell proportions across groups in a Seurat object
#'
#' @param object Seurat object
#' @param group.by Which meta.data slot to use for grouping the cells
#' @param split.by Which meta.data slot to use for splitting the groups
#' @param replicate.by Which meta.data slot to use as split.by replicates
#' @details Returns the abundances of specific groups in the Seurat object
#' across a split.by variable, using replicate.by as replicates for each
#' split.by condition
#'
#' @return Returns a data frame of group.by abundances per replicate in the
#' split.by variable
#'
#' @importFrom rlang %||%
#'
#' @export
ExtractClusterProportions <- function(object,
                                      group.by,
                                      split.by,
                                      replicate.by){

  seurat_metadata <- object@meta.data
  seurat_metadata_filtered <- seurat_metadata %>%
    dplyr::group_by(!!sym(replicate.by), !!sym(split.by), !!sym(group.by)) %>%
    dplyr::tally() %>%
    dplyr::mutate(percent = n/sum(n)) %>%
    dplyr::ungroup() %>%
    droplevels() %>%
    tidyr::complete(nesting(!!sym(replicate.by), !!sym(split.by)), !!sym(group.by), fill = list(percent = 0, n = 0))

  split.by.values <- unique(seurat_metadata_filtered[[split.by]])

  frequencies <- lapply(seq_along(split.by.values), function(i){
    seurat_metadata_filtered %>%
      dplyr::filter(!!sym(split.by) == split.by.values[i]) %>%
      dplyr::ungroup() %>%
      dplyr::select(!!sym(replicate.by), !!sym(group.by), percent) %>%
      dplyr::mutate(split_by = split.by.values[i], replicate_by = !!sym(replicate.by), group_by = !!sym(group.by))
  }) %>% do.call(rbind, .)
  return(frequencies)
}

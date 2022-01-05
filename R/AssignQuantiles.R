#' Assign cells to feature-specific quantiles.
#'
#' @param object Seurat object
#' @param feature Which feature to use for quantile calculation
#' @param assay Which assay the feature is located within
#' @param slot Which slot the to pull from within specified assay
#' @param split.by If not NULL, which metadata feature to split by before
#' calculating quantiles. Will then calculate quantiles within each split object
#' @param quantile.probs Which quantiles to split at
#' @details Determines the quantile of expression for each cell
#' for a specified feature in the Seurat object
#'
#' @return Returns a Seurat object with quantile metadata added
#'
#' @importFrom rlang %||%
#'
#' @export
AssignQuantiles <- function(object,
                            feature,
                            assay = "ADT",
                            slot = "data",
                            split.by = "run_10x",
                            quantile.probs = c(0,.25,.5,.75,1)){


  DefaultAssay(object) <- assay
  if(length(feature) != 1){
    stop("feature must be length 1")
  }

  if(!is.null(split.by)){
    split_objects <- SplitObject(object, split.by = split.by)
  } else {
    split_objects <- list(object)
  }

  split_objects <- lapply(seq_along(split_objects), function(i){
    temp_object <- split_objects[[i]]
    temp_vector <- Seurat::FetchData(temp_object, vars = feature, slot = slot)[,1]
    temp_quantiles <- quantile(temp_vector, probs = quantile.probs)
    if(length(unique(temp_quantiles)) < length(temp_quantiles)){
      stop("quantiles are non-unique, set to different values")
    }
    temp_cuts <- cut(temp_vector, breaks = temp_quantiles, include.lowest = TRUE)
    temp_cuts_levels <- levels(temp_cuts)
    temp_cuts <- forcats::fct_recode(temp_cuts, !!!setNames(temp_cuts_levels, paste0(feature, ".q", 1:length(temp_cuts_levels))))
    temp_object[[paste0(feature, ".quantile")]] <- temp_cuts
    return(temp_object)
  })

  if(!is.null(split.by)){
    merged_object <- merge(split_objects[[1]], split_objects[2:length(split_objects)])
  } else {
    merged_object <- split_objects[[1]]
  }

  object[[paste0(feature, ".quantile")]] <- merged_object[[paste0(feature, ".quantile")]]
  return(object)
}

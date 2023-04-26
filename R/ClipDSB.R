#' Clip DSB assay with quantiles per feature
#'
#' @param object Seurat object
#' @param assay Which assay, typically DSB
#' @param slot Which slot in the assay, typically data after DSB correction
#' @param split.by If not NULL, which metadata feature to split by  before
#' scaling.
#' @param low.p Which percentile below which to clip per feature
#' @param high.p Which percentile above which to clip per feature
#' @param new.assay Name of new assay to save values to
#' @details Clips the expression values for each feature in the specified assay
#' from low.p to high.p, with an optional split.by argument.
#'
#' @return Returns a Seurat object with clipped values in the data slot of the
#' specified assay
#'
#' @importFrom rlang %||%
#'
#' @export
ClipDSB <- function(object, assay = "DSB", slot = "data", split.by = NULL,
                     low.p = NULL, high.p = 0.999, new.assay = "sDSB"){

  #split the object or not...
  if(!is.null(split.by)){
    split_object <- Seurat::SplitObject(object, split.by = split.by)
  } else {
    split_object <- list(object)
  }


  #get the assay from each object
  split_assays <- lapply(seq_along(split_object), function(i){
    object_assay <- SeuratObject::GetAssayData(split_object[[i]], assay = assay, slot = slot)
    features_to_scale <- rownames(object_assay)
    clipped_matrix <- lapply(seq_along(features_to_scale), function(i){
      value_vector <- object_assay[features_to_scale[i],]
      percentile_limit_hi <- quantile(value_vector, high.p)
      if(!is.null(low.p)){
        percentile_limit_lo <- quantile(value_vector, low.p)
        value_vector[value_vector < percentile_limit_lo] <- percentile_limit_lo
      }
      value_vector[value_vector > percentile_limit_hi] <- percentile_limit_hi
      value_matrix <- data.frame(feature = value_vector)
      colnames(value_matrix) <- features_to_scale[i]
      return(t(value_matrix))
    })
    clipped_matrix <- do.call(rbind, clipped_matrix)
    return(clipped_matrix)
  })

  #combine the split object assays back into one assay and keep original order
  #of cells to add back to object
  merged_assay <- do.call(cbind, split_assays)
  merged_assay <- merged_assay[,colnames(object)]

  #return object with slot updated
  object[[new.assay]] <- CreateAssayObject(data = merged_assay,
                                           min.cells = 0,
                                           min.features = 0)

  return(object)
}

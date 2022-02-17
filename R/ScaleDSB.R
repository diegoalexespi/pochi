#' Scale DSB assay to a global max per feature
#'
#' @param object Seurat object
#' @param assay Which assay, typically DSB
#' @param slot Which slot in the assay, typically data after DSB correction
#' @param split.by If not NULL, which metadata feature to split by  before
#' scaling.
#' @param high.p Which percentile to set as the global max per object
#' @param scale.to.one Whether to set global max per feature as 1 (helpful when
#' needing equally weighted features)
#' @param new.assay Name of new assay to save values to
#' @details Scales the expression values for each feature in the specified assay
#' from 0 to a global max, with an optional split.by argument. The 0 is assumed
#' to be the 0 values in the slot, in line with DSB's procedure. We only scale
#' the high percentile.
#'
#' @return Returns a Seurat object with scaled values in the data slot of the
#' specified assay
#'
#' @importFrom rlang %||%
#'
#' @export
ScaleDSB <- function(object, assay = "DSB", slot = "data", split.by = NULL,
                     high.p = 0.999, scale.to.one = TRUE, new.assay = "sDSB"){

  #split the object or not...
  if(!is.null(split.by)){
    split_object <- Seurat::SplitObject(object, split.by = split.by)
  } else {
    split_object <- list(object)
  }

  #get global maximum first
  if(scale.to.one){
    feature_global_maxs <- rep(1, length(rownames(object)))
    names(feature_global_maxs) <- rownames(object)
  } else {
    object_global_maxs_mat <- lapply(split_object, function(x){
      object_assay <- SeuratObject::GetAssayData(x, assay = assay, slot = slot)
      features_to_scale <- rownames(object_assay)
      global_maxs <- lapply(seq_along(features_to_scale), function(i){
        value_vector <- object_assay[features_to_scale[i],]
        percentile_limit_hi <- quantile(value_vector, high.p)
        return(percentile_limit_hi)
      }) %>% unlist()
      names(global_maxs) <- features_to_scale
      return(global_maxs)
    }) %>% do.call(rbind, .)
    feature_global_maxs <- matrixStats::colMaxs(object_global_maxs_mat)
    names(feature_global_maxs) <- colnames(object_global_maxs_mat)
  }


  #get the assay from each object
  split_assays <- lapply(seq_along(split_object), function(i){
    object_assay <- SeuratObject::GetAssayData(split_object[[i]], assay = assay, slot = slot)
    features_to_scale <- rownames(object_assay)
    scaled_matrix <- lapply(seq_along(features_to_scale), function(i){
      value_vector <- object_assay[features_to_scale[i],]
      value_vector_nonzero <- value_vector[value_vector > 0]
      nonzero_min <- min(value_vector_nonzero)
      percentile_limit_hi <- quantile(value_vector_nonzero, high.p)
      hi_scale_limit <- feature_global_maxs[features_to_scale[i]]
      value_vector_nonzero[value_vector_nonzero > percentile_limit_hi] <- percentile_limit_hi
      value_vector[value_vector > 0] <- scales::rescale(value_vector_nonzero, to = c(nonzero_min, hi_scale_limit))
      value_matrix <- data.frame(feature = value_vector)
      colnames(value_matrix) <- features_to_scale[i]
      return(t(value_matrix))
    })
    scaled_matrix <- do.call(rbind, scaled_matrix)
    return(scaled_matrix)
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

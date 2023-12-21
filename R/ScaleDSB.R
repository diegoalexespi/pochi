#' Scale DSB assay to [0-1] per feature
#'
#' @param object Seurat object
#' @param assay Which assay, typically DSB
#' @param slot Which slot in the assay, typically data after DSB correction
#' @param split.by If not NULL, which metadata feature to split by  before
#' scaling.
#' @param high.p Which percentile to set as the global 1 per feature
#' @param scale.to.one Whether to scale the high.p value to 1 for each feature.
#' If not, will scale to median of the high.p quantiles across the split.by
#' groups.
#' @param new.assay Name of new assay to save values to
#' @details Scales the expression values for each feature in the specified assay
#' from 0 to 1, with an optional split.by argument.
#'
#' @return Returns a Seurat object with scaled values in the data slot of the
#' specified assay
#'
#' @importFrom rlang %||%
#'
#' @export
ScaleDSB <- function(object,
                     assay = "DSB",
                     slot = "data",
                     split.by = NULL,
                     high.p = 0.999,
                     scale.to.one = FALSE,
                     new.assay = "sDSB"){

  #split the object or not...
  if(!is.null(split.by)){
    split_object <- Seurat::SplitObject(object, split.by = split.by)
  } else {
    split_object <- list(object)
  }

  if(!scale.to.one){
    #get the median DSB high quantile value from each feature across objects
    max_dsb_values <- lapply(seq_along(split_object), function(i){
      object_assay <- SeuratObject::GetAssayData(split_object[[i]], assay = assay, slot = slot)
      object_maxs <- matrixStats::rowQuantiles(as.matrix(object_assay), probs = high.p)
      return(data.frame(object_maxs))
    }) %>% do.call(cbind, .)
    feature_medians <- matrixStats::rowMedians(as.matrix(max_dsb_values))
  }

  #quantile-cap and scale the assay from each object
  split_assays <- lapply(seq_along(split_object), function(i){
    object_assay <- SeuratObject::GetAssayData(split_object[[i]], assay = assay, slot = slot)
    features_to_scale <- rownames(object_assay)
    scaled_matrix <- lapply(seq_along(features_to_scale), function(i){
      value_vector <- object_assay[features_to_scale[i],]
      value_vector[value_vector < 0] <- 0
      percentile_limit_hi <- quantile(value_vector, high.p)
      value_vector[value_vector > percentile_limit_hi] <- percentile_limit_hi
      if(scale.to.one){
        scaled_vector <- scales::rescale(value_vector, to = c(0, 1))
        value_matrix <- data.frame(feature = scaled_vector)
      } else {
        high_value <- feature_medians[features_to_scale[i]]
        scaled_vector <- scales::rescale(value_vector, to = c(0, high_value))
        value_matrix <- data.frame(feature = value_vector)
      }
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

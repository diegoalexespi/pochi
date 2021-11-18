#' Scale ADT assay from 0 to 1 per feature
#'
#' @param object Seurat object
#' @param assay Which assay, typically ADT
#' @param slot Which slot in the assay, typically data after CLR
#' @param split.by If not NULL, which metadata feature to split by  before
#' scaling.
#' @param low.p Which percentile to set as 0
#' @param high.p Which percentile to set as 1
#' @param new.assay Name of new assay to save values to
#' @details Scales the expression values for each feature in the specified assay
#' from 0 to 1, with an optional split.by argument
#'
#' @return Returns a Seurat object with scaled values in the data slot of the
#' specified assay
#'
#' @importFrom rlang %||%
#'
#' @export
ScaleADT <- function(object, assay = "ADT", slot = "data", split.by = NULL,
                     low.p = 0.001,high.p = 0.999,new.assay = "sADT"){

  #split the object or not...
  if(!is.null(split.by)){
    split_object <- Seurat::SplitObject(object, split.by = split.by)
  } else {
    split_object <- list(object)
  }

  #get the assay from each object
  split_assays <- lapply(split_object, function(x){
    object_assay <- SeuratObject::GetAssayData(x, assay = assay, slot = slot)
    features_to_scale <- rownames(object_assay)
    scaled_matrix <- lapply(seq_along(features_to_scale), function(i){
      value_vector <- object_assay[features_to_scale[i],]
      percentile_limits <- quantile(value_vector, c(low.p, high.p))
      value_vector[value_vector < percentile_limits[1] ] <- percentile_limits[1]
      value_vector[value_vector > percentile_limits[2] ] <- percentile_limits[2]
      value_vector <- scales::rescale(value_vector, to = c(0,1))
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

#' Scale ADT assay from 0 to 1 per feature
#'
#' @param object Seurat object
#' @param assay Which assay, typically ADT
#' @param layer Which layer in the assay, typically data after CLR
#' @param split.by If not NULL, which metadata feature to split by before
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
ScaleADT <- function(object,
                     assay = "ADT",
                     layer = "data",
                     split.by = NULL,
                     low.p = 0.001,
                     high.p = 0.999,
                     scale.to.one = FALSE,
                     new.assay = "sADT"){
    
    #split the object or not...
    if(!is.null(split.by)){
        split_object <- Seurat::SplitObject(object, split.by = split.by)
    } else {
        split_object <- list(object)
    }
    
    if(!scale.to.one){
        #get the median high quantile value from each feature across objects
        max_dsb_values <- lapply(seq_along(split_object), function(i){
            object_assay <- SeuratObject::LayerData(split_object[[i]][[assay]], layer = layer)
            object_maxs <- matrixStats::rowQuantiles(as.matrix(object_assay), probs = high.p)
            return(data.frame(object_maxs))
        }) %>% do.call(cbind, .)
        max_feature_medians <- matrixStats::rowMedians(as.matrix(max_dsb_values))
        
        #get the median low quantile value from each feature across objects
        min_dsb_values <- lapply(seq_along(split_object), function(i){
            object_assay <- SeuratObject::LayerData(split_object[[i]][[assay]], layer = layer)
            object_mins <- matrixStats::rowQuantiles(as.matrix(object_assay), probs = low.p)
            return(data.frame(object_mins))
        }) %>% do.call(cbind, .)
        min_feature_medians <- matrixStats::rowMedians(as.matrix(min_dsb_values))
    }
    
    #get the assay from each object
    split_assays <- lapply(split_object, function(x){
        object_assay <- SeuratObject::LayerData(x[[assay]], layer = layer)
        features_to_scale <- rownames(object_assay)
        scaled_matrix <- lapply(seq_along(features_to_scale), function(i){
            value_vector <- object_assay[features_to_scale[i],]
            percentile_limits <- quantile(value_vector, c(low.p, high.p))
            value_vector[value_vector < percentile_limits[1]] <- percentile_limits[1]
            value_vector[value_vector > percentile_limits[2]] <- percentile_limits[2]
            if(scale.to.one){
                scaled_vector <- scales::rescale(value_vector, to = c(0, 1))
                value_matrix <- data.frame(feature = scaled_vector)
            } else {
                low_value <- min_feature_medians[features_to_scale[i]]
                high_value <- max_feature_medians[features_to_scale[i]]
                scaled_vector <- scales::rescale(value_vector, to = c(low_value, high_value))
                value_matrix <- data.frame(feature = scaled_vector)
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
    object[[new.assay]] <- SeuratObject::CreateAssay5Object(counts = merged_assay,
                                                            data = merged_assay,
                                                            min.cells = 0,
                                                            min.features = 0)
    
    return(object)
}

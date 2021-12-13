#' Plot a metadata selection in a biaxial FeatureScatter plot
#'
#' @param object Seurat object
#' @param feature1 The x-axis feature
#' @param feature2 The y-axis feature
#' @param metadata.col The column within the object meta.data to select from
#' @param metadata.selection The entry within metadata.col to highlight
#' @param assay Which assay
#' @param slot Which slot in the assay, typically data
#' @param label.cols 2-d vector specifying colors for labelling
#' @param pt.size Size of geom_point
#' @param add.density Whether to add geom_density_2d estimation to figure
#' @param density.col The color of the density contour (if added)
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
BackgatePlot <- function(object,
                         feature1,
                         feature2,
                         metadata.col,
                         metadata.selection,
                         assay = "RNA",
                         slot = "data",
                         label.cols = c("red", "grey"),
                         pt.size = 1,
                         add.density = FALSE,
                         density.col = "black"){
  DefaultAssay(seurat_object) <- assay
  plotting_data <- FetchData(seurat_object, vars = c(feature1, feature2, metadata.col), slot = slot)
  plotting_data$selected_cells <- ifelse(plotting_data[[metadata.col]] == metadata.selection,
                                         metadata.selection,
                                         "non-selected")
  plotting_data <- plotting_data %>%
    dplyr::mutate(selected_cells = factor(selected_cells, levels = c(metadata.selection, "non-selected"))) %>%
    dplyr::arrange(desc(selected_cells))

  g <- ggplot(plotting_data, aes(x = !!sym(feature1), y = !!sym(feature2), color = selected_cells, group = selected_cells ))+
    geom_point(size = pt.size)+
    scale_color_manual(values = c(label.cols))+
    theme_classic()
  if(add.density){
    g <- g + geom_density_2d(data = plotting_data[plotting_data$selected_cells == metadata.selection,], color = density.col)
  }
  g

}

#' Plot Seurat metadata distribution as a bar graph for indicated groups.
#'
#' @param object Seurat object
#' @param group.by Which meta.data slot to use for assessment
#' @param split.by Which meta.data slot to use for group splitting
#' @param as.freq Whether to plot the group.by metadata as a frequency or
#' proportion
#' @param label.count Whether to add count on top of bar graph
#' @param angled.text Whether to angle the text to 45 degrees
#' @details Plots the group.by metadata as blocks on the y axis split across the
#' specified metadata in the split.by argument.
#'
#' @return Returns a metadata plot
#'
#' @importFrom rlang %||%
#'
#' @export
MetaDataPlot <- function(object,
                         group.by,
                         split.by,
                         as.freq = TRUE,
                         label.count = FALSE,
                         text.size = 3,
                         angled.text = TRUE,
                         label.angle = 45){

  seurat_metadata <- object@meta.data
  grouped_count <- seurat_metadata %>%
    dplyr::group_by(!!sym(split.by), !!sym(group.by)) %>%
    dplyr::tally() %>%
    dplyr::mutate(y = n)
  if(as.freq){
    grouped_count <- dplyr::mutate(grouped_count, y = y/sum(y))
  }
  g <- ggplot(grouped_count, aes(x = !!sym(split.by), y = y, fill = !!sym(group.by)))+
    geom_col(color = "black")+
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle  = 45, hjust = 1, vjust = 1))


  if(as.freq){
    g <- g + ylab("Frequency")
  } else {
    g <- g + ylab("Count")
    if(label.count){
      if(angled.text) {
        g <- g + geom_text(aes(label = stat(y), group = !!sym(split.by)),
                           stat = 'summary', fun = sum, vjust = 0.5,
                           size = text.size, angle = 90, hjust = -0.1)
      } else {
        g <- g + geom_text(aes(label = stat(y), group = !!sym(split.by)),
                           stat = 'summary', fun = sum, vjust = -1,
                           size = text.size)
      }
    }
  }
  g
}

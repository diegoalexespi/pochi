#' Calculate percentage of mitochondrial, hemoglobin, and ribosomal genes
#'
#' @param object Seurat object
#' @param mt.pattern Which pattern to select mitochondrial genes
#' @param mt.pattern Which pattern to select hemoglobin genes
#' @param ribo.pattern Which pattern to select ribosomal genes
#' @details Calculate percentage of mitochondrial, hemoglobin, and ribosomal 
#' genes, show plots and returns the Seurat object with information stored
#'
#' @return Returns a Seurat object containing percentage of mitochondrial, 
#' hemoglobin, and ribosomal genes
#'
#' @importFrom rlang %||%
#'
#' @export
Calc.Perc.Features <- function(object,
                               mt.pattern = "^MT-",
                               hb.pattern = "^HB[^(P)]",
                               ribo.pattern = ("RPS|RPL"), 
                               plot.name="") {
  object@misc$cell.recovered = ncol(object)
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = mt.pattern)
  object[["percent.hb"]] <- PercentageFeatureSet(object, pattern = hb.pattern)
  object[["percent.ribo"]] <- PercentageFeatureSet(object, pattern = ribo.pattern)
  plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(object, feature1 = "percent.ribo", feature2 = "percent.hb")
  plot4 <- VlnPlot(object, features = c("percent.ribo", "percent.hb", "percent.mt"), ncol = 3)
  plot <- ((plot1 + plot2) / (plot3 + plot4))
  plot = plot + plot_annotation(title = plot.name, theme = theme(plot.title = element_text(hjust = 0.5)))
  print(plot)
  return(object)
}

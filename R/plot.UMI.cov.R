#' Plot UMIs distribution/counts
#'
#' @param object Seurat object
#' @param split.by Which meta.data slot to use for splitting the plots
#' @details Plots the UMIs distribution/counts in the Seurat object
#' across a split.by variable
#'
#' @return Print a plot of separate ggplots per split.by variable
#'
#' @importFrom rlang %||%
#'
#' @export
plot.UMI.cov <- function(object, split.by=NULL) {
  
  if(is.null(split.by)==FALSE) {
	Idents(object) <- split.by
  }
  libs = unique(Idents(object))
  
  #Create an object for each metadata var
  for (i in 1:length(libs)) {
    lib.seurat = subset(object, idents = libs[i])

    colSums.counts = as.data.frame(colSums(lib.seurat@assays$RNA@counts))
    colSums.counts$colSums <- colSums.counts$`colSums(lib.seurat@assays$RNA@counts)`
    colSums.counts$`colSums(lib.seurat@assays$RNA@counts)` <- NULL
    a=ggplot(colSums.counts, aes(x=colSums))+ geom_histogram(color="black", fill="white", bins=50) + ggtitle(libs[[i]])
    b=ggplot(colSums.counts, aes(x=colSums))+ geom_histogram(color="black", fill="white", bins=50, aes(y=..density..))+geom_density(alpha=.2, fill="#FF6666") + ggtitle(libs[[i]])
    print(a+b)
  }
}

#' Plot Seurat VlnPlot showing median of population
#'
#' @param object Seurat object
#' @param features Which features to plot
#' @param median.colour Which colour set the median dot
#' @param median.pt.size Which size is median dot
#' @param cell.pt.size Which size is each cell/dot
#' @param ncol How many columns to organize the plot if more than 1 plot.
#' @param legend Whether to show the legend or not.
#' @details Plot Seurat VlnPlot showing median of population
#'
#' @return Returns a patchwork plot
#'
#' @importFrom rlang %||%
#'
#' @export
VlnPlot.median <- function(object, features, median.colour = "black", median.pt.size=1, cell.pt.size = 0, ncol = NULL, legend = TRUE) {
	#Tris function creates VlnPlot plus a median 
	#VlnPlot.median(object, c("IL1B","CLEC10A"))
	
	myplots <- vector("list")
	
	#Create a plot for each gene and add median
	for (gene in features) {
		if (legend == TRUE) {
		myplots[[gene]] <- Seurat::VlnPlot(object = object, features = gene, pt.size = cell.pt.size) +
		stat_summary(fun = median.stat, geom='point', size = median.pt.size, colour = median.colour)
		} else if (legend == FALSE) {
		myplots[[gene]] <- Seurat::VlnPlot(object = object, features = gene, pt.size = cell.pt.size) +
		stat_summary(fun = median.stat, geom='point', size = median.pt.size, colour = median.colour) + NoLegend()
		}
	}
	#patchwork function to combine multiple plots
	patchwork::wrap_plots(myplots, ncol=ncol)
	}
	
#This is part of VlnPlot.median function
median.stat <- function(x){
	out <- quantile(x, probs = c(0.5))
	names(out) <- c("ymed")
	return(out) 
	}

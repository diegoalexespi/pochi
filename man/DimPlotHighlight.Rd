% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SeuratPlotEnhancements.R
\name{DimPlotHighlight}
\alias{DimPlotHighlight}
\title{Seurat DimPlot enhancement III}
\usage{
DimPlotHighlight(
  object,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt.size = NULL,
  reduction = NULL,
  group.by = NULL,
  split.by = NULL,
  shape.by = NULL,
  order = NULL,
  shuffle = FALSE,
  seed = 1,
  label = FALSE,
  label.size = 4,
  label.color = "black",
  label.box = FALSE,
  repel = FALSE,
  cols.highlight = "#DE2D26",
  sizes.highlight = 1,
  na.value = "grey50",
  ncol = NULL,
  combine = TRUE,
  raster = NULL,
  metadata.col,
  metadata.selection
)
}
\arguments{
\item{object}{Seurat object}

\item{dims}{Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions}

\item{cells}{Vector of cells to plot (default is all cells)}

\item{cols}{Vector of colors, each color corresponds to an identity class. This may also be a single character
or numeric value corresponding to a palette as specified by \code{\link[RColorBrewer]{brewer.pal.info}}.
By default, ggplot2 assigns colors. We also include a number of palettes from the pals package.
See \code{\link[Seurat]{DiscretePalette}} for details.}

\item{pt.size}{Adjust point size for plotting}

\item{reduction}{Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca}

\item{group.by}{Name of one or more metadata columns to group (color) cells by
(for example, orig.ident); pass 'ident' to group by identity class}

\item{split.by}{A factor in object metadata to split the plot by, pass 'ident'
to split by cell identity}

\item{shape.by}{If NULL, all points are circles (default). You can specify any
cell attribute (that can be pulled with FetchData) allowing for both
different colors and different shapes on cells.  Only applicable if \code{raster = FALSE}.}

\item{order}{Specify the order of plotting for the idents. This can be
useful for crowded plots if points of interest are being buried. Provide
either a full list of valid idents or a subset to be plotted last (on top)}

\item{shuffle}{Whether to randomly shuffle the order of points. This can be
useful for crowded plots if points of interest are being buried. (default is FALSE)}

\item{seed}{Sets the seed if randomly shuffling the order of points.}

\item{label}{Whether to label the clusters}

\item{label.size}{Sets size of labels}

\item{label.color}{Sets the color of the label text}

\item{label.box}{Whether to put a box around the label text (geom_text vs
geom_label)}

\item{repel}{Repel labels}

\item{cols.highlight}{A vector of colors to highlight the cells as; will
repeat to the length groups in cells.highlight}

\item{sizes.highlight}{Size of highlighted cells; will repeat to the length
groups in cells.highlight.  If \code{sizes.highlight = TRUE} size of all
points will be this value.}

\item{na.value}{Color value for NA points when using custom scale}

\item{ncol}{Number of columns for display when combining plots}

\item{combine}{Combine plots into a single \code{\link[patchwork]{patchwork}ed}
ggplot object. If \code{FALSE}, return a list of ggplot objects}

\item{raster}{Convert points to raster format, default is \code{NULL} which
automatically rasterizes if plotting more than 100,000 cells}

\item{metadata.col}{The column within the object meta.data to select from}

\item{metadata.selection}{The entry within metadata.col to highlight in the
Seurat DimPlot}
}
\value{
Returns a DimPlot
}
\description{
Seurat DimPlot enhancement III
}
\details{
Identical to DimPlot with the additional step of highlighting a
specific piece of given metadata
}

# Seurat DimPlot enhancement III

Seurat DimPlot enhancement III

## Usage

``` r
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
```

## Arguments

- object:

  Seurat object

- dims:

  Dimensions to plot, must be a two-length numeric vector specifying x-
  and y-dimensions

- cells:

  Vector of cells to plot (default is all cells)

- cols:

  Vector of colors, each color corresponds to an identity class. This
  may also be a single character or numeric value corresponding to a
  palette as specified by
  [`brewer.pal.info`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html).
  By default, ggplot2 assigns colors. We also include a number of
  palettes from the pals package. See
  [`DiscretePalette`](https://satijalab.org/seurat/reference/DiscretePalette.html)
  for details.

- pt.size:

  Adjust point size for plotting

- reduction:

  Which dimensionality reduction to use. If not specified, first
  searches for umap, then tsne, then pca

- group.by:

  Name of one or more metadata columns to group (color) cells by (for
  example, orig.ident); pass 'ident' to group by identity class

- split.by:

  A factor in object metadata to split the plot by, pass 'ident' to
  split by cell identity

- shape.by:

  If NULL, all points are circles (default). You can specify any cell
  attribute (that can be pulled with FetchData) allowing for both
  different colors and different shapes on cells. Only applicable if
  `raster = FALSE`.

- order:

  Specify the order of plotting for the idents. This can be useful for
  crowded plots if points of interest are being buried. Provide either a
  full list of valid idents or a subset to be plotted last (on top)

- shuffle:

  Whether to randomly shuffle the order of points. This can be useful
  for crowded plots if points of interest are being buried. (default is
  FALSE)

- seed:

  Sets the seed if randomly shuffling the order of points.

- label:

  Whether to label the clusters

- label.size:

  Sets size of labels

- label.color:

  Sets the color of the label text

- label.box:

  Whether to put a box around the label text (geom_text vs geom_label)

- repel:

  Repel labels

- cols.highlight:

  A vector of colors to highlight the cells as; will repeat to the
  length groups in cells.highlight

- sizes.highlight:

  Size of highlighted cells; will repeat to the length groups in
  cells.highlight. If `sizes.highlight = TRUE` size of all points will
  be this value.

- na.value:

  Color value for NA points when using custom scale

- ncol:

  Number of columns for display when combining plots

- combine:

  Combine plots into a single
  [`patchwork`](https://patchwork.data-imaginist.com/reference/patchwork-package.html)`ed`
  ggplot object. If `FALSE`, return a list of ggplot objects

- raster:

  Convert points to raster format, default is `NULL` which automatically
  rasterizes if plotting more than 100,000 cells

- metadata.col:

  The column within the object meta.data to select from

- metadata.selection:

  The entry within metadata.col to highlight in the Seurat DimPlot

## Value

Returns a DimPlot

## Details

Identical to DimPlot with the additional step of highlighting a specific
piece of given metadata

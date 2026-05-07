# Heatmap of average expression values per group

Heatmap of average expression values per group

## Usage

``` r
DoClusteredHeatmap(
  object,
  features = NULL,
  assay = "RNA",
  layer = "data",
  group.by = "seurat_clusters",
  max_zed = 3,
  y_text_size = 10,
  scale_rows = TRUE,
  plot_dendro = FALSE,
  dot.scale = 5,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  viridis_option = "B",
  viridis_direction = 1
)
```

## Arguments

- object:

  Seurat object

- features:

  Which features to plot

- assay:

  Which assay to get expression values from

- layer:

  Which layer to pull the assay data from

- group.by:

  Which metadata layer to group cells by

- max_zed:

  Cutoff for maximum z-scores (absolute value)

- y_text_size:

  Size of text on y axis (features)

- scale_rows:

  Whether to scale the expression per row

- plot_dendro:

  Whether to plot accompanying x-axis and y-axis dendrograms

- dot.scale:

  Value by which to scale the dots in DotPlot

- cluster_cols:

  Whether to cluster the columns

- cluster_rows:

  Whether to cluster the rows

- viridis_option:

  The viridis color scale to use

- viridis_direction:

  The direction of the viridis color scale

## Value

Returns a patchwork object of a heatmap and dendrogram if specified

## Details

Plots the average expression value of the specified genes/features for
each group in the Seurat object, optionally using hierarchical
clustering for prettier visualization

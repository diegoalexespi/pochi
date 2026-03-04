# Seurat DotPlot enhancement with hierarchical clustering

Seurat DotPlot enhancement with hierarchical clustering

## Usage

``` r
DoClusteredDotPlot(
  object,
  features = NULL,
  assay = "RNA",
  slot = "data",
  group.by = "seurat_clusters",
  max_zed = 3,
  y_text_size = 10,
  scale_rows = TRUE,
  plot_dendro = FALSE,
  auc_choice = 0.6,
  dot.scale = 5,
  cluster_cols = TRUE
)
```

## Arguments

- object:

  Seurat object

- features:

  Which features to plot

- assay:

  Which assay to get expression values from

- slot:

  Which slot to pull the assay data from

- group.by:

  Which metadata slot to group cells by

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

  Whether to cluster the columns as well as the rows

## Value

Returns a patchwork object of a DotPlot and dendrogram if specified

## Details

Identical to DotPlot with the exception that the genes/features are
clustered using hierarchical clustering for prettier visualization

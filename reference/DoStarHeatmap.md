# Heatmap of average expression values per group with stars

Heatmap of average expression values per group with stars

## Usage

``` r
DoStarHeatmap(
  object,
  diff_exp_results = NULL,
  assay = "ADT",
  layer = "data",
  group.by = "seurat_clusters",
  max_zed = 3,
  plot_rownames = FALSE,
  label_features = TRUE,
  y_text_size = 10,
  scale_rows = TRUE,
  p_val_choice = 0.01,
  logFC_choice = 0.4,
  auc_choice = NULL,
  plot_all = FALSE,
  plot_dendro = FALSE,
  subset_features = NULL,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  viridis_option = "C",
  viridis_direction = 1,
  star_size = 8
)
```

## Arguments

- object:

  Seurat object

- viridis_option:

  The viridis color scale to use

- viridis_direction:

  The direction of the viridis color scale

## Value

Returns a patchwork object of a heatmap and dendrogram if specified

## Details

TBD!!!

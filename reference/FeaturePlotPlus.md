# Seurat FeaturePlot enhancement

Seurat FeaturePlot enhancement

## Usage

``` r
FeaturePlotPlus(
  object,
  feature,
  min.value = 0,
  assay = "RNA",
  viridis_choice = "D",
  ...
)
```

## Arguments

- object:

  Seurat object

- feature:

  Which feature to plot (singular)

- min.value:

  Which threshold value to set as not-expressed/grey in the FeaturePlot

- assay:

  Which assay the feature is located in

- viridis_choice:

  Which viridis palette to choose

- ...:

  Other arguments to Seurat::FeaturePlot

## Value

Returns a FeaturePlot

## Details

Identical to FeaturePlot with the exception that 0-values are set to
grey, which replicates the Monocle-type UMAP plots.

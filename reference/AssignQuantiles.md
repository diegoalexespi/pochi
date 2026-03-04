# Assign cells to feature-specific quantiles.

Assign cells to feature-specific quantiles.

## Usage

``` r
AssignQuantiles(
  object,
  feature,
  assay = "ADT",
  slot = "data",
  split.by = "run_10x",
  quantile.probs = c(0, 0.25, 0.5, 0.75, 1)
)
```

## Arguments

- object:

  Seurat object

- feature:

  Which feature to use for quantile calculation

- assay:

  Which assay the feature is located within

- slot:

  Which slot the to pull from within specified assay

- split.by:

  If not NULL, which metadata feature to split by before calculating
  quantiles. Will then calculate quantiles within each split object

- quantile.probs:

  Which quantiles to split at

## Value

Returns a Seurat object with quantile metadata added

## Details

Determines the quantile of expression for each cell for a specified
feature in the Seurat object

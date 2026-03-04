# Derive an SNN from a kNN Graph object in Seurat

Derive an SNN from a kNN Graph object in Seurat

## Usage

``` r
BuildSNNfromKNN(
  object,
  knn_choice = "RNA_nn",
  snn_name = paste0(knn_choice, "_snn"),
  prune_snn = 1/15
)
```

## Arguments

- object:

  Seurat object

- knn_choice:

  Name of Graph slot to use for SNN calculation

- snn_name:

  Name of Graph slot to store SNN

- prune_snn:

  Threshold at which to prune SNN (set to 0 below this)

## Value

Returns a Seurat object with the SNN Graph stored.

## Details

code is adapted from @immunogenomics/singlecellmethods. Builds a shared
nearest-neighbors graph from an input kNN using matrix multiplication so
pretty quick.

## References

code is adapted from @immunogenomics/singlecellmethods

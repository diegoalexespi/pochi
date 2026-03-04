# Average feature expression across clustered samples in a Seurat object using fast sparse matrix methods

Average feature expression across clustered samples in a Seurat object
using fast sparse matrix methods

## Usage

``` r
TrueAverageExpression(
  object,
  group.by = "seurat_clusters",
  group.by.delim = "_",
  assay = "RNA",
  slot = "data",
  verbose = TRUE
)
```

## Arguments

- object:

  Seurat object

- group.by:

  Ident with sample clustering information (default is seurat_clusters).
  Can also input a vector with multiple idents specified and
  TrueAverageExpression will return specified idents concatenated into
  new a new output ident with group.by.delim as the separator.

- group.by.delim:

  The delimiter used if group.by is a vector of length \> 1.

- assay:

  Assay to average (default is the active assay)

- slot:

  Slot to average (default is counts)

- verbose:

  Boolean or integer, show progress bar (default is TRUE)

## Value

dense matrix of cluster averages

## References

credit to GitHub user @zdebruine on Seurat issues!

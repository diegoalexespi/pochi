# Extract cell proportions and counts across groups in a Seurat object

Extract cell proportions and counts across groups in a Seurat object

## Usage

``` r
ExtractClusterProportions(object, group.by, split.by, replicate.by)
```

## Arguments

- object:

  Seurat object

- group.by:

  Which meta.data slot to use for grouping the cells

- split.by:

  Which meta.data slot to use for splitting the groups

- replicate.by:

  Which meta.data slot to use as split.by replicates

## Value

Returns a data frame of group.by abundances per replicate in the
split.by variable

## Details

Returns the abundances of specific groups in the Seurat object across a
split.by variable, using replicate.by as replicates for each split.by
condition

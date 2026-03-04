# Perform Shannon Component Analysis on a Seurat object

Perform Shannon Component Analysis on a Seurat object

## Usage

``` r
RunSCA(
  object,
  assay = "RNA",
  slot = "data",
  n.comps = 50L,
  iters = 3L,
  reduction.name = "sca",
  reduction.key = "sca_"
)
```

## Arguments

- object:

  Seurat object

- assay:

  Which assay to use for SCA calculation

- slot:

  Which slot in the assay to use for SCA calculation

- n.comps:

  Number of SCs to calculate

- iters:

  How many iterations to perform for SCA

- reduction.name:

  Name of the Reduction object to store the SCA as

- reduction.key:

  Key preceding the numbers of SCA components in reduction.name

## Value

Returns a Seurat object with an SCA calculated

## References

Demeo and Berger, \_biorxiv\_ 2021. shannonca.readthedocs.io

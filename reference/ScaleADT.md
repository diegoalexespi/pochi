# Scale ADT assay from 0 to 1 per feature

Scale ADT assay from 0 to 1 per feature

## Usage

``` r
ScaleADT(
  object,
  assay = "ADT",
  layer = "data",
  split.by = NULL,
  low.p = 0.001,
  high.p = 0.999,
  scale.to.one = FALSE,
  new.assay = "sADT"
)
```

## Arguments

- object:

  Seurat object

- assay:

  Which assay, typically ADT

- layer:

  Which layer in the assay, typically data after CLR

- split.by:

  If not NULL, which metadata feature to split by before scaling.

- low.p:

  Which percentile to set as 0

- high.p:

  Which percentile to set as 1

- new.assay:

  Name of new assay to save values to

## Value

Returns a Seurat object with scaled values in the data slot of the
specified assay

## Details

Scales the expression values for each feature in the specified assay
from 0 to 1, with an optional split.by argument

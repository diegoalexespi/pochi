# Scale DSB assay to \[0-1\] per feature

Scale DSB assay to \[0-1\] per feature

## Usage

``` r
ScaleDSB(
  object,
  assay = "DSB",
  layer = "data",
  split.by = NULL,
  high.p = 0.999,
  scale.to.one = FALSE,
  new.assay = "sDSB"
)
```

## Arguments

- object:

  Seurat object

- assay:

  Which assay, typically DSB

- layer:

  Which layer in the assay, typically data after DSB correction

- split.by:

  If not NULL, which metadata feature to split by before scaling.

- high.p:

  Which percentile to set as the global 1 per feature

- scale.to.one:

  Whether to scale the high.p value to 1 for each feature. If not, will
  scale to median of the high.p quantiles across the split.by groups.

- new.assay:

  Name of new assay to save values to

## Value

Returns a Seurat object with scaled values in the data slot of the
specified assay

## Details

Scales the expression values for each feature in the specified assay
from 0 to 1, with an optional split.by argument.

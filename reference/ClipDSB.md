# Clip DSB assay with quantiles per feature

Clip DSB assay with quantiles per feature

## Usage

``` r
ClipDSB(
  object,
  assay = "DSB",
  slot = "data",
  split.by = NULL,
  low.p = NULL,
  high.p = 0.999,
  new.assay = "sDSB"
)
```

## Arguments

- object:

  Seurat object

- assay:

  Which assay, typically DSB

- slot:

  Which slot in the assay, typically data after DSB correction

- split.by:

  If not NULL, which metadata feature to split by before scaling.

- low.p:

  Which percentile below which to clip per feature

- high.p:

  Which percentile above which to clip per feature

- new.assay:

  Name of new assay to save values to

## Value

Returns a Seurat object with clipped values in the data slot of the
specified assay

## Details

Clips the expression values for each feature in the specified assay from
low.p to high.p, with an optional split.by argument.

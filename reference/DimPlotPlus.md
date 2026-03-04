# Seurat DimPlot enhancement II

Seurat DimPlot enhancement II

## Usage

``` r
DimPlotPlus(
  object,
  reduction = "umap",
  shape.choice = 21,
  pt.size = 1.75,
  st.size = 0.25,
  ...
)
```

## Arguments

- object:

  Seurat object

- reduction:

  Which reduction to use in the DimPlot

- shape.choice:

  Which shape to use for the points int he DimPlot

- pt.size:

  Size of DimPlot points

- st.size:

  Size of DimPlot stencil around points

- ...:

  Other arguments to Seurat::DimPlot

## Value

Returns a DimPlot

## Details

Identical to DimPlot with the exception that points are drawn as
fill-able circles with a black border.

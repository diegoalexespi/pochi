# Seurat DimPlot enhancement

Seurat DimPlot enhancement

## Usage

``` r
DimPlotEdges(
  object,
  graph.choice = "RNA_nn",
  reduction = "umap",
  line.thickness = 0.05,
  line.color = "grey",
  raster = FALSE,
  ...
)
```

## Arguments

- object:

  Seurat object

- graph.choice:

  Which Graph to utilize for plotting with the DimPlot

- reduction:

  Which reduction to use in the DimPlot

- line.thickness:

  Thickness of lines connecting individual points in the DimPlot based
  on the graph.choice Graph

- line.color:

  Color of the lines connecting individual points

- raster:

  Raseterization boolean

- ...:

  Other arguments to Seurat::DimPlot

## Value

Returns a DimPlot

## Details

Identical to DimPlot with the exception that edges can be drawn on the
DimPlot between cells connected in the specified Graph.

# Plot a metadata selection in a biaxial FeatureScatter plot

Plot a metadata selection in a biaxial FeatureScatter plot

## Usage

``` r
FlowPlot(
  object,
  feature1,
  feature2,
  assay = "ADT",
  layer = "data",
  pt.size = 1,
  noise.zero = FALSE,
  noise.zero.bound = -0.01,
  gg_bandwidth = 0.1,
  rasterize = FALSE
)
```

## Arguments

- object:

  Seurat object

- feature1:

  The x-axis feature

- feature2:

  The y-axis feature

- assay:

  Which assay

- layer:

  Which layer in the assay, typically data

- pt.size:

  Size of geom_point

- noise.zero:

  Whether to simulate values from 0 to a negative value, done for visual
  purposes ONLY and for kernel density estimation when overlaying a
  density visualization. Will simulate values using a normal
  distribution centered at 0, take absolute values, and then scale from
  noise.zero.bound to 0.

- noise.zero.bound:

  The negative value to which to simulate zeros to.

- rasterize:

  Whether to rasterize the resulting plot.

## Value

Returns a biaxial plot utilizing ggpointdensity::geom_pointdensity

## Details

Scales the expression values for each feature in the specified assay
from 0 to 1, with an optional split.by argument

# Plot a metadata selection in a biaxial FeatureScatter plot

Plot a metadata selection in a biaxial FeatureScatter plot

## Usage

``` r
BackgatePlot(
  object,
  feature1,
  feature2,
  metadata.col,
  metadata.selection,
  assay = "RNA",
  slot = "data",
  label.cols = c("red", "grey"),
  pt.size = 1,
  add.density = FALSE,
  density.col = "black",
  noise.zero = FALSE,
  noise.zero.bound = -0.01,
  raster = TRUE,
  raster.dpi = c(512, 512)
)
```

## Arguments

- object:

  Seurat object

- feature1:

  The x-axis feature

- feature2:

  The y-axis feature

- metadata.col:

  The column within the object meta.data to select from

- metadata.selection:

  The entry within metadata.col to highlight

- assay:

  Which assay

- slot:

  Which slot in the assay, typically data

- label.cols:

  2-d vector specifying colors for labelling

- pt.size:

  Size of geom_point

- add.density:

  Whether to add geom_density_2d estimation to figure

- density.col:

  The color of the density contour (if added)

- noise.zero:

  Whether to simulate values from 0 to a negative value, done for visual
  purposes ONLY and for kernel density estimation when overlaying a
  density visualization. Will simulate values using a normal
  distribution centered at 0, take absolute values, and then scale from
  noise.zero.bound to 0.

- noise.zero.bound:

  The negative value to which to simulate zeros to.

- high.p:

  Which percentile to set as 1

## Value

Returns a Seurat object with scaled values in the data slot of the
specified assay

## Details

Scales the expression values for each feature in the specified assay
from 0 to 1, with an optional split.by argument

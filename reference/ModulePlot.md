# Plot avereage feature expression per sample across groups in a Seurat object

Plot avereage feature expression per sample across groups in a Seurat
object

## Usage

``` r
ModulePlot(
  seurat_object,
  assay = "AUC",
  features = rownames(seurat_object[[assay]]),
  split.by,
  replicate.by,
  draw_paths = FALSE,
  perform.stat.test = TRUE,
  test.choice = "t.test",
  ncols = 3,
  paired = FALSE,
  set_minimum_to_zero = FALSE,
  fill_colors = NULL,
  point_colors = NULL,
  p.adjust.method = "BH",
  filter.p.above = 1,
  step.increase = 0.15,
  title_size = 10,
  p_val_size = 2.5,
  spacing_above = 1.3,
  plot_type = "box",
  sina_shift = TRUE,
  x_lab = NULL,
  y_lab = NULL,
  same_y_limit = FALSE
)
```

## Arguments

- assay:

  Location of features (typically a gene module assay, not RNA...)

- features:

  Which features to use for grouping cells

- split.by:

  Which meta.data slot to use for splitting the groups

- replicate.by:

  Which meta.data slot to use as split.by replicates

- draw_paths:

  Whether or not to draw paths in the abundance plots between replicates
  if they exist

- perform.stat.test:

  Whether or not to perform a statistical test between the splits

- test.choice:

  Which choice of test to use (supports t.test of wilcox)

- ncols:

  Number of columns to plot the groups as

- paired:

  Whether the t.test should be paired if using t.test

- set_minimum_to_zero:

  Whether to set minimum abundance on plots to 0 globally

- point_colors:

  Vector of colors to color points as

- p.adjust.method:

  Which method to use for p-value adjustment for mulitple testing

- filter.p.above:

  Plot p-values only below this filter

- step.increase:

  Tweaking the space between brackets for p-values vertically

- title_size:

  Text size of each group plot's title

- p_val_size:

  Text size for p-values

- plot_type:

  Either box or violin plot

- sina_shift:

  Whether to use ggforce::geom_sina to pretty-shift the points

- x_lab:

  The label on the x axis

- y_lab:

  The label on the y axis

- same_y_limit:

  Whether to keep the upper y limit constant across all groups.

- object:

  Seurat object

- point_size:

  Size of the points plotted

## Value

Returns a patchwork plot of separate ggplots per features specified

## Details

Plots the average expression of specific features in the Seurat object
across a split.by variable, using replicate.by as replicates for each
split.by condition. Typically used for gene module testing. Not
recommended for RNA or single-cell expression levels that are not
Gaussian.

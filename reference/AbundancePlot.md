# Plot cell abundances across groups in a Seurat object

Plot cell abundances across groups in a Seurat object

## Usage

``` r
AbundancePlot(
  object,
  group.by,
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
  plot_type = "box",
  bees = TRUE,
  x_lab = NULL,
  y_lab = NULL,
  point_size = 1,
  same_y_limit = FALSE,
  selected.groups = NULL,
  rotated_axis = FALSE,
  print_results = FALSE
)
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

- x_lab:

  The label on the x axis

- y_lab:

  The label on the y axis

- point_size:

  Size of the points plotted

- same_y_limit:

  Whether to keep the upper y limit constant across all groups.

- selected.groups:

  Optional. If not NULL, will only plot the groups specified in this
  argument.

- rotated_axis:

  Whether to apply Seurat::RotatedAxis to the plots

- print_results:

  Prints to screen the results of the abundance calculations and
  abundance statistical testing.

- sina_shift:

  Whether to use ggforce::geom_sina to pretty-shift the points

## Value

Returns a patchwork plot of separate ggplots per group.by variable

## Details

Plots the abundances of specific groups in the Seurat object across a
split.by variable, using replicate.by as replicates for each split.by
condition

# Plot Seurat metadata distribution as a bar graph for indicated groups.

Plot Seurat metadata distribution as a bar graph for indicated groups.

## Usage

``` r
MetaDataPlot(
  object,
  group.by,
  split.by,
  as.freq = TRUE,
  label.count = FALSE,
  text.size = 3,
  angled.text = TRUE,
  label.angle = 45
)
```

## Arguments

- object:

  Seurat object

- group.by:

  Which meta.data slot to use for assessment

- split.by:

  Which meta.data slot to use for group splitting

- as.freq:

  Whether to plot the group.by metadata as a frequency or proportion

- label.count:

  Whether to add count on top of bar graph

- angled.text:

  Whether to angle the text to 45 degrees

## Value

Returns a metadata plot

## Details

Plots the group.by metadata as blocks on the y axis split across the
specified metadata in the split.by argument.

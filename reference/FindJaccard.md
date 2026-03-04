# Determine overlap measure between a list of genes

Determine overlap measure between a list of genes

## Usage

``` r
FindJaccard(
  genelist,
  distance.measure = "jaccard",
  as.distance = FALSE,
  return.tidy = TRUE
)
```

## Arguments

- genelist:

  Named list of gene vectors

- distance.measure:

  Which sort of distance measure to use (Jaccard by default)

- as.distance:

  Whether to transform similarities to distances

- return.tidy:

  Whether to return a tidy dataframe

## Value

Returns a matrix of genelist distance measures

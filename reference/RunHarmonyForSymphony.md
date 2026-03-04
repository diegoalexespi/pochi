# Run Harmony with necessary internal object tinkering for Symphony

Run Harmony with necessary internal object tinkering for Symphony

## Usage

``` r
RunHarmonyForSymphony(
  object,
  group.by.vars,
  reduction.use = "pca",
  dims.use = NULL,
  reduction.save = "harmony",
  project.dim = TRUE,
  ...
)
```

## Value

Runs Harmony in the implementation of Symphony such that the Harmony
model is stored within the Seurat object for later use with Symphony,
see
https://github.com/immunogenomics/symphony/blob/main/vignettes/utils_seurat.R

## References

copy-pasted from @immunogenomics/harmony github repo

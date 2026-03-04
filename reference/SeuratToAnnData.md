# Save a Seurat object as an anndata h5ad file

Save a Seurat object as an anndata h5ad file

## Usage

``` r
SeuratToAnnData(
  object,
  outFile = NULL,
  main.assay = "RNA",
  main.layer = "counts",
  transfer.other.main.layers = NULL,
  transfer.other.assays.layers = NULL,
  transfer.graphs = TRUE
)
```

## Arguments

- object:

  Seurat object

- outFile:

  Name of the target file to write the anndata object

- main.assay:

  Assay to transfer to AnnData as X

- main.layer:

  Layer to transfer to AnnData as X

- transfer.other.main.layers:

  Which other layers to transfer into AnnData

- transfer.other.assays.layers:

  Which other assay-layer pairs to transfer into AnnData (since AnnData
  only takes one assay for now)

- transfer.graphs:

  Whether to transfer Graphs to AnnData .uns

## Value

Writes a Seurat object to an .h5ad file after anndata conversion

## References

Code adapted from the awesome sceasy github!

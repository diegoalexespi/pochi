# Using pochi with Seurat

## Introduction

### Installation

Currently, [pochi](http://github.com/diegoalexespi/pochi) is available
at Github and can be downloaded using the devtools (or remotes) package.

``` r

devtools::install_github("diegoalexespi/pochi")
```

## Loading data

### Loading required packages

The [pochi](http://github.com/diegoalexespi/pochi) package operates on
[Seurat](http://github.com/satijalab/Seurat) objects. We load the
`Seurat`, `SeuratData`, and `pochi` packages here for our analyses, as
well as the `magrittr` package in order to improve legibility of code
through using the pipe `%>%` operator.

``` r

require("Seurat")
require("SeuratData")
require("pochi")
require("ggplot2")
require("tidyr")
```

We load the pbmc3k data from the `SeuratData` package.

``` r

InstallData("pbmc3k")
```

    ## Error in `InstallData()`:
    ## ! could not find function "InstallData"

``` r

pbmc3k_seurat <- LoadData("pbmc3k")
```

    ## Error in `LoadData()`:
    ## ! could not find function "LoadData"

``` r

set.seed(789)
pbmc3k_seurat <- pbmc3k_seurat[,!is.na(pbmc3k_seurat$seurat_annotations)]
```

    ## Error:
    ## ! object 'pbmc3k_seurat' not found

``` r

pbmc3k_seurat <- NormalizeData(pbmc3k_seurat)
```

    ## Error:
    ## ! object 'pbmc3k_seurat' not found

``` r

pbmc3k_seurat$condition <- sample(c("WT", "KO"), size = ncol(pbmc3k_seurat), replace = TRUE)
```

    ## Error:
    ## ! object 'pbmc3k_seurat' not found

``` r

pbmc3k_seurat$replicate <- sample(1:3, size = ncol(pbmc3k_seurat), replace = TRUE)
```

    ## Error:
    ## ! object 'pbmc3k_seurat' not found

``` r

pbmc3k_seurat@meta.data %>% head()
```

    ## Error:
    ## ! object 'pbmc3k_seurat' not found

## Visualizations

### `AbundancePlot`

``` r

AbundancePlot(pbmc3k_seurat, group.by = "seurat_annotations", split.by = "condition", replicate.by = "replicate")
```

    ## Error:
    ## ! object 'pbmc3k_seurat' not found

``` r

AbundancePlot(pbmc3k_seurat, group.by = "seurat_annotations", split.by = "condition", replicate.by = "replicate", paired = TRUE, draw_paths = TRUE, sina_shift = FALSE)
```

    ## Error in `AbundancePlot()`:
    ## ! unused argument (sina_shift = FALSE)

### `AssignQuantiles`

``` r

pbmc3k_seurat <- AssignQuantiles(pbmc3k_seurat, feature = "B2M", assay = "RNA", slot = "data", split.by = "condition", quantile.probs = c(0,0.4,0.8,1))
```

    ## Error in `AssignQuantiles()`:
    ## ! unused argument (slot = "data")

``` r

RidgePlot(pbmc3k_seurat, group.by = "B2M.quantile", features = "B2M")
```

    ## Error:
    ## ! object 'pbmc3k_seurat' not found

### `BackGatePlot`

``` r

BackgatePlot(pbmc3k_seurat, feature1 = "CD4", feature2 = "CD8A", metadata.col = "seurat_annotations", metadata.selection = "Memory CD4 T")
```

    ## Error:
    ## ! object 'pbmc3k_seurat' not found

``` r

BackgatePlot(pbmc3k_seurat, feature1 = "CD4", feature2 = "CD8A", metadata.col = "seurat_annotations", metadata.selection = "CD8 T")
```

    ## Error:
    ## ! object 'pbmc3k_seurat' not found

### `Heatmaps`

``` r

rna_markers <- presto::wilcoxauc(pbmc3k_seurat, group_by = "seurat_annotations")
```

    ## Error in `loadNamespace()`:
    ## ! there is no package called 'presto'

``` r

top_rna_markers <- rna_markers %>% 
    dplyr::filter(padj < 0.01, logFC > 0) %>%
    dplyr::group_by(group) %>% 
    dplyr::slice_min(padj, with_ties = FALSE, n = 3)
```

    ## Error:
    ## ! object 'rna_markers' not found

``` r

DoStarHeatmap(pbmc3k_seurat, diff_exp_results = rna_markers %>% dplyr::filter(group != "Platelet"), assay = "RNA", slot = "data", group.by = "seurat_annotations", p_val_choice = 0.01, logFC_choice = 2)
```

    ## Error in `DoStarHeatmap()`:
    ## ! unused argument (slot = "data")

``` r

DoClusteredHeatmap(pbmc3k_seurat, features = top_rna_markers$feature, assay = "RNA", group.by = "seurat_annotations")
```

    ## Error:
    ## ! object 'pbmc3k_seurat' not found

### `MetaDataPlot`

``` r

MetaDataPlot(pbmc3k_seurat, group.by = "seurat_annotations", split.by = "replicate")
```

    ## Error:
    ## ! object 'pbmc3k_seurat' not found

### `ModulePlot`

``` r

pbmc3k_seurat <- CellCycleScoring(pbmc3k_seurat, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
```

    ## Error:
    ## ! object 'pbmc3k_seurat' not found

``` r

ModulePlot(pbmc3k_seurat, features = "S.Score", assay = "RNA", split.by = "condition", replicate.by = "replicate")
```

    ## Error:
    ## ! object 'pbmc3k_seurat' not found

## Session Info

``` r

sessionInfo()
```

    ## R version 4.6.0 (2026-04-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] tidyr_1.3.2        ggplot2_4.0.3      pochi_0.1.0        Seurat_5.5.0      
    ## [5] SeuratObject_5.4.0 sp_2.2-1          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3         
    ##   [4] rlang_1.2.0            magrittr_2.0.5         RcppAnnoy_0.0.23      
    ##   [7] otel_0.2.0             spatstat.geom_3.7-3    matrixStats_1.5.0     
    ##  [10] ggridges_0.5.7         compiler_4.6.0         png_0.1-9             
    ##  [13] systemfonts_1.3.2      vctrs_0.7.3            reshape2_1.4.5        
    ##  [16] stringr_1.6.0          pkgconfig_2.0.3        fastmap_1.2.0         
    ##  [19] promises_1.5.0         rmarkdown_2.31         ragg_1.5.2            
    ##  [22] purrr_1.2.2            xfun_0.57              cachem_1.1.0          
    ##  [25] jsonlite_2.0.0         goftest_1.2-3          later_1.4.8           
    ##  [28] spatstat.utils_3.2-2   irlba_2.3.7            parallel_4.6.0        
    ##  [31] cluster_2.1.8.2        R6_2.6.1               ica_1.0-3             
    ##  [34] spatstat.data_3.1-9    bslib_0.10.0           stringi_1.8.7         
    ##  [37] RColorBrewer_1.1-3     reticulate_1.46.0      spatstat.univar_3.1-7 
    ##  [40] parallelly_1.47.0      lmtest_0.9-40          jquerylib_0.1.4       
    ##  [43] scattermore_1.2        Rcpp_1.1.1-1.1         knitr_1.51            
    ##  [46] tensor_1.5.1           future.apply_1.20.2    zoo_1.8-15            
    ##  [49] sctransform_0.4.3      httpuv_1.6.17          Matrix_1.7-5          
    ##  [52] splines_4.6.0          igraph_2.3.1           tidyselect_1.2.1      
    ##  [55] abind_1.4-8            yaml_2.3.12            spatstat.random_3.4-5 
    ##  [58] spatstat.explore_3.8-0 codetools_0.2-20       miniUI_0.1.2          
    ##  [61] listenv_0.10.1         lattice_0.22-9         tibble_3.3.1          
    ##  [64] plyr_1.8.9             withr_3.0.2            shiny_1.13.0          
    ##  [67] S7_0.2.2               ROCR_1.0-12            evaluate_1.0.5        
    ##  [70] Rtsne_0.17             future_1.70.0          fastDummies_1.7.6     
    ##  [73] desc_1.4.3             survival_3.8-6         polyclip_1.10-7       
    ##  [76] fitdistrplus_1.2-6     pillar_1.11.1          KernSmooth_2.23-26    
    ##  [79] plotly_4.12.0          generics_0.1.4         RcppHNSW_0.6.0        
    ##  [82] scales_1.4.0           globals_0.19.1         xtable_1.8-8          
    ##  [85] glue_1.8.1             lazyeval_0.2.3         tools_4.6.0           
    ##  [88] data.table_1.18.4      RSpectra_0.16-2        RANN_2.6.2            
    ##  [91] fs_2.1.0               dotCall64_1.2          cowplot_1.2.0         
    ##  [94] grid_4.6.0             nlme_3.1-169           patchwork_1.3.2       
    ##  [97] cli_3.6.6              spatstat.sparse_3.1-0  textshaping_1.0.5     
    ## [100] spam_2.11-3            viridisLite_0.4.3      dplyr_1.2.1           
    ## [103] uwot_0.2.4             gtable_0.3.6           sass_0.4.10           
    ## [106] digest_0.6.39          progressr_0.19.0       ggrepel_0.9.8         
    ## [109] htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.9       
    ## [112] pkgdown_2.2.0          lifecycle_1.0.5        httr_1.4.8            
    ## [115] mime_0.13              MASS_7.3-65

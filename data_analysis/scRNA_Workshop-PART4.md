---
title: "Introduction to Single Cell RNAseq Part 4"
author: "UCD Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

Last Updated: July 15, 2022

# Part 4: PCA and choice in number of PCS

## Load libraries

```r
library(Seurat)
library(biomaRt)
library(knitr)
library(ggplot2)
```

## Load the Seurat object

```r
load(file="pre_sample_corrected.RData")
experiment.aggregate
```

```
## An object of class Seurat 
## 21005 features across 10595 samples within 1 assay 
## Active assay: RNA (21005 features, 5986 variable features)
```

## Scale the data

ScaleData - Scales and centers genes in the dataset. If variables are provided in vars.to.regress, they are individually regressed against each gene, and the resulting residuals are then scaled and centered unless otherwise specified. Here we regress out cell cycle results S.Score and G2M.Score, percentage mitochondria (percent.mito) and the number of features (nFeature_RNA).


```r
experiment.aggregate <- ScaleData(
  object = experiment.aggregate,
  vars.to.regress = c("S.Score", "G2M.Score", "percent.mito", "nFeature_RNA"))
```

## Dimensionality reduction with PCA

Next we perform PCA (principal components analysis) on the scaled data.  


```r
?RunPCA
```


```r
experiment.aggregate <- RunPCA(object = experiment.aggregate, npcs=100)
```

Seurat then provides a number of ways to visualize the PCA results

Visualize PCA loadings

```r
VizDimLoadings(experiment.aggregate, dims = 1, ncol = 1) + theme_minimal(base_size = 8)
```

![](scRNA_Workshop-PART4_files/figure-html/viz_pca-1.png)<!-- -->

```r
VizDimLoadings(experiment.aggregate, dims = 2, ncol = 1) + theme_minimal(base_size = 8)
```

![](scRNA_Workshop-PART4_files/figure-html/viz_pca-2.png)<!-- -->

Principal components plot

```r
DimPlot(object = experiment.aggregate, reduction = "pca")
```

![](scRNA_Workshop-PART4_files/figure-html/plot_pca-1.png)<!-- -->

Draws a heatmap focusing on a principal component. Both cells and genes are sorted by their principal component scores. Allows for nice visualization of sources of heterogeneity in the dataset.


```r
DimHeatmap(object = experiment.aggregate, dims = 1:6, cells = 500, balanced = TRUE)
```

![](scRNA_Workshop-PART4_files/figure-html/heatmap_pca-1.png)<!-- -->

```r
DimHeatmap(object = experiment.aggregate, dims = 7:12, cells = 500, balanced = TRUE)
```

![](scRNA_Workshop-PART4_files/figure-html/heatmap_pca-2.png)<!-- -->

#### Questions

1. Go back to the original data (rerun the load RData section) and then try modifying the ScaleData vars.to.regres, remove some variables, try adding in orig.ident? See how choices effect the pca plot

### Selecting which PCs to use
To overcome the extensive technical noise in any single gene, Seurat clusters cells based on their PCA scores, with each PC essentially representing a metagene that combines information across a correlated gene set. Determining how many PCs to include downstream is therefore an important step.

ElbowPlot plots the standard deviations (or approximate singular values if running PCAFast) of the principle components for easy identification of an elbow in the graph. This elbow often corresponds well with the significant PCs and is much faster to run.  This is the traditional approach to selecting principal components.


```r
ElbowPlot(experiment.aggregate, ndims = 100)
```

![](scRNA_Workshop-PART4_files/figure-html/elbow-1.png)<!-- -->

The JackStraw function randomly permutes a subset of data, and calculates projected PCA scores for these 'random' genes, then compares the PCA scores for the 'random' genes with the observed PCA scores to determine statistical signifance. End result is a p-value for each gene's association with each principal component. We identify significant PCs as those who have a strong enrichment of low p-value genes.


```r
experiment.aggregate <- JackStraw(object = experiment.aggregate, dims = 100)
```


```r
experiment.aggregate <- ScoreJackStraw(experiment.aggregate, dims = 1:100)
JackStrawPlot(object = experiment.aggregate, dims = 1:100) + theme(legend.position="bottom")
```

![](scRNA_Workshop-PART4_files/figure-html/plot_jackstraw-1.png)<!-- -->

## Finally, lets save the filtered and normalized data

```r
save(experiment.aggregate, file="pca_sample_corrected.RData")
```

## Get the next Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2022-July-Single-Cell-RNA-Seq-Analysis/main/data_analysis/scRNA_Workshop-PART5.Rmd", "scRNA_Workshop-PART5.Rmd")
```

## Session Information

```r
sessionInfo()
```

```
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Catalina 10.15.7
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggplot2_3.3.5      knitr_1.37         biomaRt_2.50.2     SeuratObject_4.0.4
## [5] Seurat_4.1.0      
## 
## loaded via a namespace (and not attached):
##   [1] BiocFileCache_2.2.1    plyr_1.8.6             igraph_1.2.11         
##   [4] lazyeval_0.2.2         splines_4.1.2          listenv_0.8.0         
##   [7] scattermore_0.7        GenomeInfoDb_1.30.0    digest_0.6.29         
##  [10] htmltools_0.5.2        fansi_1.0.2            magrittr_2.0.2        
##  [13] memoise_2.0.1          tensor_1.5             cluster_2.1.2         
##  [16] ROCR_1.0-11            globals_0.14.0         Biostrings_2.62.0     
##  [19] matrixStats_0.61.0     spatstat.sparse_2.1-0  prettyunits_1.1.1     
##  [22] colorspace_2.0-2       rappdirs_0.3.3         blob_1.2.2            
##  [25] ggrepel_0.9.1          xfun_0.29              dplyr_1.0.8           
##  [28] crayon_1.5.0           RCurl_1.98-1.5         jsonlite_1.8.0        
##  [31] spatstat.data_2.1-2    survival_3.2-13        zoo_1.8-9             
##  [34] glue_1.6.2             polyclip_1.10-0        gtable_0.3.0          
##  [37] zlibbioc_1.40.0        XVector_0.34.0         leiden_0.3.9          
##  [40] future.apply_1.8.1     BiocGenerics_0.40.0    abind_1.4-5           
##  [43] scales_1.1.1           DBI_1.1.2              miniUI_0.1.1.1        
##  [46] Rcpp_1.0.8.3           progress_1.2.2         viridisLite_0.4.0     
##  [49] xtable_1.8-4           reticulate_1.24        spatstat.core_2.3-2   
##  [52] bit_4.0.4              stats4_4.1.2           htmlwidgets_1.5.4     
##  [55] httr_1.4.2             RColorBrewer_1.1-2     ellipsis_0.3.2        
##  [58] ica_1.0-2              farver_2.1.0           pkgconfig_2.0.3       
##  [61] XML_3.99-0.8           dbplyr_2.1.1           sass_0.4.0            
##  [64] uwot_0.1.11            deldir_1.0-6           utf8_1.2.2            
##  [67] labeling_0.4.2         tidyselect_1.1.2       rlang_1.0.2           
##  [70] reshape2_1.4.4         later_1.3.0            AnnotationDbi_1.56.2  
##  [73] munsell_0.5.0          tools_4.1.2            cachem_1.0.6          
##  [76] cli_3.2.0              generics_0.1.2         RSQLite_2.2.9         
##  [79] ggridges_0.5.3         evaluate_0.14          stringr_1.4.0         
##  [82] fastmap_1.1.0          yaml_2.3.5             goftest_1.2-3         
##  [85] bit64_4.0.5            fitdistrplus_1.1-6     purrr_0.3.4           
##  [88] RANN_2.6.1             KEGGREST_1.34.0        pbapply_1.5-0         
##  [91] future_1.23.0          nlme_3.1-155           mime_0.12             
##  [94] xml2_1.3.3             compiler_4.1.2         rstudioapi_0.13       
##  [97] filelock_1.0.2         curl_4.3.2             plotly_4.10.0         
## [100] png_0.1-7              spatstat.utils_2.3-0   tibble_3.1.6          
## [103] bslib_0.3.1            stringi_1.7.6          highr_0.9             
## [106] lattice_0.20-45        Matrix_1.4-0           vctrs_0.3.8           
## [109] pillar_1.7.0           lifecycle_1.0.1        spatstat.geom_2.3-1   
## [112] lmtest_0.9-39          jquerylib_0.1.4        RcppAnnoy_0.0.19      
## [115] data.table_1.14.2      cowplot_1.1.1          bitops_1.0-7          
## [118] irlba_2.3.5            httpuv_1.6.5           patchwork_1.1.1       
## [121] R6_2.5.1               promises_1.2.0.1       KernSmooth_2.23-20    
## [124] gridExtra_2.3          IRanges_2.28.0         parallelly_1.30.0     
## [127] codetools_0.2-18       MASS_7.3-55            assertthat_0.2.1      
## [130] withr_2.4.3            sctransform_0.3.3      S4Vectors_0.32.3      
## [133] GenomeInfoDbData_1.2.7 hms_1.1.1              mgcv_1.8-38           
## [136] parallel_4.1.2         grid_4.1.2             rpart_4.1.16          
## [139] tidyr_1.2.0            rmarkdown_2.11         Rtsne_0.15            
## [142] Biobase_2.54.0         shiny_1.7.1
```

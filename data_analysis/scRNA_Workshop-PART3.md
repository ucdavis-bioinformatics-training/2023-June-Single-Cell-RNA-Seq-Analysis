---
title: "Introduction to Single Cell RNAseq Part 3"
author: "UCD Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

Last Updated: December 7, 2022

# Part 3: PCA and choice in number of PCS

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

![](scRNA_Workshop-PART3_files/figure-html/viz_pca-1.png)<!-- -->

```r
VizDimLoadings(experiment.aggregate, dims = 2, ncol = 1) + theme_minimal(base_size = 8)
```

![](scRNA_Workshop-PART3_files/figure-html/viz_pca-2.png)<!-- -->

Principal components plot

```r
DimPlot(object = experiment.aggregate, reduction = "pca")
```

![](scRNA_Workshop-PART3_files/figure-html/plot_pca-1.png)<!-- -->

Draws a heatmap focusing on a principal component. Both cells and genes are sorted by their principal component scores. Allows for nice visualization of sources of heterogeneity in the dataset.


```r
DimHeatmap(object = experiment.aggregate, dims = 1:6, cells = 500, balanced = TRUE)
```

![](scRNA_Workshop-PART3_files/figure-html/heatmap_pca-1.png)<!-- -->

```r
DimHeatmap(object = experiment.aggregate, dims = 7:12, cells = 500, balanced = TRUE)
```

![](scRNA_Workshop-PART3_files/figure-html/heatmap_pca-2.png)<!-- -->

#### Questions

1. Go back to the original data (rerun the load RData section) and then try modifying the ScaleData vars.to.regres, remove some variables, try adding in orig.ident? See how choices effect the pca plot

### Selecting which PCs to use
To overcome the extensive technical noise in any single gene, Seurat clusters cells based on their PCA scores, with each PC essentially representing a metagene that combines information across a correlated gene set. Determining how many PCs to include downstream is therefore an important step.

ElbowPlot plots the standard deviations (or approximate singular values if running PCAFast) of the principle components for easy identification of an elbow in the graph. This elbow often corresponds well with the significant PCs and is much faster to run.  This is the traditional approach to selecting principal components.


```r
ElbowPlot(experiment.aggregate, ndims = 100)
```

![](scRNA_Workshop-PART3_files/figure-html/elbow-1.png)<!-- -->

The JackStraw function randomly permutes a subset of data, and calculates projected PCA scores for these 'random' genes, then compares the PCA scores for the 'random' genes with the observed PCA scores to determine statistical signifance. End result is a p-value for each gene's association with each principal component. We identify significant PCs as those who have a strong enrichment of low p-value genes.


```r
experiment.aggregate <- JackStraw(object = experiment.aggregate, dims = 100)
```


```r
experiment.aggregate <- ScoreJackStraw(experiment.aggregate, dims = 1:100)
JackStrawPlot(object = experiment.aggregate, dims = 1:100) + theme(legend.position="bottom")
```

![](scRNA_Workshop-PART3_files/figure-html/plot_jackstraw-1.png)<!-- -->

## Finally, lets save the filtered and normalized data

```r
save(experiment.aggregate, file="pca_sample_corrected.RData")
```

## Get the next Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2022-December-Single-Cell-RNA-Seq-Analysis/main/data_analysis/scRNA_Workshop-PART4.Rmd", "scRNA_Workshop-PART4.Rmd")
```

## Session Information

```r
sessionInfo()
```

```
## R version 4.2.2 (2022-10-31)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Catalina 10.15.7
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggplot2_3.4.0      knitr_1.41         biomaRt_2.54.0     SeuratObject_4.1.3
## [5] Seurat_4.3.0      
## 
## loaded via a namespace (and not attached):
##   [1] BiocFileCache_2.6.0    plyr_1.8.8             igraph_1.3.5          
##   [4] lazyeval_0.2.2         sp_1.5-1               splines_4.2.2         
##   [7] listenv_0.8.0          scattermore_0.8        GenomeInfoDb_1.34.3   
##  [10] digest_0.6.30          htmltools_0.5.3        fansi_1.0.3           
##  [13] magrittr_2.0.3         memoise_2.0.1          tensor_1.5            
##  [16] cluster_2.1.4          ROCR_1.0-11            globals_0.16.2        
##  [19] Biostrings_2.66.0      matrixStats_0.63.0     spatstat.sparse_3.0-0 
##  [22] prettyunits_1.1.1      colorspace_2.0-3       rappdirs_0.3.3        
##  [25] blob_1.2.3             ggrepel_0.9.2          xfun_0.35             
##  [28] dplyr_1.0.10           crayon_1.5.2           RCurl_1.98-1.9        
##  [31] jsonlite_1.8.3         progressr_0.11.0       spatstat.data_3.0-0   
##  [34] survival_3.4-0         zoo_1.8-11             glue_1.6.2            
##  [37] polyclip_1.10-4        gtable_0.3.1           zlibbioc_1.44.0       
##  [40] XVector_0.38.0         leiden_0.4.3           future.apply_1.10.0   
##  [43] BiocGenerics_0.44.0    abind_1.4-5            scales_1.2.1          
##  [46] DBI_1.1.3              spatstat.random_3.0-1  miniUI_0.1.1.1        
##  [49] Rcpp_1.0.9             progress_1.2.2         viridisLite_0.4.1     
##  [52] xtable_1.8-4           reticulate_1.26        bit_4.0.5             
##  [55] stats4_4.2.2           htmlwidgets_1.5.4      httr_1.4.4            
##  [58] RColorBrewer_1.1-3     ellipsis_0.3.2         ica_1.0-3             
##  [61] farver_2.1.1           pkgconfig_2.0.3        XML_3.99-0.12         
##  [64] dbplyr_2.2.1           sass_0.4.4             uwot_0.1.14           
##  [67] deldir_1.0-6           utf8_1.2.2             labeling_0.4.2        
##  [70] tidyselect_1.2.0       rlang_1.0.6            reshape2_1.4.4        
##  [73] later_1.3.0            AnnotationDbi_1.60.0   munsell_0.5.0         
##  [76] tools_4.2.2            cachem_1.0.6           cli_3.4.1             
##  [79] generics_0.1.3         RSQLite_2.2.19         ggridges_0.5.4        
##  [82] evaluate_0.18          stringr_1.4.1          fastmap_1.1.0         
##  [85] yaml_2.3.6             goftest_1.2-3          bit64_4.0.5           
##  [88] fitdistrplus_1.1-8     purrr_0.3.5            RANN_2.6.1            
##  [91] KEGGREST_1.38.0        pbapply_1.6-0          future_1.29.0         
##  [94] nlme_3.1-160           mime_0.12              xml2_1.3.3            
##  [97] compiler_4.2.2         rstudioapi_0.14        filelock_1.0.2        
## [100] curl_4.3.3             plotly_4.10.1          png_0.1-8             
## [103] spatstat.utils_3.0-1   tibble_3.1.8           bslib_0.4.1           
## [106] stringi_1.7.8          highr_0.9              lattice_0.20-45       
## [109] Matrix_1.5-3           vctrs_0.5.1            pillar_1.8.1          
## [112] lifecycle_1.0.3        spatstat.geom_3.0-3    lmtest_0.9-40         
## [115] jquerylib_0.1.4        RcppAnnoy_0.0.20       bitops_1.0-7          
## [118] data.table_1.14.6      cowplot_1.1.1          irlba_2.3.5.1         
## [121] httpuv_1.6.6           patchwork_1.1.2        R6_2.5.1              
## [124] promises_1.2.0.1       KernSmooth_2.23-20     gridExtra_2.3         
## [127] IRanges_2.32.0         parallelly_1.32.1      codetools_0.2-18      
## [130] MASS_7.3-58.1          assertthat_0.2.1       withr_2.5.0           
## [133] sctransform_0.3.5      GenomeInfoDbData_1.2.9 S4Vectors_0.36.0      
## [136] hms_1.1.2              parallel_4.2.2         grid_4.2.2            
## [139] tidyr_1.2.1            rmarkdown_2.18         Rtsne_0.16            
## [142] spatstat.explore_3.0-5 Biobase_2.58.0         shiny_1.7.3
```

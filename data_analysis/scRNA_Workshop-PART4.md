---
title: "Introduction to Single Cell RNAseq Part 4"
author: "UCD Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

Last Updated: March 20 2020

# Part 4: PCA and choice in number of PCS

## Load libraries

```r
library(Seurat)
library(ggplot2)
```

## Load the Seurat object

```r
load(file="pre_sample_corrected.RData")
experiment.aggregate
```

```
## An object of class Seurat 
## 36601 features across 3343 samples within 1 assay 
## Active assay: RNA (36601 features, 3941 variable features)
```


## Now doing so for 'real'

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
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2022-March-Single-Cell-RNA-Seq-Analysis/master/data_analysis/scRNA_Workshop-PART5.Rmd", "scRNA_Workshop-PART5.Rmd")
```

## Session Information

```r
sessionInfo()
```

```
## R version 4.1.2 (2021-11-01)
## Platform: aarch64-apple-darwin20 (64-bit)
## Running under: macOS Monterey 12.0.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggplot2_3.3.5      SeuratObject_4.0.4 Seurat_4.1.0      
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.15            colorspace_2.0-3      deldir_1.0-6         
##   [4] ellipsis_0.3.2        ggridges_0.5.3        rstudioapi_0.13      
##   [7] spatstat.data_2.1-2   farver_2.1.0          leiden_0.3.9         
##  [10] listenv_0.8.0         ggrepel_0.9.1         fansi_1.0.2          
##  [13] codetools_0.2-18      splines_4.1.2         knitr_1.37           
##  [16] polyclip_1.10-0       jsonlite_1.8.0        ica_1.0-2            
##  [19] cluster_2.1.2         png_0.1-7             uwot_0.1.11          
##  [22] shiny_1.7.1           sctransform_0.3.3     spatstat.sparse_2.1-0
##  [25] compiler_4.1.2        httr_1.4.2            assertthat_0.2.1     
##  [28] Matrix_1.4-0          fastmap_1.1.0         lazyeval_0.2.2       
##  [31] cli_3.2.0             later_1.3.0           htmltools_0.5.2      
##  [34] tools_4.1.2           igraph_1.2.11         gtable_0.3.0         
##  [37] glue_1.6.2            RANN_2.6.1            reshape2_1.4.4       
##  [40] dplyr_1.0.8           Rcpp_1.0.8.3          scattermore_0.8      
##  [43] jquerylib_0.1.4       vctrs_0.3.8           nlme_3.1-155         
##  [46] lmtest_0.9-39         spatstat.random_2.1-0 xfun_0.30            
##  [49] stringr_1.4.0         globals_0.14.0        mime_0.12            
##  [52] miniUI_0.1.1.1        lifecycle_1.0.1       irlba_2.3.5          
##  [55] goftest_1.2-3         future_1.24.0         MASS_7.3-55          
##  [58] zoo_1.8-9             scales_1.1.1          spatstat.core_2.4-0  
##  [61] promises_1.2.0.1      spatstat.utils_2.3-0  parallel_4.1.2       
##  [64] RColorBrewer_1.1-2    yaml_2.3.5            reticulate_1.24      
##  [67] pbapply_1.5-0         gridExtra_2.3         sass_0.4.0           
##  [70] rpart_4.1.16          stringi_1.7.6         highr_0.9            
##  [73] rlang_1.0.2           pkgconfig_2.0.3       matrixStats_0.61.0   
##  [76] evaluate_0.15         lattice_0.20-45       ROCR_1.0-11          
##  [79] purrr_0.3.4           tensor_1.5            labeling_0.4.2       
##  [82] patchwork_1.1.1       htmlwidgets_1.5.4     cowplot_1.1.1        
##  [85] tidyselect_1.1.2      parallelly_1.30.0     RcppAnnoy_0.0.19     
##  [88] plyr_1.8.6            magrittr_2.0.2        R6_2.5.1             
##  [91] generics_0.1.2        DBI_1.1.2             withr_2.5.0          
##  [94] mgcv_1.8-39           pillar_1.7.0          fitdistrplus_1.1-8   
##  [97] survival_3.3-1        abind_1.4-5           tibble_3.1.6         
## [100] future.apply_1.8.1    crayon_1.5.0          KernSmooth_2.23-20   
## [103] utf8_1.2.2            spatstat.geom_2.3-2   plotly_4.10.0        
## [106] rmarkdown_2.13        grid_4.1.2            data.table_1.14.2    
## [109] digest_0.6.29         xtable_1.8-4          tidyr_1.2.0          
## [112] httpuv_1.6.5          munsell_0.5.0         viridisLite_0.4.0    
## [115] bslib_0.3.1
```

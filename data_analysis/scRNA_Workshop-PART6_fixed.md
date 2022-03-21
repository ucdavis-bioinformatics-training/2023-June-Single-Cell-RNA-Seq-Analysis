---
title: "Single Cell RNAseq Part 6"
author: "Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

# Part 6: Enrichment, Model-Based DE, and Cell-Type Identification



## Load libraries

```r
library(Seurat)
library(ggplot2)
library(limma)
library(topGO)
```

## Load the Seurat object

```r
load("clusters_seurat_object.RData")
experiment.merged
```

<div class='r_output'> An object of class Seurat 
 36601 features across 3343 samples within 1 assay 
 Active assay: RNA (36601 features, 3941 variable features)
  3 dimensional reductions calculated: pca, tsne, umap
</div>
```r
Idents(experiment.merged) <- "finalcluster"
```

# 1. Gene Ontology (GO) Enrichment of Genes Expressed in a Cluster
[Gene Ontology](http://geneontology.org/docs/ontology-documentation/) provides a controlled vocabulary for describing gene products.  Here we use enrichment analysis to identify GO terms that are overrepresented among the gene expressed in cells in a given cluster. 


```r
cluster0 <- subset(experiment.merged, idents = '0')
expr <- as.matrix(GetAssayData(cluster0))
# Filter out genes that are 0 for every cell in this cluster
bad <- which(rowSums(expr) == 0)
expr <- expr[-bad,]

# Select genes that are expressed > 0 in at least half of cells
n.gt.0 <- apply(expr, 1, function(x)length(which(x > 0)))
expressed.genes <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.5)]
all.genes <- rownames(expr)

# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
names(geneList) <- all.genes

# Create topGOdata object
	GOdata <- new("topGOdata",
		ontology = "BP", # use biological process ontology
		allGenes = geneList,
		geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")
# Test for enrichment using Fisher's Exact Test
	resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
	GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
```

<div class='r_output'>         GO.ID                                                            Term Annotated Significant Expected  Fisher
 1  GO:0002181                                         cytoplasmic translation       142         108    17.73 < 1e-30
 2  GO:0006364                                                 rRNA processing       215          97    26.84 6.2e-19
 3  GO:0006120            mitochondrial electron transport, NADH to ubiquinone        41          28     5.12 1.3e-16
 4  GO:0032981              mitochondrial respiratory chain complex I assembly        51          30     6.37 4.7e-15
 5  GO:0000027                                ribosomal large subunit assembly        24          17     3.00 5.8e-11
 6  GO:0000398                                  mRNA splicing, via spliceosome       265         103    33.08 7.8e-11
 7  GO:1904874 positive regulation of telomerase RNA localization to Cajal ...        15          13     1.87 1.4e-10
 8  GO:0042776            mitochondrial ATP synthesis coupled proton transport        16          13     2.00 6.6e-10
 9  GO:0001732         formation of cytoplasmic translation initiation complex        16          13     2.00 6.6e-10
 10 GO:0006457                                                 protein folding       173          68    21.60 9.7e-10
 11 GO:0000028                                ribosomal small subunit assembly        17          13     2.12 2.5e-09
 12 GO:0006412                                                     translation       570         226    71.16 2.7e-09
 13 GO:0048025           negative regulation of mRNA splicing, via spliceosome        18          13     2.25 7.9e-09
 14 GO:0006446                          regulation of translational initiation        72          30     8.99 2.9e-08
 15 GO:0006123        mitochondrial electron transport, cytochrome c to oxygen        12          10     1.50 4.6e-08
 16 GO:1904851 positive regulation of establishment of protein localization...        10           9     1.25 6.4e-08
 17 GO:0045727                              positive regulation of translation       115          36    14.36 8.1e-08
 18 GO:1904871       positive regulation of protein localization to Cajal body        11           9     1.37 3.1e-07
 19 GO:0010498                           proteasomal protein catabolic process       429          91    53.56 4.9e-07
 20 GO:0042273                              ribosomal large subunit biogenesis        69          45     8.61 5.2e-07
</div>* Annotated: number of genes (out of all.genes) that are annotated with that GO term
* Significant: number of genes that are annotated with that GO term and meet our criteria for "expressed"
* Expected: Under random chance, number of genes that would be expected to be annotated with that GO term and meeting our criteria for "expressed"
* Fisher: (Raw) p-value from Fisher's Exact Test

# 2. Model-based DE analysis in limma
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html) is an R package for differential expression analysis of bulk RNASeq and microarray data.  We apply it here to single cell data.

Limma can be used to fit any linear model to expression data and is useful for analyses that go beyond two-group comparisons.  A detailed tutorial of model specification in limma is available [here](https://ucdavis-bioinformatics-training.github.io/2021-June-RNA-Seq-Analysis/data_analysis/DE_Analysis_mm_with_quizzes) and in the [limma User's Guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf).


```r
# filter genes to those expressed in at least 10% of cells
keep <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.1)]
expr2 <- expr[keep,]

# Set up "design matrix" with statistical model
mm <- model.matrix(~0 + orig.ident, data = cluster0[[]])
head(mm)
```

<div class='r_output'>                             orig.identconv_COVID orig.identconv_MMR orig.identconv_Tdap orig.identnorm_COVID
 AAACGGGTCTGGGCCA-conv_COVID                    1                  0                   0                    0
 AAACGGGTCTTAGAGC-conv_COVID                    1                  0                   0                    0
 AAAGCAACATCTACGA-conv_COVID                    1                  0                   0                    0
 AACACGTCAATGGAAT-conv_COVID                    1                  0                   0                    0
 AACCATGAGCCACGTC-conv_COVID                    1                  0                   0                    0
 AACGTTGAGGCCATAG-conv_COVID                    1                  0                   0                    0
</div>
```r
tail(mm)
```

<div class='r_output'>                             orig.identconv_COVID orig.identconv_MMR orig.identconv_Tdap orig.identnorm_COVID
 TTGCGTCGTGAGGCTA-conv_Tdap                     0                  0                   1                    0
 TTTACTGCAGACGCAA-conv_Tdap                     0                  0                   1                    0
 TTTGCGCCAGTGGGAT-conv_Tdap                     0                  0                   1                    0
 TTTGTCAGTCCAAGTT-conv_Tdap                     0                  0                   1                    0
 AACCATGCATTCCTCG-norm_COVID                    0                  0                   0                    1
 CACCTTGAGATCCCGC-norm_COVID                    0                  0                   0                    1
</div>
```r
# Fit model in limma
fit <- lmFit(expr2, mm)
head(coef(fit))
```

<div class='r_output'>           orig.identconv_COVID orig.identconv_MMR orig.identconv_Tdap orig.identnorm_COVID
 LINC01128            0.1599865         0.13309920          0.14106590            0.0000000
 NOC2L                0.8371142         0.82052021          0.85604706            0.0000000
 ISG15                0.2404222         0.04410456          0.05140408            0.0000000
 TNFRSF18             0.2219353         0.23514976          0.25065392            0.0000000
 TNFRSF4              1.1613225         1.16327044          1.14487791            0.9122914
 SDF4                 1.0265788         0.98423684          0.99927764            0.0000000
</div>
```r
# Test convMMR - conv_COVID
contr <- makeContrasts(orig.identconv_MMR - orig.identconv_COVID, levels = colnames(coef(fit)))
contr
```

<div class='r_output'>                       Contrasts
 Levels                 orig.identconv_MMR - orig.identconv_COVID
   orig.identconv_COVID                                        -1
   orig.identconv_MMR                                           1
   orig.identconv_Tdap                                          0
   orig.identnorm_COVID                                         0
</div>
```r
fit2 <- contrasts.fit(fit, contrasts = contr)
fit2 <- eBayes(fit2)
out <- topTable(fit2, n = Inf, sort.by = "P")
head(out, 30)
```

<div class='r_output'>             logFC    AveExpr          t      P.Value    adj.P.Val           B
 MX1    -0.5898165 0.21663819 -14.550979 1.650358e-42 8.926786e-39 83.47076328
 IFITM1 -0.9798207 0.89658334 -14.164616 1.327366e-40 3.589860e-37 79.29776589
 IFI6   -0.4061648 0.15156917 -12.684540 1.372467e-33 2.474557e-30 63.92754305
 IRF7   -0.4552179 0.23304852 -11.302995 1.699649e-27 2.298350e-24 50.57210277
 ISG20  -0.3508907 0.23182357  -8.683973 2.305719e-17 2.494327e-14 28.38349260
 PLSCR1 -0.2805733 0.17871718  -8.313932 4.238105e-16 3.820652e-13 25.62070049
 SP100  -0.3518974 0.40242985  -7.616077 7.752661e-14 5.990592e-11 20.68417353
 IFI35  -0.2659176 0.20875367  -7.107738 2.714792e-12 1.835539e-09 17.32123061
 TRIM22 -0.2685543 0.20393164  -7.055481 3.867453e-12 2.324339e-09 16.98693562
 CD38   -0.2434399 0.17399476  -7.028966 4.624242e-12 2.501252e-09 16.81814032
 ISG15  -0.1963177 0.10704512  -6.942417 8.254235e-12 4.058833e-09 16.27103395
 SP110  -0.3270012 0.47765636  -6.647878 5.665578e-11 2.553759e-08 14.45387486
 EPSTI1 -0.1709894 0.11077193  -6.208515 8.781664e-10 3.653848e-07 11.87361276
 NAPA   -0.3170928 0.77339244  -5.763465 1.196254e-08 4.621813e-06  9.42237923
 IRF9   -0.2623210 0.43120260  -5.624120 2.617931e-08 9.440260e-06  8.68907289
 TYMP   -0.2297078 0.24187115  -5.610358 2.825924e-08 9.553391e-06  8.61754072
 BST2   -0.2305673 0.28082039  -5.476842 5.883230e-08 1.871905e-05  7.93192857
 CHMP5  -0.2785508 0.64583170  -5.369012 1.051833e-07 3.160759e-05  7.38932953
 PPM1K  -0.2299323 0.36600636  -5.172345 2.957540e-07 8.419650e-05  6.42543599
 GBP1   -0.2144269 0.31978732  -4.509679 7.514667e-06 2.032342e-03  3.42565479
 PARP9  -0.1040153 0.09555163  -4.463625 9.275191e-06 2.389024e-03  3.23155998
 LY6E   -0.2520333 0.76701973  -4.245856 2.446380e-05 6.014759e-03  2.33940161
 CD47    0.2744001 1.16229929   4.077628 5.027916e-05 1.182435e-02  1.67930352
 IFNGR2  0.1355285 0.20026671   3.935771 9.051168e-05 2.039907e-02  1.14251377
 AKT2   -0.1611336 0.31175504  -3.752394 1.884503e-04 4.077311e-02  0.47564290
 NRDC   -0.1394348 0.24128247  -3.738598 1.988979e-04 4.137842e-02  0.42670986
 MT-CO1  0.2330131 2.42863370   3.694673 2.359113e-04 4.726089e-02  0.27206630
 YWHAQ   0.2260659 2.09639042   3.655522 2.742715e-04 5.298338e-02  0.13571358
 GZMB   -0.2487144 0.31212605  -3.587107 3.557047e-04 6.634506e-02 -0.09918947
 HNRNPM  0.1955155 1.18172787   3.483723 5.227108e-04 9.424476e-02 -0.44602711
</div>
### Output columns:
* logFC: log2 fold change
* AveExpr: Average expression across all cells in expr2
* t: logFC divided by its standard error
* P.Value: Raw p-value (based on t) from test that logFC differs from 0
* adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
* B: log-odds that gene is DE 



# BONUS: Cell type identification with scMRMA
[scMRMA]([https://academic.oup.com/nar/article/50/2/e7/6396893]) (single cell Multi-Resolution Marker-based Annotation Algorithm) classifies cells by iteratively clustering them then annotating based on a hierarchical external database.

The databases included with the current version are only for use with human and mouse, but a user-constructed hierarchichal database can be used. 

The package can be installed from [Github](https://github.com/JiaLiVUMC/scMRMA):


```r
# Remove hashes to run
# install.packages("devtools")
# devtools::install_github("JiaLiVUMC/scMRMA")
```


```r
suppressPackageStartupMessages(library(scMRMA))
result <- scMRMA(input = experiment.merged,
                 species = "Hs",
                 db = "TcellAI")
```

<div class='r_output'> Pre-defined cell type database TcellAI will be used.
 Multi Resolution Annotation Started. 
 Level 1 annotation started. 
 Level 2 annotation started. 
 Uniform Resolution Annotation Started.
</div>
```r
table(result$uniformR$annotationResult)
```

<div class='r_output'> 
     CD4_naive           Th1 Gamma_delta_T         nTreg 
          2014          1136           147            46
</div>
```r
## Add cell types to metadata
experiment.merged <- AddMetaData(experiment.merged, result$uniformR$annotationResult, col.name = "CellType")
table(experiment.merged$CellType, experiment.merged$orig.ident)
```

<div class='r_output'>                
                 conv_COVID conv_MMR conv_Tdap norm_COVID
   CD4_naive            413      625       465        511
   Th1                  367      398       371          0
   Gamma_delta_T          0        0         0        147
   nTreg                 13       10        23          0
</div>
```r
table(experiment.merged$CellType, experiment.merged$finalcluster)
```

<div class='r_output'>                
                   0   1   2   3   4   5   7   8  10  12  13
   CD4_naive     701 414 338 198 295  23  13   0   1   2  29
   Th1            53 295  23   0   1 156 177 186 174  71   0
   Gamma_delta_T   0   0   0 135   2   0   0   0   0   0  10
   nTreg           6   5   1   0   0  34   0   0   0   0   0
</div>
```r
DimPlot(experiment.merged, group.by = "CellType", label = TRUE)
```

![](scRNA_Workshop-PART6_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

                 
## Session Information

```r
sessionInfo()
```

<div class='r_output'> R version 4.1.3 (2022-03-10)
 Platform: x86_64-w64-mingw32/x64 (64-bit)
 Running under: Windows 10 x64 (build 19044)
 
 Matrix products: default
 
 locale:
 [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                           LC_TIME=English_United States.1252    
 
 attached base packages:
 [1] stats4    stats     graphics  grDevices datasets  utils     methods   base     
 
 other attached packages:
  [1] scMRMA_1.0           networkD3_0.4        data.tree_1.0.0      tidyr_1.2.0          RANN_2.6.1           plyr_1.8.6           irlba_2.3.5          Matrix_1.4-0         org.Hs.eg.db_3.14.0  topGO_2.46.0         SparseM_1.81         GO.db_3.14.0         AnnotationDbi_1.56.2 IRanges_2.28.0       S4Vectors_0.32.3     Biobase_2.54.0       graph_1.72.0         BiocGenerics_0.40.0  limma_3.50.1         ggplot2_3.3.5        SeuratObject_4.0.4  
 [22] Seurat_4.1.0        
 
 loaded via a namespace (and not attached):
   [1] igraph_1.2.11          lazyeval_0.2.2         splines_4.1.3          listenv_0.8.0          scattermore_0.8        usethis_2.1.5          GenomeInfoDb_1.30.1    digest_0.6.29          htmltools_0.5.2        fansi_1.0.2            magrittr_2.0.2         memoise_2.0.1          tensor_1.5             cluster_2.1.2          ROCR_1.0-11            remotes_2.4.2          globals_0.14.0         Biostrings_2.62.0      matrixStats_0.61.0    
  [20] spatstat.sparse_2.1-0  prettyunits_1.1.1      colorspace_2.0-3       blob_1.2.2             ggrepel_0.9.1          xfun_0.30              dplyr_1.0.8            callr_3.7.0            crayon_1.5.0           RCurl_1.98-1.6         jsonlite_1.8.0         spatstat.data_2.1-2    survival_3.2-13        zoo_1.8-9              glue_1.6.2             polyclip_1.10-0        gtable_0.3.0           zlibbioc_1.40.0        XVector_0.34.0        
  [39] leiden_0.3.9           pkgbuild_1.3.1         future.apply_1.8.1     abind_1.4-5            scales_1.1.1           DBI_1.1.2              spatstat.random_2.1-0  miniUI_0.1.1.1         Rcpp_1.0.8             viridisLite_0.4.0      xtable_1.8-4           reticulate_1.24        spatstat.core_2.4-0    bit_4.0.4              htmlwidgets_1.5.4      httr_1.4.2             RColorBrewer_1.1-2     ellipsis_0.3.2         ica_1.0-2             
  [58] farver_2.1.0           pkgconfig_2.0.3        sass_0.4.0             uwot_0.1.11            deldir_1.0-6           utf8_1.2.2             labeling_0.4.2         tidyselect_1.1.2       rlang_1.0.2            reshape2_1.4.4         later_1.3.0            munsell_0.5.0          tools_4.1.3            cachem_1.0.6           cli_3.2.0              generics_0.1.2         RSQLite_2.2.10         devtools_2.4.3         ggridges_0.5.3        
  [77] evaluate_0.15          stringr_1.4.0          fastmap_1.1.0          yaml_2.3.5             goftest_1.2-3          processx_3.5.2         fs_1.5.2               knitr_1.37             bit64_4.0.5            fitdistrplus_1.1-8     purrr_0.3.4            KEGGREST_1.34.0        pbapply_1.5-0          future_1.24.0          nlme_3.1-155           mime_0.12              brio_1.1.3             compiler_4.1.3         rstudioapi_0.13       
  [96] plotly_4.10.0          png_0.1-7              testthat_3.1.2         spatstat.utils_2.3-0   tibble_3.1.6           bslib_0.3.1            stringi_1.7.6          highr_0.9              ps_1.6.0               desc_1.4.1             lattice_0.20-45        vctrs_0.3.8            pillar_1.7.0           lifecycle_1.0.1        spatstat.geom_2.3-2    lmtest_0.9-39          jquerylib_0.1.4        RcppAnnoy_0.0.19       data.table_1.14.2     
 [115] cowplot_1.1.1          bitops_1.0-7           httpuv_1.6.5           patchwork_1.1.1        R6_2.5.1               promises_1.2.0.1       renv_0.15.4            KernSmooth_2.23-20     gridExtra_2.3          parallelly_1.30.0      sessioninfo_1.2.2      codetools_0.2-18       pkgload_1.2.4          MASS_7.3-55            rprojroot_2.0.2        withr_2.5.0            sctransform_0.3.3      GenomeInfoDbData_1.2.7 mgcv_1.8-39           
 [134] parallel_4.1.3         grid_4.1.3             rpart_4.1.16           rmarkdown_2.13         Rtsne_0.15             shiny_1.7.1
</div>
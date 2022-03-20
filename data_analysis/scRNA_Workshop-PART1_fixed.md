---
title: "Introduction to Single Cell RNAseq Part 1"
author: "UCD Bioinformatics Core"
output:
    html_document:
      keep_md: TRUE
---

Last Updated: March 20 2022

# Part 1: Loading data from CellRanger into R

Our first Markdown document concentrates on getting data into R and setting up our initial object.

## Single Cell Analysis with Seurat and some custom code!

[Seurat](http://satijalab.org/seurat/) (now Version 4) is a popular R package that is designed for QC, analysis, and exploration of single cell data. Seurat aims to enable users to identify and interpret sources of heterogeneity from single cell transcriptomic measurements, and to integrate diverse types of single cell data. Further, the authors provide several [tutorials](https://satijalab.org/seurat/vignettes.html), on their website.

The **expression_data_cellranger.zip** file contains the single cell matrix files and HDF5 files for four T cell samples from [Mysore et al., 2021](https://www.cell.com/med/fulltext/S2666-6340(21)00289-0).

We start each markdown document with loading needed libraries for R:



```r
# must have Seurat
library(Seurat)
library(kableExtra)
library(ggplot2)
```


### Setup the experiment folder and data info

```r
experiment_name = "Covid Example"
dataset_loc <- "./expression_data_cellranger"
ids <- c("conv_COVID", "conv_MMR", "conv_Tdap", "norm_COVID")
```


## Read in the cellranger sample metrics csv files

Usually, there would be a single metrics summary file for each sample. In this case, we are working with a pool of samples that has been demultiplexed using antibody capture data prior to the workshop, and only have one metrics summary file generated for the entire pool.

The code in the box below reads in a single metrics summary file. **To read in a summary file for each sample in your experiment, you should use the code in the next box instead.**


```r
experiment.metrics <- read.csv(file.path(dataset_loc, "metrics_summary.csv"))
sequencing.metrics <- data.frame(t(experiment.metrics[,c(1:19)]))
rownames(sequencing.metrics) <- gsub("\\.", " ", rownames(sequencing.metrics))
colnames(sequencing.metrics) <- "All samples"
sequencing.metrics %>%
  kable(caption = 'Cell Ranger Results') %>%
  pack_rows("Overview", 1, 3, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Sequencing Characteristics", 4, 9, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Mapping Characteristics", 10, 19, label_row_css = "background-color: #666; color: #fff;") %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>Cell Ranger Results</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> All samples </th>
  </tr>
 </thead>
<tbody>
  <tr grouplength="3"><td colspan="2" style="background-color: #666; color: #fff;"><strong>Overview</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Estimated Number of Cells </td>
   <td style="text-align:left;"> 7,072 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Mean Reads per Cell </td>
   <td style="text-align:left;"> 89,005 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Median Genes per Cell </td>
   <td style="text-align:left;"> 1,880 </td>
  </tr>
  <tr grouplength="6"><td colspan="2" style="background-color: #666; color: #fff;"><strong>Sequencing Characteristics</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Number of Reads </td>
   <td style="text-align:left;"> 629,444,626 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Valid Barcodes </td>
   <td style="text-align:left;"> 94.4% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Sequencing Saturation </td>
   <td style="text-align:left;"> 84.0% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Q30 Bases in Barcode </td>
   <td style="text-align:left;"> 94.3% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Q30 Bases in RNA Read </td>
   <td style="text-align:left;"> 89.5% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Q30 Bases in UMI </td>
   <td style="text-align:left;"> 93.8% </td>
  </tr>
  <tr grouplength="10"><td colspan="2" style="background-color: #666; color: #fff;"><strong>Mapping Characteristics</strong></td></tr>
<tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped to Genome </td>
   <td style="text-align:left;"> 90.0% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped Confidently to Genome </td>
   <td style="text-align:left;"> 85.4% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped Confidently to Intergenic Regions </td>
   <td style="text-align:left;"> 1.7% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped Confidently to Intronic Regions </td>
   <td style="text-align:left;"> 3.5% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped Confidently to Exonic Regions </td>
   <td style="text-align:left;"> 80.2% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped Confidently to Transcriptome </td>
   <td style="text-align:left;"> 76.5% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Reads Mapped Antisense to Gene </td>
   <td style="text-align:left;"> 1.7% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Fraction Reads in Cells </td>
   <td style="text-align:left;"> 75.5% </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Total Genes Detected </td>
   <td style="text-align:left;"> 19,923 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Median UMI Counts per Cell </td>
   <td style="text-align:left;"> 5,671 </td>
  </tr>
</tbody>
</table>


```r
d10x.metrics <- lapply(ids, function(i){
  metrics <- read.csv(file.path(dataset_loc,paste0(i,"/outs"),"metrics_summary.csv"), colClasses = "character")
})
experiment.metrics <- do.call("rbind", d10x.metrics)
rownames(experiment.metrics) <- ids

sequencing.metrics <- data.frame(t(experiment.metrics[,c(1:19)]))

row.names(sequencing.metrics) <- gsub("\\."," ", rownames(sequencing.metrics))

sequencing.metrics %>%
  kable(caption = 'Cell Ranger Results') %>%
  pack_rows("Overview", 1, 3, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Sequencing Characteristics", 4, 9, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Mapping Characteristics", 10, 19, label_row_css = "background-color: #666; color: #fff;") %>%
  kable_styling("striped")
```


## Load the Cell Ranger Matrix Data and create the base Seurat object.
Cell Ranger provides a function `cellranger aggr` that will combine multiple samples into a single matrix file. However, when processing data in R this is unnecessary and we can quickly aggregate them in R.

Seurat provides a function `Read10X` and `Read10X_h5` to read in 10X data folder. First we read in data from each individual sample folder. 

Later, we initialize the Seurat object (`CreateSeuratObject`) with the raw (non-normalized data). Keep all cells with at least 200 detected genes. Also extracting sample names, calculating and adding in the metadata mitochondrial percentage of each cell. Adding in the metadata batchid and cell cycle. Finally, saving the raw Seurat object.

## Load the Cell Ranger Matrix Data (hdf5 file) and create the base Seurat object.

```r
d10x.data <- lapply(ids, function(i){
  d10x <- Read10X_h5(file.path(dataset_loc, i, "/outs","raw_feature_bc_matrix.h5"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
names(d10x.data) <- ids

str(d10x.data)
```

<div class='r_output'> List of 4
  $ conv_COVID:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
   .. ..@ i       : int [1:3341791] 24 42 43 44 53 59 61 62 78 83 ...
   .. ..@ p       : int [1:22472] 0 0 0 0 3471 3471 3471 3471 3478 3478 ...
   .. ..@ Dim     : int [1:2] 36601 22471
   .. ..@ Dimnames:List of 2
   .. .. ..$ : chr [1:36601] "MIR1302-2HG" "FAM138A" "OR4F5" "AL627309.1" ...
   .. .. ..$ : chr [1:22471] "AAACCTGAGAACAACT-conv_COVID" "AAACCTGAGAAGGCCT-conv_COVID" "AAACCTGAGACACGAC-conv_COVID" "AAACCTGAGACTAGAT-conv_COVID" ...
   .. ..@ x       : num [1:3341791] 3 2 11 6 2 1 1 1 1 2 ...
   .. ..@ factors : list()
  $ conv_MMR  :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
   .. ..@ i       : int [1:4372343] 9776 13612 4881 29788 17208 18735 2096 13949 20008 24103 ...
   .. ..@ p       : int [1:39830] 0 1 1 2 2 2 4 4 4 4 ...
   .. ..@ Dim     : int [1:2] 36601 39829
   .. ..@ Dimnames:List of 2
   .. .. ..$ : chr [1:36601] "MIR1302-2HG" "FAM138A" "OR4F5" "AL627309.1" ...
   .. .. ..$ : chr [1:39829] "AAACCTGAGAATCTCC-conv_MMR" "AAACCTGAGACCGGAT-conv_MMR" "AAACCTGAGACCTTTG-conv_MMR" "AAACCTGAGACGCACA-conv_MMR" ...
   .. ..@ x       : num [1:4372343] 1 1 1 1 1 1 1 1 1 1 ...
   .. ..@ factors : list()
  $ conv_Tdap :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
   .. ..@ i       : int [1:3801450] 4325 28917 28335 32724 5414 31575 184 187 190 219 ...
   .. ..@ p       : int [1:23839] 0 2 2 4 6 821 822 822 823 4513 ...
   .. ..@ Dim     : int [1:2] 36601 23838
   .. ..@ Dimnames:List of 2
   .. .. ..$ : chr [1:36601] "MIR1302-2HG" "FAM138A" "OR4F5" "AL627309.1" ...
   .. .. ..$ : chr [1:23838] "AAACCTGAGAACTCGG-conv_Tdap" "AAACCTGAGACTGGGT-conv_Tdap" "AAACCTGAGAGGACGG-conv_Tdap" "AAACCTGAGAGTCTGG-conv_Tdap" ...
   .. ..@ x       : num [1:3801450] 1 1 1 1 1 1 1 1 1 2 ...
   .. ..@ factors : list()
  $ norm_COVID:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
   .. ..@ i       : int [1:966879] 13909 25538 19314 31392 32724 11613 33848 2429 11932 24 ...
   .. ..@ p       : int [1:17821] 0 0 0 0 2 2 5 5 5 5 ...
   .. ..@ Dim     : int [1:2] 36601 17820
   .. ..@ Dimnames:List of 2
   .. .. ..$ : chr [1:36601] "MIR1302-2HG" "FAM138A" "OR4F5" "AL627309.1" ...
   .. .. ..$ : chr [1:17820] "AAACCTGAGAGGGATA-norm_COVID" "AAACCTGAGATAGCAT-norm_COVID" "AAACCTGAGATCCTGT-norm_COVID" "AAACCTGAGCCACCTG-norm_COVID" ...
   .. ..@ x       : num [1:966879] 1 1 1 1 1 1 1 1 1 1 ...
   .. ..@ factors : list()
</div>
If you don't have the needed hdf5 libraries you can read in the matrix files like such


```r
d10x.data <- sapply(ids, function(i){
  d10x <- Read10X(file.path(dataset_loc, i, "/outs","raw_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x), split="-"), '[[', 1L), i, sep="-")
  d10x
})
names(d10x.data) <- ids
```


Lets recreate the barcode rank plot from the Cell Ranger web summary file.

In this case, we only have estimated cell counts for the pool as a whole. To color each barcode rank plot using the number of cells estimated by Cell Ranger, skip the code in the box below and use the next one instead.

```r
plot_cellranger_cells <- function(ind){
  xbreaks = c(1,1e1,1e2,1e3,1e4,1e5,1e6)
  xlabels = c("1","10","100","1000","10k","100K","1M")
  ybreaks = c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000)
  ylabels = c("1","2","5","10","2","5","100","2","5","1000","2","5","10k","2","5","100K","2","5","1M")

  pl1 <- data.frame(index=seq.int(1,ncol(d10x.data[[ind]])),
                    nCount_RNA = sort(Matrix:::colSums(d10x.data[[ind]])+1,decreasing=T),
                    nFeature_RNA = sort(Matrix:::colSums(d10x.data[[ind]]>0)+1,decreasing=T)) %>%
    ggplot() + 
    scale_color_manual(values=c("red2","blue4"), labels=c("Features", "UMI"), name=NULL) +
    ggtitle(paste("CellRanger filltered cells:",ids[ind],sep=" ")) + xlab("Barcodes") + ylab("counts (UMI or Features") + 
    scale_x_continuous(trans = 'log2', breaks=xbreaks, labels = xlabels) + 
    scale_y_continuous(trans = 'log2', breaks=ybreaks, labels = ylabels) +
    geom_line(aes(x=index, y=nCount_RNA, color = "UMI"), size=1.75) +
    geom_line(aes(x=index, y=nFeature_RNA, color = "Features"), size=1.25)

  return(pl1)
}

plot_cellranger_cells(1)
```

![](scRNA_Workshop-PART1_files/figure-html/fig_barcode_umi-1.png)<!-- -->

```r
plot_cellranger_cells(2)
```

![](scRNA_Workshop-PART1_files/figure-html/fig_barcode_umi-2.png)<!-- -->

```r
plot_cellranger_cells(3)
```

![](scRNA_Workshop-PART1_files/figure-html/fig_barcode_umi-3.png)<!-- -->

```r
plot_cellranger_cells(4)
```

![](scRNA_Workshop-PART1_files/figure-html/fig_barcode_umi-4.png)<!-- -->


```r
cr_filtered_cells <- as.numeric(gsub(",","",as.character(unlist(sequencing.metrics["Estimated Number of Cells",]))))

plot_cellranger_cells <- function(ind){
  xbreaks = c(1,1e1,1e2,1e3,1e4,1e5,1e6)
  xlabels = c("1","10","100","1000","10k","100K","1M")
  ybreaks = c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000)
  ylabels = c("1","2","5","10","2","5","100","2","5","1000","2","5","10k","2","5","100K","2","5","1M")

  pl1 <- data.frame(index=seq.int(1,ncol(d10x.data[[ind]])), nCount_RNA = sort(Matrix:::colSums(d10x.data[[ind]])+1,decreasing=T), nFeature_RNA = sort(Matrix:::colSums(d10x.data[[ind]]>0)+1,decreasing=T)) %>% ggplot() + 
    scale_color_manual(values=c("grey50","red2","blue4"), labels=c("UMI_Background", "Features", "UMI_Cells"), name=NULL) +
    ggtitle(paste("CellRanger filltered cells:",ids[ind],sep=" ")) + xlab("Barcodes") + ylab("counts (UMI or Features") + 
    scale_x_continuous(trans = 'log2', breaks=xbreaks, labels = xlabels) + 
    scale_y_continuous(trans = 'log2', breaks=ybreaks, labels = ylabels) +
    geom_line(aes(x=index, y=nCount_RNA, color=index<=cr_filtered_cells[ind] , group=1), size=1.75) +
    geom_line(aes(x=index, y=nFeature_RNA, color="Features", group=1), size=1.25)

  return(pl1)
}

plot_cellranger_cells(1)
plot_cellranger_cells(2)
plot_cellranger_cells(3)
plot_cellranger_cells(4)
```


### Create the Seurat object

Filter criteria: remove genes that do not occur in a minimum of 0 cells and remove cells that don't have a minimum of 200 features


```r
experiment.data <- do.call("cbind", d10x.data)

experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = experiment_name,
  min.cells = 0,
  min.features = 300,
  names.field = 2,
  names.delim = "\\-")

experiment.aggregate
```

<div class='r_output'> An object of class Seurat 
 36601 features across 5803 samples within 1 assay 
 Active assay: RNA (36601 features, 0 variable features)
</div>
```r
str(experiment.aggregate)
```

<div class='r_output'> Formal class 'Seurat' [package "SeuratObject"] with 13 slots
   ..@ assays      :List of 1
   .. ..$ RNA:Formal class 'Assay' [package "SeuratObject"] with 8 slots
   .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
   .. .. .. .. .. ..@ i       : int [1:12282940] 24 42 43 44 53 59 61 62 78 83 ...
   .. .. .. .. .. ..@ p       : int [1:5804] 0 3471 7135 9695 13117 13746 14071 17010 20128 22256 ...
   .. .. .. .. .. ..@ Dim     : int [1:2] 36601 5803
   .. .. .. .. .. ..@ Dimnames:List of 2
   .. .. .. .. .. .. ..$ : chr [1:36601] "MIR1302-2HG" "FAM138A" "OR4F5" "AL627309.1" ...
   .. .. .. .. .. .. ..$ : chr [1:5803] "AAACCTGAGACTAGAT-conv_COVID" "AAACCTGAGCTACCTA-conv_COVID" "AAACCTGCAGACTCGC-conv_COVID" "AAACCTGGTAAATGTG-conv_COVID" ...
   .. .. .. .. .. ..@ x       : num [1:12282940] 3 2 11 6 2 1 1 1 1 2 ...
   .. .. .. .. .. ..@ factors : list()
   .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
   .. .. .. .. .. ..@ i       : int [1:12282940] 24 42 43 44 53 59 61 62 78 83 ...
   .. .. .. .. .. ..@ p       : int [1:5804] 0 3471 7135 9695 13117 13746 14071 17010 20128 22256 ...
   .. .. .. .. .. ..@ Dim     : int [1:2] 36601 5803
   .. .. .. .. .. ..@ Dimnames:List of 2
   .. .. .. .. .. .. ..$ : chr [1:36601] "MIR1302-2HG" "FAM138A" "OR4F5" "AL627309.1" ...
   .. .. .. .. .. .. ..$ : chr [1:5803] "AAACCTGAGACTAGAT-conv_COVID" "AAACCTGAGCTACCTA-conv_COVID" "AAACCTGCAGACTCGC-conv_COVID" "AAACCTGGTAAATGTG-conv_COVID" ...
   .. .. .. .. .. ..@ x       : num [1:12282940] 3 2 11 6 2 1 1 1 1 2 ...
   .. .. .. .. .. ..@ factors : list()
   .. .. .. ..@ scale.data   : num[0 , 0 ] 
   .. .. .. ..@ key          : chr "rna_"
   .. .. .. ..@ assay.orig   : NULL
   .. .. .. ..@ var.features : logi(0) 
   .. .. .. ..@ meta.features:'data.frame':	36601 obs. of  0 variables
   .. .. .. ..@ misc         : list()
   ..@ meta.data   :'data.frame':	5803 obs. of  3 variables:
   .. ..$ orig.ident  : Factor w/ 4 levels "conv_COVID","conv_MMR",..: 1 1 1 1 1 1 1 1 1 1 ...
   .. ..$ nCount_RNA  : num [1:5803] 17521 16250 7321 16235 813 ...
   .. ..$ nFeature_RNA: int [1:5803] 3471 3664 2560 3422 629 325 2939 3118 2128 2964 ...
   ..@ active.assay: chr "RNA"
   ..@ active.ident: Factor w/ 4 levels "conv_COVID","conv_MMR",..: 1 1 1 1 1 1 1 1 1 1 ...
   .. ..- attr(*, "names")= chr [1:5803] "AAACCTGAGACTAGAT-conv_COVID" "AAACCTGAGCTACCTA-conv_COVID" "AAACCTGCAGACTCGC-conv_COVID" "AAACCTGGTAAATGTG-conv_COVID" ...
   ..@ graphs      : list()
   ..@ neighbors   : list()
   ..@ reductions  : list()
   ..@ images      : list()
   ..@ project.name: chr "Covid Example"
   ..@ misc        : list()
   ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
   .. ..$ : int [1:3] 4 0 4
   ..@ commands    : list()
   ..@ tools       : list()
</div>
### The percentage of reads that map to the mitochondrial genome

* Low-quality / dying cells often exhibit extensive mitochondrial contamination.
* We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features.
* We use the set of all genes, in mouse these genes can be identified as those that begin with 'mt', in human data they begin with MT.


```r
experiment.aggregate$percent.mito <- PercentageFeatureSet(experiment.aggregate, pattern = "^MT-")
summary(experiment.aggregate$percent.mito)
```

<div class='r_output'>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.0000  0.5711  1.1171  2.7418  3.8196 48.2639
</div>
### Lets spend a little time getting to know the Seurat object.

The Seurat object is the center of each single cell analysis. It stores __all__ information associated with the dataset, including data, annotations, analyses, etc. The R function slotNames can be used to view the slot names within an object.


```r
slotNames(experiment.aggregate)
```

<div class='r_output'>  [1] "assays"       "meta.data"    "active.assay" "active.ident" "graphs"      
  [6] "neighbors"    "reductions"   "images"       "project.name" "misc"        
 [11] "version"      "commands"     "tools"
</div>

```r
head(experiment.aggregate[[]])
```

<div class='r_output'>                             orig.ident nCount_RNA nFeature_RNA percent.mito
 AAACCTGAGACTAGAT-conv_COVID conv_COVID      17521         3471    0.5878660
 AAACCTGAGCTACCTA-conv_COVID conv_COVID      16250         3664    0.4307692
 AAACCTGCAGACTCGC-conv_COVID conv_COVID       7321         2560    0.3824614
 AAACCTGGTAAATGTG-conv_COVID conv_COVID      16235         3422    0.8192177
 AAACCTGTCTTCCTTC-conv_COVID conv_COVID        813          629    1.7220172
 AAACGGGAGAGTCGGT-conv_COVID conv_COVID        452          325    2.4336283
</div>
#### Question(s)

1. What slots are empty, what slots have data?
2. What columns are available in meta.data?
3. Look up the help documentation for subset?


## Finally, save the original object and view the object.

Original data set in Seurat class, with no filtering

```r
save(experiment.aggregate,file="original_seurat_object.RData")
```

## Get the next Rmd file

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2022-March-Single-Cell-RNA-Seq-Analysis/main/data_analysis/scRNA_Workshop-PART2.Rmd", "scRNA_Workshop-PART2.Rmd")
```

## Session Information

```r
sessionInfo()
```

<div class='r_output'> R version 4.1.2 (2021-11-01)
 Platform: aarch64-apple-darwin20 (64-bit)
 Running under: macOS Monterey 12.0.1
 
 Matrix products: default
 BLAS:   /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/lib/libRblas.0.dylib
 LAPACK: /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/lib/libRlapack.dylib
 
 locale:
 [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
 
 attached base packages:
 [1] stats     graphics  grDevices utils     datasets  methods   base     
 
 other attached packages:
 [1] ggplot2_3.3.5      kableExtra_1.3.4   SeuratObject_4.0.4 Seurat_4.1.0      
 
 loaded via a namespace (and not attached):
   [1] Rtsne_0.15            colorspace_2.0-3      deldir_1.0-6         
   [4] ellipsis_0.3.2        ggridges_0.5.3        rstudioapi_0.13      
   [7] spatstat.data_2.1-2   farver_2.1.0          leiden_0.3.9         
  [10] listenv_0.8.0         bit64_4.0.5           ggrepel_0.9.1        
  [13] fansi_1.0.2           xml2_1.3.3            codetools_0.2-18     
  [16] splines_4.1.2         knitr_1.37            polyclip_1.10-0      
  [19] jsonlite_1.8.0        ica_1.0-2             cluster_2.1.2        
  [22] png_0.1-7             uwot_0.1.11           shiny_1.7.1          
  [25] sctransform_0.3.3     spatstat.sparse_2.1-0 compiler_4.1.2       
  [28] httr_1.4.2            assertthat_0.2.1      Matrix_1.4-0         
  [31] fastmap_1.1.0         lazyeval_0.2.2        cli_3.2.0            
  [34] later_1.3.0           htmltools_0.5.2       tools_4.1.2          
  [37] igraph_1.2.11         gtable_0.3.0          glue_1.6.2           
  [40] RANN_2.6.1            reshape2_1.4.4        dplyr_1.0.8          
  [43] Rcpp_1.0.8.3          scattermore_0.8       jquerylib_0.1.4      
  [46] vctrs_0.3.8           svglite_2.1.0         nlme_3.1-155         
  [49] lmtest_0.9-39         spatstat.random_2.1-0 xfun_0.30            
  [52] stringr_1.4.0         globals_0.14.0        rvest_1.0.2          
  [55] mime_0.12             miniUI_0.1.1.1        lifecycle_1.0.1      
  [58] irlba_2.3.5           goftest_1.2-3         future_1.24.0        
  [61] MASS_7.3-55           zoo_1.8-9             scales_1.1.1         
  [64] spatstat.core_2.4-0   promises_1.2.0.1      spatstat.utils_2.3-0 
  [67] parallel_4.1.2        RColorBrewer_1.1-2    yaml_2.3.5           
  [70] reticulate_1.24       pbapply_1.5-0         gridExtra_2.3        
  [73] sass_0.4.0            rpart_4.1.16          stringi_1.7.6        
  [76] highr_0.9             systemfonts_1.0.4     rlang_1.0.2          
  [79] pkgconfig_2.0.3       matrixStats_0.61.0    evaluate_0.15        
  [82] lattice_0.20-45       ROCR_1.0-11           purrr_0.3.4          
  [85] tensor_1.5            patchwork_1.1.1       htmlwidgets_1.5.4    
  [88] bit_4.0.4             cowplot_1.1.1         tidyselect_1.1.2     
  [91] parallelly_1.30.0     RcppAnnoy_0.0.19      plyr_1.8.6           
  [94] magrittr_2.0.2        R6_2.5.1              generics_0.1.2       
  [97] DBI_1.1.2             withr_2.5.0           mgcv_1.8-39          
 [100] pillar_1.7.0          fitdistrplus_1.1-8    survival_3.3-1       
 [103] abind_1.4-5           tibble_3.1.6          future.apply_1.8.1   
 [106] hdf5r_1.3.5           crayon_1.5.0          KernSmooth_2.23-20   
 [109] utf8_1.2.2            spatstat.geom_2.3-2   plotly_4.10.0        
 [112] rmarkdown_2.13        grid_4.1.2            data.table_1.14.2    
 [115] webshot_0.5.2         digest_0.6.29         xtable_1.8-4         
 [118] tidyr_1.2.0           httpuv_1.6.5          munsell_0.5.0        
 [121] viridisLite_0.4.0     bslib_0.3.1
</div>
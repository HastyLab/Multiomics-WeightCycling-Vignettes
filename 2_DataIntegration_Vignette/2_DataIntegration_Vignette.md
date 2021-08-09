---
title: "2_DataIntegration_Vignette"
author: "Matthew A. Cottam"
date: "8/8/2021"
output:
  html_document: 
    toc: yes
    number_sections: yes
    keep_md: true
  html_notebook: 
    theme: united
    toc_float: yes
    toc: yes
    fig_height: 6
    number_sections: yes
---
# In this vignette
In this data integration vignette, we will be taking the output from the Preprocessing vignette and merging our unique samples together into one unified Seurat object. We will also perform dimensional reduction on the integrated Seurat object. Herein, we will do the following:

1. Normalize the RNA and ADT assays
2. Identify integration features
3. Run PCA on each sample independently
4. Integrate data using rPCA
5. Perform non-linear dimensional reduction
6. Visualize unclustered dimensional reduction

There are multiple methods for data integration, some of which have been incorporated into Seurat V4. In this experiment, all four samples were prepared simultaneously and ultimately run on the 10X Chromium platform on the same chip. Therefore, we expect batch effects to be very small. For computational efficiency, we have integrated these data sets using the `rPCA` (reciprocal PCR) approach.

# Prepare the R environment
## Load in the necessary libraries

```r
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(SeuratWrappers)
library(plotly)
library(RColorBrewer)
```

## Load in our Seurat object list
This `.rds` file comes from Vignette #1 - Data Preprocessing.

```r
Seurat_list <- readRDS("../1_Preprocessing_Vignette/1_Preprocessing.rds")
```

# Prepare Seurat objects for integration
## Normalize RNA assay
First, we will normalize and find variable features for the RNA assay containing our gene expression data. The variable features will be used to determine integration features.


```r
Seurat_list <- lapply(X = Seurat_list, FUN = function(x) {
    x <- NormalizeData(x, assay="RNA", normalization.method = "LogNormalize", verbose = FALSE)
    x <- FindVariableFeatures(x, assay="RNA", selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})
```

## Normalize ADT assay
Like for the `HTO` assay in the previous vignette, we will use CLR normalization for the `ADT` assay.

```r
Seurat_list <- lapply(X = Seurat_list, FUN = function(x) {
    x <- NormalizeData(x, assay="ADT", normalization.method="CLR", verbose= FALSE)
    x <- ScaleData(x, assay= "ADT", verbose= FALSE)
})
```

## Select integration features
The rPCA approach requires scaling of the data and subsequent PCA of each Seurat object prior to integration. Here, we are scaling on the integration features to produce our PCA results.

```r
features <- SelectIntegrationFeatures(object.list = Seurat_list)
Seurat_list <- lapply(X = Seurat_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
```

# Perform integration using rPCA
Data integration retains metaData for each individual Seurat object. The “project name”, or diet in this experiment, is stored in the `orig.ident` metaData column.

## Determine integration anchors

```r
options(future.globals.maxSize= 5242880000) #Integration requires a larger global size limit than default. I've set it to 5 Gb, but ~2 Gb was required.
anchors <- FindIntegrationAnchors(object.list = Seurat_list, anchor.features = features, reduction = "rpca", verbose=FALSE)
```

```
## Warning in CheckDuplicateCellNames(object.list = object.list): Some cell names
## are duplicated across objects provided. Renaming to enforce unique cell names.
```

## Integrate objects and generated an integrated assay
A new assay, called `integrated` will be generated upon data integration. We will use this assay for clustering, but will always revert to the `RNA` assay for differential expression.

```r
data.integrated <- IntegrateData(anchorset = anchors)
```

```
## Merging dataset 2 into 3
```

```
## Extracting anchors for merged samples
```

```
## Finding integration vectors
```

```
## Finding integration vector weights
```

```
## Integrating data
```

```
## Merging dataset 4 into 3 2
```

```
## Extracting anchors for merged samples
```

```
## Finding integration vectors
```

```
## Finding integration vector weights
```

```
## Integrating data
```

```
## Merging dataset 1 into 3 2 4
```

```
## Extracting anchors for merged samples
```

```
## Finding integration vectors
```

```
## Finding integration vector weights
```

```
## Integrating data
```

# Dimensional Reduction
As mentioned before, the `integrated` data assay will be used for clustering and visualization.

## Scale the data and perform PCA on the integrated assay
PCA is a linear dimensional reduction approach that seeks to identify the most variable sets of features. I've chosen to run 50 PCs for this experiment.

```r
# Change default assay
DefaultAssay(data.integrated) <- "integrated"

# Scale the data and perform PCA on the integrated assay
data.integrated <- ScaleData(data.integrated, assay="integrated", verbose=FALSE)
data.integrated <- RunPCA(data.integrated, assay="integrated", npcs = 50, verbose=FALSE)
```

## Identify how many dimensions to use
Multiple methods can be used to determine how many PCs to use for downstream analysis. Here, we use both a statistical approach `JackStraw` and a heuristic approach `ElbowPlot` as recommended in Seurat’s PBMC3k tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html.

### JackStraw
First, we will calculate and score via JackStraw. Warning: Jackstraw takes a long time to run.

```r
# Calculate the JackStraw visualization
data.integrated <- JackStraw(data.integrated, assay= "integrated", dims= 50, num.replicate= 100, verbose= FALSE)

# Score JackStraw 
data.integrated <- ScoreJackStraw(data.integrated, dims = 1:50)

# Plot JackStraw
JackStrawPlot(data.integrated, dims=1:50) + theme(legend.position="bottom")
```

```
## Warning: Removed 70000 rows containing missing values (geom_point).
```

![](2_DataIntegration_Vignette_files/figure-html/jackstraw-1.png)<!-- -->

### ElbowPlot
The `ElbowPlot` simply ranks PCs by percentage of variance.

```r
ElbowPlot(data.integrated, ndims=50)
```

![](2_DataIntegration_Vignette_files/figure-html/elbow-1.png)<!-- -->

## Nonlinear Dimensional Reduction
Uniform Manifold Projection (UMAP) is the preferred visualization for large datasets. It runs more efficiently than tSNE and preserves global structure much better. However, both work well for visualization and we run both to store them as embeddings in case we need to look at both. We use all 50 PCAs, since all appear to be significant. 

Note: We calculate here both 2-dimensional and 3-dimensional embeddings. 3-dimensional embeddings can be useful to capture interactions between cell types that are otherwise very different when only looking at two dimensions.

### Calculate UMAP

```r
data.integrated <- RunUMAP(data.integrated, dims=1:50, verbose=FALSE)
```

```
## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
## This message will be shown once per session
```

```r
data.integrated <- RunUMAP(data.integrated, dims=1:50, n.components=3, reduction.name="umap3D", reduction.key="UMAP3D_", verbose=FALSE)
```

### Calculate tSNE

```r
data.integrated <- RunTSNE(data.integrated, dims=1:50, verbose=FALSE)
data.integrated <- RunTSNE(data.integrated, dims=1:50, dim.embed=3, reduction.name="tsne3D", reduction.key="tSNE3D_")
```

# Data Projection
## Replace embeddings
Warning: UMAP in particular seems extremely sensitive to changes in OS and to any updates in the underlying software. These changes seem to result in rotations or flips across axes, even when setting the same seed. This was the case when rerunning this vignette and specifically seems to be the number of edges defined when determining the epochs (1503980 here vs. 1501094). I've fixed this by simply replacing the embeddings with those calculated previously. To do so, I've saved the original embeddings as an `rds` file and have placed them in the folder containing the second vignette. You could easily skip this and continue on, but the exact numbers of the subclusters (but not the associated genes) in the next vignette may be different and therefore annotations could be wrong. 

If you'd like to try your hand at this from the start or would like to reannotate, feel free to comment out this chunk. All steps will still be the same. For transparency sake, I've plotted the 2D visualization for UMAP and tSNE calculated before replacing the embeddings for comparison.

```r
d1 <- DimPlot(data.integrated, reduction="umap", label=F) + NoLegend() + labs(subtitle="UMAP")
d2 <- DimPlot(data.integrated, reduction="tsne", label=F) + NoLegend() + labs(subtitle="tSNE")
d1+d2
```

![](2_DataIntegration_Vignette_files/figure-html/2dcomp-1.png)<!-- -->


```r
#Replace reductions with previously calculated embeddings
embeddings <- readRDS("embeddings.rds")
data.integrated@reductions <- embeddings
```

## 2D visualization

```r
d1 <- DimPlot(data.integrated, reduction="umap", label=F) + NoLegend() + labs(subtitle="UMAP")
d2 <- DimPlot(data.integrated, reduction="tsne", label=F) + NoLegend() + labs(subtitle="tSNE")
d1+d2
```

![](2_DataIntegration_Vignette_files/figure-html/2dumap-1.png)<!-- -->
## 3D visualization
For 3-dimensional visualization, we use the `plotly` packaged to produce an interactive plot.
### UMAP

```r
plot.data <- FetchData(object = data.integrated, vars = c("UMAP3D_1", "UMAP3D_2", "UMAP3D_3"))
plot_ly(data = plot.data, #3D interactive plot of the data
        x = ~UMAP3D_1, y = ~UMAP3D_2, z = ~UMAP3D_3, 
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 0.5, width=0.5),
        hoverinfo="text")
```

```{=html}
<div id="htmlwidget-8708243e357302e0acf6" style="width:672px;height:480px;" class="plotly html-widget"></div>
```

### tSNE

```r
plot.data2 <- FetchData(object = data.integrated, vars = c("tSNE3D_1", "tSNE3D_2", "tSNE3D_3"))
plot_ly(data = plot.data2, #3D interactive plot of the data
        x = ~tSNE3D_1, y = ~tSNE3D_2, z = ~tSNE3D_3,
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 0.5, width=0.5),
        hoverinfo="text")
```

```{=html}
<div id="htmlwidget-54276fa55e6295f3f412" style="width:672px;height:480px;" class="plotly html-widget"></div>
```
# Save Progress

```r
saveRDS(data.integrated, file="2_DataIntegration.rds")
```
# Session Information

```r
sessionInfo()
```

```
## R version 4.1.0 (2021-05-18)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.2 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] RColorBrewer_1.1-2   plotly_4.9.3         SeuratWrappers_0.3.0
## [4] cowplot_1.1.1        dplyr_1.0.6          ggplot2_3.3.3       
## [7] SeuratObject_4.0.1   Seurat_4.0.2.9004   
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.15            colorspace_2.0-1      deldir_0.2-10        
##   [4] ellipsis_0.3.2        ggridges_0.5.3        spatstat.data_2.1-0  
##   [7] farver_2.1.0          leiden_0.3.8          listenv_0.8.0        
##  [10] remotes_2.4.0         ggrepel_0.9.1         RSpectra_0.16-0      
##  [13] fansi_0.5.0           codetools_0.2-18      splines_4.1.0        
##  [16] knitr_1.33            polyclip_1.10-0       jsonlite_1.7.2       
##  [19] ica_1.0-2             cluster_2.1.2         png_0.1-7            
##  [22] uwot_0.1.10           shiny_1.6.0           sctransform_0.3.2    
##  [25] spatstat.sparse_2.0-0 BiocManager_1.30.15   compiler_4.1.0       
##  [28] httr_1.4.2            assertthat_0.2.1      Matrix_1.3-4         
##  [31] fastmap_1.1.0         lazyeval_0.2.2        later_1.2.0          
##  [34] htmltools_0.5.1.1     tools_4.1.0           rsvd_1.0.5           
##  [37] igraph_1.2.6          gtable_0.3.0          glue_1.4.2           
##  [40] RANN_2.6.1            reshape2_1.4.4        Rcpp_1.0.6           
##  [43] scattermore_0.7       jquerylib_0.1.4       vctrs_0.3.8          
##  [46] nlme_3.1-152          crosstalk_1.1.1       lmtest_0.9-37        
##  [49] xfun_0.23             stringr_1.4.0         globals_0.14.0       
##  [52] mime_0.10             miniUI_0.1.1.1        lifecycle_1.0.0      
##  [55] irlba_2.3.3           goftest_1.2-2         future_1.21.0        
##  [58] MASS_7.3-54           zoo_1.8-9             scales_1.1.1         
##  [61] spatstat.core_2.1-2   promises_1.2.0.1      spatstat.utils_2.1-0 
##  [64] parallel_4.1.0        yaml_2.2.1            reticulate_1.20      
##  [67] pbapply_1.4-3         gridExtra_2.3         sass_0.4.0           
##  [70] rpart_4.1-15          stringi_1.6.2         highr_0.9            
##  [73] rlang_0.4.11          pkgconfig_2.0.3       matrixStats_0.59.0   
##  [76] evaluate_0.14         lattice_0.20-44       ROCR_1.0-11          
##  [79] purrr_0.3.4           tensor_1.5            labeling_0.4.2       
##  [82] patchwork_1.1.1       htmlwidgets_1.5.3     tidyselect_1.1.1     
##  [85] parallelly_1.25.0     RcppAnnoy_0.0.18      plyr_1.8.6           
##  [88] magrittr_2.0.1        R6_2.5.0              generics_0.1.0       
##  [91] DBI_1.1.1             pillar_1.6.1          withr_2.4.2          
##  [94] mgcv_1.8-36           fitdistrplus_1.1-5    survival_3.2-11      
##  [97] abind_1.4-5           tibble_3.1.2          future.apply_1.7.0   
## [100] crayon_1.4.1          KernSmooth_2.23-20    utf8_1.2.1           
## [103] spatstat.geom_2.1-0   rmarkdown_2.8         grid_4.1.0           
## [106] data.table_1.14.0     digest_0.6.27         xtable_1.8-4         
## [109] tidyr_1.1.3           httpuv_1.6.1          munsell_0.5.0        
## [112] viridisLite_0.4.0     bslib_0.2.5.1
```

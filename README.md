# Multiomics in weight cycling - Generating the integrated data file
This repository provides code and the minimally necessary files to generate the integrated Seurat v4 data object used in the following publication: **Under Review**. If you use code accessed through this github resource, please remember to cite our work:

**CITATION INFORMATION**

## What is contained in this repository?
This repository has three vignettes. We start with sequencing data that has already been processed through 10X Genomic's CellRanger V3 software. We provide the clustering information, the filtered barcode matrix, and the raw barcode matrix for each sample. The associated datafiles also include both gene expression and feature barcoding information necessary to reproduce our integrated data object.

1) Data Preprocessing - Quality control, ambient mRNA removal, and hashtag demultiplexing.
2) Data Integrated - Using reciprocal PCA to integrate our samples into one object and perform dimensional reduction.
3) Cluster Annotation - Applying multiple strategies to accurately annotate our cell clusters.

## Instructions 
To access our code for the three vignettes in the form of R Notebooks AND to download the minimally necessary files, clone this github repository or download the files as a .zip using this <a href="https://github.com/HastyLab/Multiomics-WeightCycling-Vignettes/archive/refs/heads/main.zip">download link</a>.

Alternatively, clone this repository directly using the command:
```
git clone https://github.com/HastyLab/Multiomics-WeightCycling-Vignettes.git
```

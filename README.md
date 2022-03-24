# Analysis of 10X Genomics scRNA-seq data
This Github describes the pipeline used by the Cebola Lab to analyse single-cell RNA-seq data (scRNA-seq) from 10X Genomics. 

## Table of contents

1. [Pipeline overview](#1-Pipeline-overview)  
2. [Background](#2-Background)
3. [Workspace setup](#3-Workspace-setup)
4. [Generate a reference transcriptome](#2-Generate-a-reference-transcriptome)
5. [Cellranger count](#5-Cellranger-count)
6. [Secondary analysis](#6-Secondary-analysis)
7. [Cell type identity](#7-cell-type-identity)

## 1. Pipeline overview

1. **Generate count matrix**: `CellRanger count` is used to generate count matrices with some initial filtering to remove empty droplets. 
2. **Pre-processing**: the count matrix is read into R using `Seurat`.  Initial pre-processing is carried out to prepare the data for the next steps, including normalization with a negative binominal model (`SCTransform`), merging of technical replicates and initial dimensionality reduction and clustering using `RunPCA`, `RunUMAP`, `FindNeighbors` and `FindClusters`.
3. **Correct for ambient gene expression** `SoupX` is used to correct for ambient gene expression. 
4. **QC filtering**: repeat the normalization and clustering with the corrected data. Identify and remove clusters of low-quality cells. Several rounds of pre-processing, clustering and filtering may be required.
5. **DoubletFinder**: identify and remove droplets with doublets i.e. two (or more) cells using `doubletFinder`. 
6. **Final clustering**: the cleaned data is processed for a final time and clusters are labelled using known marker genes. (Supervised and/or unsupervised clustering may be carried out). 

## 2. Background

This pipeline has been developed by carefully reviewing current tools and best practises used in the analysis of 10X Genomics scRNA-seq data, as of March 2022. This Github will first present an overview of various available tools, followed by the Cebola Lab pipeline. Resources used are shown in the [References](#references) at the bottom of this page.

### Alignment, demultiplexing and quantification

There are several alignment algorithms to choose from, including CellRanger, STARsolo, Alevin, Alevin-fry and Kallisto; these are compared in [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741). This pipeline will use the 10X Genomics toolbox, **CellRanger**. If memory requirement is an issue, [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) suggest STARsolo as an alternative.

### Secondary analysis

CellRanger includes the demultiplexing of reads from individual droplets, or 'GEMs' (Gel Beads in EMulsion), based on unique barcodes. These should *in theory* correspond to unique single cells. However, droplets can contain more or less than one cell, or contain damaged or low quality cells. Therefore, the secondary analysis will assess the quality of the dataset and filter the data to retain only droplets with high-quality single cells. 

The pipeline for secondary analysis, including references, is discussed below.

> **Normalization (within-sample)**  

Within-sample normalization aims to normalise counts across cells which can differ due to sequencing depth, RNA content, and efficiency of lysis and reverse transcription [(Saket Choudhary & Rahul Satija, 2022; ](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02584-9)[Grün, Kester and van Oudenaarden. 2014)](http://scholar.google.com/scholar_lookup?&title=Validation%20of%20noise%20models%20for%20single-cell%20transcriptomics&journal=Nat%20Methods&volume=11&issue=6&pages=637-40&publication_year=2014&author=Grün%2CD&author=Kester%2CL&author=van%20Oudenaarden%2CA). Some methods normalize based on the total expression detected per-cell, however this approach is more appropriate for bulk RNA-seq as it assumes that each cell started with the same number of RNA molecules (e.g **Seurat** `NormalizeData` which normalizes by the total expression, multiplied by a scale factor and log-transformed). More recent methods apply downsampling or statistical models. Tools include Linnorm [(Yip et al., 2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7019105/#B42), SCnorm [(Bacher et al., 2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7019105/#B3), scran [(Lun et al., 2016)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7019105/#B25), and more recently, Normalisr [(Wang, 2021)](https://www.nature.com/articles/s41467-021-26682-1), sctransform [(Hafemeister and Satija, 2019)](https://doi.org/10.1186%2Fs13059-019-1874-1), bayNorm [(Tang et al. 2020)](https://www.nature.com/articles/s41467-021-26682-1#ref-CR15), and Sanity [(Breda et al. 2019, *preprint*)](https://www.biorxiv.org/content/10.1101/2019.12.28.889956v1). Other normalization may be considered depending on the specific study design. For example, normalising for biological covariates (such as cell cycle stage) may be useful for trajectory inference. See discussion in [(Leucken and Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746).

Others?? *tbc* **Variance stabilization**. The aim here is to correct for the relationship between gene expression and variation in expression (a well-known effect which is corrected for in bulk RNA-seq pipelines, for example by DESeq2. Between-sample normalization, imputation

> **Empty droplets and ambient gene expression**  

A proportion of the barcodes in the count matrix will correspond to empty droplets. A simple method to identify such droplets may use a minimum threshold of counts. However, this risks removing cells with low levels of gene expression, such as quiescent cells. More sophisticated methods identify empty droplets through *ambient gene expression*, which refers to the background expression detected in all droplets from contaminating cell-free RNA. Modelling ambient gene expression in empty droplets can also allow other droplets to be corrected for contaminating counts, thus avoiding problems in downstream analysis and clustering. Several sophisticated methods exist with which to remove empty droplets and correct the remaining data for ambient or 'background' gene expression. These include SoupX [(Young and Behjati, 2020)](https://doi.org/10.1093/gigascience/giaa151), EmptyDrops [(Lun et al. 2019)](https://doi.org/10.1186%2Fs13059-019-1662-y) and DEIM [(Alvarez et al. 2020](https://www.nature.com/articles/s41598-020-67513-5) (note DEIM was designed for single-nuclei data, although shown to work with scRNA as well). This pipeline will use **SoupX**. 

> **Filtering low-quality and uninformative barcodes**  

GEMs containing damaged cells can be identified using commonly-used metrics: **(1) the number of read counts per barcode, (2) the number of genes detected per barcode & (3) the proportion of mitochondrial DNA.** For example, very high mtDNA expression can indicate cells where the cytoplasmic RNA leaked from the cell membrane. Hard cut-offs may be used to filter the data, however this may remove meaningful biology, such as quiesent cells with low gene counts or highly active cells with a high level of mitochondria gene expression. Prior knowledge of the expected cell types may inform the use of thresholds. 

**QC thresholds should be specific to your data**. For example, [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) retained cells with gene counts \>200 and \<2,500 and a mitochondrial content \<10%. [Xu et al. (2021)](https://academic.oup.com/hmg/article/30/5/370/6131713?login=false#supplementary-data) retained cells with >30,000 UMIs, 200-5,000 genes, and less than 50% mitochondrial expression. It is challenging to determine optimal QC thresholds *a priori*, particularly in samples with heterogenous cell-types; thus it is advisable to first apply permisive threshold and then revist QC steps one or more times to optimise the distribution of QC metrics and clusters (*not* to alter or improve the results and outcome statistics, such as differential gene expression, which falls into the trap of 'cherry picking' or 'data peeking'!) [(Leucken and Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746). It is also advisable to carry out QC filtering seperately for independent samples, due to differences in sample quality [(Plasschaert et al. 2018; ](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6108322/)[Chen, Ning and Shi. 2019)](https://doi.org/10.3389/fgene.2019.00317).

> **Doublets**   

Doublets have been removed in the past by simply removing cells with particularly high gene counts (e.g. with significant deviation). More sophisticated methods typically identify doublets based on their similarity with simulated doublets. Currently available tools include: doubletCells (Lun et al., 2016), Scrublet [(Wolock et al., 2019)](https://www.sciencedirect.com/science/article/pii/S2405471220301952#bib32), cxds (Bais and Kostka, 2020), bcds (Bais and Kostka, 2020), hybrid (Bais and Kostka, 2020), Solo [(Bernstein et al. 2020)](https://doi.org/10.1016/j.cels.2020.05.010), DoubletDetection (Gayoso and Shor, 2018), DoubletFinder ([McGinnis et al., 2019a](https://www.sciencedirect.com/science/article/pii/S2405471220304592#bib44), [2019b](https://www.sciencedirect.com/science/article/pii/S2405471220304592#bib45)), and DoubletDecon (DePasquale et al., 2019). These methods were recently compared by [Xi and Li (2021)](https://www.sciencedirect.com/science/article/pii/S2405471220304592), who report **DoubletFinder** and **Solo** as the top two performing methods. Briefly, **DoubletFinder** uses a *k*-nearest neighbors (kNN) algorithm to identify doublets based on their clustering with simulated doublets in principal component space. **Solo** (included in the scvi-tools suite from the Yosef Lab at UC Berkeley) uses a semi-supervised deep learning approach and claims improvements over DoubletFinder by not assuming linear gene expression circuits (note [Xi and Li (2021)](https://www.sciencedirect.com/science/article/pii/S2405471220304592) reported DoubletFinder as the top method).


### Downstream analysis

> **Feature selection, dimensionality reduction and visualisation**    

Including *clustering and cell identity*
*TBC*

#### Dataset integration

*TBC*

#### Differential expression and pseudobulk expression  
*TBC*

## 3. Workspace setup

This pipeline will process data from Chromium Single Cell A Chip Kit (10X Genomics) using Seurat v4 and Cell Ranger v6.1.

A conda environment for the analysis can be created as follows:
```bash
conda create -n scRNA -c conda-forge r-base r-essentials
conda install -n scRNA -c r r-irkernel
conda install -n scRNA -c bioconda r-seurat
conda install -c r r-devtools
#conda install -n scRNA2 -c conda-forge scvi-tools
```

You can install the latest version of Cell Ranger. To install Cell Ranger, you will need to register at [this link](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest). 

Create a folder to save the new CellRanger code, `cd` and run the `wget` or `curl` command provided following the registration. Then:

```bash
tar -zxvf cellranger-6.1.2.tar.gz
```
Next add the path with the CellRanger executable to your PATH. **NOTE: you will have to run this command every time you want to use CellRanger**.
```bash
export PATH=/rds/general/user/hm1412/home/anaconda3/envs/scRNA2/bin/cellranger-6.1.2:$PATH
```

## 4. Generate a reference transcriptome

The first steps will follow the recommended pipeline from [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ov).

`cellranger count` will require a transcriptome reference. You can download the prebuilt reference. However, here we will compile our own using our preferred versions of the reference genome and gene annotation:

```bash
# Genome metadata
genome="GRCh38"
version="2021-Mar"

# Set up source and build directories
build="GRCh38-gencode_v36_build"
mkdir -p "$build"

fasta_in=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
gtf_in=gencode.v36.annotation.gtf

# Create reference package
cellranger mkref --ref-version="$version" --genome="$genome" --fasta="$fasta_in" --genes="$gtf_in"
```

## 5. Cellranger count

To run `cellranger count`, make sure your files are in the `bcl2fastq` naming convention e.g. `SRR10009414_S1_L00X_R1_001.fastq.gz` (and the corresponding `I1` and `R2`). The below command should be run, where `<ID>` is the sample ID at the start of the filename (e.g. SRR10009414) and the `<PATH>` should direct to the reference directory created by the previous command. **Technical replicates should be combined here**.

```bash
#Run cellranger count with the sampleID and cellranger reference directory
cellranger count --id <ID> --transcriptome <PATH>

#If working with public data i.e. pre-computed clusters:
cellranger count --nosecondary --id <ID> --transcriptome <PATH>

#Example for donor SAMN12614700, with four sets of fastq files from four runs in the SAMN12614700 directory:
#cellranger count --nosecondary --id SAMN12614700 --sample SRR10009414,SRR10009415,SRR10009416,SRR10009417 --transcriptome $GENOMEDIR/GRCh38 --fastqs SAMN12614700/
``` 

## 5. Secondary analysis

An overview of the secondary analysis steps:

1. Read data into R
2. Pre-processing and preliminary clustering
2. `SoupX` to correct for ambient gene expression (vignette available [here](https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html))
4. QC and removal of outlying clusters
5. Repeat pre-processing and clustering
3. `DoubletFinder` to identify and remove doublets
4. Final round of processing and clustering 


First, load the required R libraries:

```R
#Install SoupX using CRAN
#install.packages("SoupX")
#Or the latest developers version directly from Github
#devtools::install_github("constantAmateur/SoupX",ref='devel')
library(SoupX)
library(purrr)
library(data.table)
library(SoupX)
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(ggplot2)
library(limma)
library(DoubletFinder)
```

Read in the data filtered by CellRanger, here shown for the example donor SAMN12614700 from [Ramachandran et al. (2019)](https://www.nature.com/articles/s41586-019-1631-3), which includes four technical replicates. All barcodes which passed the cellranger filtering will be included in the initial pre-processing, in order to generate preliminary clusters for input to SoupX. 

```R
#Read in data for four technical replicates
#Initialize the Seurat object with the raw (non-normalized data).
SAMN12614700.data=Read10X("SAMN12614700_male_healthy/SAMN12614700/outs/filtered_feature_bc_matrix/")
SAMN12614700 <- CreateSeuratObject(counts = SAMN12614700.data, project = "SAMN12614700", 
                                          min.cells = 1, min.features = 1)

#If technical replicates are being merged at this stage instead, the following code can be used:
#Technical replicate 1
#rep1.data=Read10X("SAMN12614700_male_healthy/SRR10009414_control/outs/filtered_feature_bc_matrix/")
#rep1 <- CreateSeuratObject(counts = SRR10009414.data, project = "SAMN12614700", min.cells = 1, min.features = 1)
#Technical replicate 2
#rep2.data=Read10X("SAMN12614700_male_healthy/SRR10009415_control/outs/filtered_feature_bc_matrix/")
#rep2 <- CreateSeuratObject(counts = SRR10009415.data, project = "SAMN12614700", min.cells = 1, min.features = 1)
#To combine technical replicates:
#merged <- merge(rep1, y = rep2,add.cell.ids = c("rep1", "rep2"),project = "SAMN12614700")

```

The output shows that there is `26276 features across 5560 samples within 1 assay`, meaning 26,276 expressed genes and 5,560 cells. 

The following round of pre-processing will include normalization, clustering and marker gene identification. The normalization will be carried out using [`SCTransform`](https://satijalab.org/seurat/articles/sctransform_vignette.html), which normalizes using a negative binominal in place of a scale factor. SCTransform has been reported to recover improved biological meaning compared to log-normalization based on a scale factor. 

```R
#SCTransform can be used in place of the NormalizeData, FindVariableFeatures, ScaleData workflow.
SAMN12614700 <- SCTransform(SAMN12614700, conserve.memory=TRUE,return.only.var.genes=TRUE)

#Carry out dimensionality reduction and clustering:
SAMN12614700 <- RunPCA(object = SAMN12614700, verbose = FALSE)
SAMN12614700 <- RunUMAP(object = SAMN12614700, dims = 1:20, verbose = FALSE)
SAMN12614700 <- FindNeighbors(object = SAMN12614700, dims = 1:20, verbose = FALSE)
SAMN12614700 <- FindClusters(object = SAMN12614700, verbose = FALSE)
```

Identify cluster marker genes:
```R
#Identify marker genes
SAMN12614700.markers <- FindAllMarkers(SAMN12614700, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Subset the top 5 marker genes per cluster
markers=SAMN12614700.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
```

These initial clusters can be visualised in a UMAP plot. Note that these clusters form the temporary input to the next stage of analysis (`SoupX`) and are *not* the final clusters.

```R
#UMAP plot
DimPlot(object = SAMN12614700, label = TRUE, reduction = "umap") + NoLegend() + ggtitle("preliminary clustering - sctransform")
```

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/new_UMAP1.png" width="400">

Add the UMAP coordinates to the metaData for compatability with SoupX:

```R
#Add the UMAP coordinates to the metaData for later compatability with SoupX
umapCoord <- as.data.frame(Embeddings(object = SAMN12614700[["umap"]]))

#You can then add the PC columns to the meta.data as you please (cbind, sapply, paste(..,collapse=),...)
#For example, if you only need to add the PCs components, without modifying them or doing other processing, you can do a simple merge
metaData <- SAMN12614700@meta.data

metaData$UMI_id <- rownames(metaData)
umapCoord$UMI_id <- rownames(umapCoord)

metaData <- merge.data.frame(metaData,umapCoord,by = "UMI_id")
metaData$UMI_id <- NULL
rownames(metaData) = rownames(SAMN12614700@meta.data)

SAMN12614700@meta.data <- metaData
head(SAMN12614700@meta.data)
```

The UMAP coordinates for each barcode are in the UMAP_1 and UMAP_2 columns:
<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/new_metaData.png" width="800">

### SoupX

`SoupX` will be used to correct ambient gene expression. SoupX will read the data from CellRanger, including both the gene count matrices for the raw data (all droplets) and the filtered droplets. The raw data is used to estimate the contaminating gene counts, the *ambient gene expression*. The 'soup' model is created as part of the `load10X` command and the raw counts are removed afterwards. SoupX will be run on each of the four technical replicates prior to merging the background-corrected gene count matrices:

```R
#Read in data for SoupX
sc = load10X("SAMN12614700_male_healthy/SAMN12614700/outs/")

#Add the barcode clusters to the soupChannel object
sc = setClusters(sc, setNames(metaData[colnames(sc$toc),]$seurat_clusters, colnames(sc$toc)))

#Set the UMAP coordinates:
sc = setDR(sc, metaData[colnames(sc$toc), c("UMAP_1", "UMAP_2")])

#Estimate the contaminating fraction
sc = autoEstCont(sc) 
# Clean the data
out = adjustCounts(sc)
```

The output data includes the genes whose contaminating reads make a significant contribution to the "soup", i.e. contribute the *ambient gene expression*. These will most likely include several mitochondrial genes.

```R
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
plotMarkerDistribution(sc)
```

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/soup.png" height="400">

The UMAP plot can be used to explore the ambient gene expression. For example, you can highlight cells where a gene is expressed *at all* (i.e. at least one count detected) vs cells where the gene is significantly expressed *above* the background. We can explore this below for the highly expressed hepatocyte gene, *ALB*. (Note that this example data is from non-parenchymal cells, or NPCs, and should not contain hepatocytes).

```R
#Plot the expression of an example gene
dd$ALB = sc$toc["ALB", ] 
ggplot(dd, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = ALB > 0))
plotMarkerMap(sc, "ALB",DR=sc$metaData[,c('UMAP_1','UMAP_2')])
```

Cells where ALB is expressed at >1 count:  
<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/ALB.png" width="400">

Cells where ALB has significant expression, above that of the background "soup" contamination in all cells:  
<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/ALB2.png" width="400">

Create a Seurat data object from the four corrected replicates. This step will involve some initial filtering, to excluded genes detected in less than 3 cells and to exclude cells with less than 200 detected genes. 

### QC and filtering 

```R 
SAMN12614700.noSoup <- CreateSeuratObject(counts = out, project = "SAMN12614700",  min.cells = 3, min.features = 200)
SAMN12614700.noSoup[["percent.mt"]] <- PercentageFeatureSet(SAMN12614700.noSoup, pattern = "^MT-")
#Look at the QC distribution
VlnPlot(SAMN12614700, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/new_QC_plots1.png"> 

```R
#Remove droplets with %mtDNA>50
SAMN12614700.noSoup <- subset(SAMN12614700.noSoup, subset = percent.mt < 50)
#Repeat the processing
SAMN12614700.noSoup <- SCTransform(SAMN12614700.noSoup, conserve.memory=TRUE,return.only.var.genes=TRUE)
SAMN12614700.noSoup <- RunPCA(object = SAMN12614700.noSoup, verbose = FALSE)
SAMN12614700.noSoup <- RunUMAP(object = SAMN12614700.noSoup, dims = 1:20, verbose = FALSE)
SAMN12614700.noSoup <- FindNeighbors(object = SAMN12614700.noSoup, dims = 1:20, verbose = FALSE)
SAMN12614700.noSoup <- FindClusters(object = SAMN12614700.noSoup, verbose = FALSE)
#Plot the QC measures for the new clusters
VlnPlot(SAMN12614700.noSoup, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/new_QC_plots2.png">

Here, we can already see a suspicous cluster, cluster 2, which has a higher distribution of %mtDNA expression and lower number of counts and features than other clusters. We repeat the processing and define marker genes, then we can check what marker genes have been defined for cluster 2. 

```R
#Find cluster markers
SAMN12614700.markers.noSoup <- FindAllMarkers(SAMN12614700.noSoup, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers.noSoup=SAMN12614700.markers.noSoup %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

#What are the marker genes of cluster 2?
markers.noSoup[markers.noSoup$cluster==2,]
```

We can see that the marker genes for cluster 2 are mitochondrial genes, a clear indication that this cluster contains low-quality, possibly ruptured, cells.  

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/cluster2.png" height="200" >

We remove cluster 2 and repeat the processing. Notably, we can also see some other suspicious clusters, such as cluster 13. However, the marker genes for cluster 13 are more informative (*ALB*, *APOC3*, *APOA1*) and so we will keep it (for now!). Cluster 13 also has a higher distribution of counts so may contain a higher number of doublets, i.e. drolets with two (or more!) cells. These will be investigated in the next step.

```R
#Remove cluster 2
SAMN12614700.filtered=subset(SAMN12614700.noSoup,idents=2,invert=TRUE)
#Repeat the processing
SAMN12614700.filtered <- SCTransform(SAMN12614700.filtered, conserve.memory=TRUE,return.only.var.genes=TRUE)
SAMN12614700.filtered <- RunPCA(object = SAMN12614700.filtered, verbose = FALSE)
SAMN12614700.filtered <- RunUMAP(object = SAMN12614700.filtered, dims = 1:20, verbose = FALSE)
SAMN12614700.filtered <- FindNeighbors(object = SAMN12614700.filtered, dims = 1:20, verbose = FALSE)
SAMN12614700.filtered <- FindClusters(object = SAMN12614700.filtered, verbose = FALSE)
```

## Remove doublets

Next, doublets will be identified and removed by `doubletFinder`. (Note, you can check the number of informative dimensions with an ElbowPlot `ElbowPlot(SRR10009414_control.noSoup.SCT)`).

```R
## pK Identification (no ground-truth) 
sweep.res.list_SAMN12614700 <- paramSweep_v3(SAMN12614700.filtered, sct = TRUE, PCs = 1:15)

## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(SAMN12614700.filtered@meta.data$seurat_clusters)  ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(SAMN12614700.filtered@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies 
SAMN12614700.filtered <- doubletFinder_v3(SAMN12614700.filtered, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

#Get the column name from the metaData e.g. DF.classifications_0.25_0.09_90
head(SAMN12614700.filtered@meta.data) 

#How many doublets?
table(SAMN12614700.filtered@meta.data$DF.classifications_0.25_0.09_90)

#Where are the doublets?
DimPlot(SAMN12614700.filtered, group.by = "DF.classifications_0.25_0.09_90", pt.size = 0.01, cols = c("red", "azure3"))

# check the nUMI for doublet and singlet
VlnPlot(SAMN12614700.filtered,
        features = "nCount_RNA",
        pt.size = 0,
        group.by = "DF.classifications_0.25_0.09_90") + NoLegend()
```

Which cells were marked as doublets?  

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/doublets_UMAP.png" height="500" >

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/doublets_counts.png" height="400" >

This code from the scrublet vignette can show how many doublets were present per-cluster:
```R
#https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/scrublet-doublet-validation.html
df <- data.table(SAMN12614700.filtered@meta.data)
sel.meta <- c("DF.classifications_0.25_0.09_90", "seurat_clusters", "orig.ident")
df <- df[, sel.meta, with = FALSE]

df[, 2:3] %>% map( ~ {
  freq1 <- df[, .N, keyby = .(.x, DF.classifications_0.25_0.09_90)]
  freq1[, total := sum(N), by = .(.x)]
  freq1[, ratio := N / total]
  
  linesize = .35
  fontsize = 8

  ggplot(freq1, aes(fill=DF.classifications_0.25_0.09_90, y=ratio, x= .x)) + 
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values = c("Doublet" = 'red', "Singlet" = "grey")) +
    xlab('Cluster') +
    scale_y_continuous(breaks = seq(0,1,0.1), expand = c(0,0), name = 'Percentage')+
    theme_bw()+
    theme( panel.grid.major.x = element_blank(), 
           panel.grid.major.y = element_blank(),
           panel.grid.minor = element_blank(),
           strip.background = element_blank(),panel.border = element_rect(size = linesize),
           axis.ticks = element_blank(), 
           axis.text.x = element_text(size = 5))
  
})
```

We can see here that, of the new clusters, cluster 6 was made up of almost 75% doublets:  

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/doublets_clusters.png" height="400" >

Remove doublets and repeat the clustering:

```R
#Remove the doublets, keeping only Singlets
SAMN12614700.filtered=subset(x = SAMN12614700.filtered, subset = DF.classifications_0.25_0.09_90 == 'Singlet')

#Repeat the processing with SCTransform, dimensionality reduction and clustering.
#This time, SCTransform will be run in the full mode, i.e will return all genes 
SAMN12614700.filtered <- SCTransform(SAMN12614700.filtered, conserve.memory=FALSE)
SAMN12614700.filtered <- RunPCA(object = SAMN12614700.filtered, verbose = FALSE)
SAMN12614700.filtered <- RunUMAP(object = SAMN12614700.filtered, dims = 1:15, verbose = FALSE)
SAMN12614700.filtered <- FindNeighbors(object = SAMN12614700.filtered, dims = 1:15, verbose = FALSE)
SAMN12614700.filtered <- FindClusters(object = SAMN12614700.filtered, verbose = FALSE)
DimPlot(object = SAMN12614700.filtered, label = TRUE, reduction = "umap") + NoLegend() + ggtitle("sctransform. noSoup. noDoublet.")
```

To check how many cells are in each cluster:
```R
#How many cells are in each cluster?
table(Idents(object = SAMN12614700.filtered))
#Explore the clusters and QC metrics of the filtered data
DimPlot(object = SAMN12614700.filtered, label = TRUE, reduction = "umap") + NoLegend() + ggtitle("sctransform")
VlnPlot(SAMN12614700.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)
```

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/filtered_UMAP.png" height="400" >  

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/filtered_violin.png">


Any threshold for filtering genes should be informed by your experimental design, including the number of cells in the dataset and the number of cells in the smallest cluster of interest [(Leucken and Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746).
* * * * * * * * * * * * * * * * * * * * * * * * * *

## 6. Cell type identity

Next, the clusters will be identified as specific cell-types based on curated marker gene sets. Listed below are several methods used by published papers to annotate clusters of liver scRNA-seq data:

- [Ramachandran et al. (2019)](https://www.nature.com/articles/s41586-019-1631-3): removed contaminating circulatory cells based on clustering with scRNA from PBMCs. Obtained a signature score across a curated  list of known marker genes per cell lineage in the liver. The score was defined as the mean expression of the signature marker genes. 
- [MacParland et al. (2018)](https://www.nature.com/articles/s41467-018-06318-7): "the cell-type identities for each cluster were determined manually using a compiled panel of available known hepatocyte/immune cell transcripts."
- [Wang et al. (2021)](https://www.nature.com/articles/s41598-021-98806-y): "we compared out manual annotation of clusters based on known cell-type-specific marker genes to that produced through automated classification using [SingleR](https://www.nature.com/articles/s41598-021-98806-y#ref-CR33)".
- [Payen et al. (2021)](https://www.sciencedirect.com/science/article/pii/S2589555921000549#sec2) "The expression of different combinations of genes was used to define scores and signatures using the Seurat PercentageFeatureSet function".
- [Aizarani et al. (2019)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6687507/). cells from select clusters were reanalyzed with [RaceID3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6687507/#R4) (see paper for parameters). [StemID](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6687507/#R24) was run on resulting, filtered and feature-selected expression matrix, with target clusters "inferred by FateID using ASGR1 plus ALB and CXCL8 plus MMP7 as markers for hepatocyte and cholangiocyte lineage target clusters. Using the KRT19 and CFTR as mature cholangiocyte markers yields highly similar results."

Here, we will calculate a cell-type score based on curated marker gene sets obtained by integrating marker genes from PanglaoDB and in-house curated gene sets. The curated gene set lists are available as a file in this Github (xxxxxx), with the pipeline used to generate this file described below:

```sh
#Download the PangloaDB markers 2020 file
wget https://www.panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz 

#Extract markers for liver cell types
awk -v FS="\t" '{if(($1=="Mm Hs" || $1=="Hs") && $3~/(Hepatocytes|Macrophages|Cholangiocytes|Endothelial cells|Hepatic stellate cells|Hepatoblasts|Hepatocytes|Kupffer cells|Macrophages|Mast cells|NK cells|T cells|B cells|Monocytes|Myoblasts|Myocytes)/) print $0}' PanglaoDB_markers_27_Mar_2020.tsv | cut -f 2,3 | grep -Ev "blood|aorta" > liver.markers
```

Here, we will use Seurat `PercentageFeatureSet`:

*"This function enables you to easily calculate the percentage of all the counts belonging to a subset of the possible features for each cell. This is useful when trying to compute the percentage of transcripts that map to mitochondrial genes for example. The calculation here is simply the column sum of the matrix present in the counts slot for features belonging to the set divided by the column sum for all features times 100."*

```R
markers=read.table('liver.markers',sep='\t')
LSEC.markers=read.table('LSEC.markers')
#subset the marker genes dataframe for the genes reported
markers=subset(markers,markers[,1] %in% rownames(SAMN12614700.filtered@assays$SCT@counts))
LSEC.markers=subset(LSEC.markers,LSEC.markers[,1] %in% rownames(SAMN12614700.filtered@assays$SCT@counts))
#Combine the database and in-house LSEC marker list
markers=rbind(markers,cbind(V1=LSEC.markers,V2='LSEC'))

#PercentageFeatureSet
for(x in unique(markers[,2])){
    name=print(gsub(' ','.',x)) #replace spaces in the name with .
    features=as.character(subset(markers,markers[,2]==x)[,1])
    SAMN12614700.filtered[[name]]=PercentageFeatureSet(SAMN12614700.filtered,features = features)
}
```

You can colour the UMAP plot according to the features, which are the specific cell-type % expression scores stored in the metaData.

```R
FeaturePlot(object = SAMN12614700.filtered, features = "Endothelial.cells") #unique(markers[,2]))
```

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/Endothelial_score_UMAP.png" height="300">

We can also plot a heatmap of the cell-type scores by creating a matrix with the mean score per cluster:

```R
metaData=SAMN12614700.filtered@meta.data
cell.types=gsub(' ','.',unique(markers[,2]))
cell_type.scores=aggregate(metaData[,cell.types], list(metaData$seurat_clusters), mean)

heatmap(as.matrix(cell_type.scores[,-1]),scale="row",margins=c(10,6))
```


`Seurat PercentageFeatureSet`

## Integrating biological replicates

Also calculated each cell lineage, calculated Pearson correlation across replicates).

# References

- [Leucken and Theis (2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746): A 2019 effort to compile current best practises in scRNA-seq. This paper is very useful for "newbies" and gives an excellent overview of the essential steps of scRNA-seq analysis (note that some specific tools mentioned are superceeded by more recent published tools).
- https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

e.g. [Xu et al. (2021)](https://academic.oup.com/hmg/article/30/5/370/6131713?login=false#supplementary-data) use `DoubletFinder` and then `decontX` of `celda` to correct for "cross containment".
3. **CellRanger**: cellranger also provides it's own tools for secondary analysis, the "`cellranger reanalyze` command reruns secondary analysis performed on the feature-barcode matrix (dimensionality reduction, clustering and visualization) using different parameter settings."

# Extra notes

[`DropletUtils`](https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html) is an R package for filtering droplet-based scRNA-seq data. (See example use in [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741)). This package can be used to read in the count matrix to R, filter empty droplets and visualise QC plots for the distribution of counts across barcodes, as well as downsample the count matrix or raw reads etc.


 **Recommendations**: *do not* run `DoubletFinder` on aggregated scRNA-seq data representing multiple distinct samples and filter your data first to remove low-quality cell clusters. The authors recommends (1) manually threshold raw gene expression matrices according to RNA nUMIs, (2) pre-process data using [standard workflow](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html), including data normalization and scaling (3) identify clusters with low RNA UMIs, high % mitochondrial reads and/or uninformative marker genes, (4) remove clusters, pre-process again, and run DoubletFinder.



Note...sctransform was used by the Satija lab in their 2022 [publication](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02584-9#Sec8)

"In particular, two recent studies proposed to use generalized linear models (GLMs), where cellular sequencing depth was included as a covariate, as part of scRNA-seq preprocessing workflows. Our sctransform [9] approach utilizes the Pearson residuals from negative binomial regression as input to standard dimensional reduction techniques, while GLM-PCA [10] focuses on a generalized version of principal component analysis (PCA) for data with Poisson-distributed errors. More broadly, multiple techniques aim to learn a latent state that captures biologically relevant cellular heterogeneity using either matrix factorization or neural networks [11–13], alongside a defined error model that describes the variation that is not captured by the latent space.



- [(Brüning et al. 2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) filter cells using the R packages DropletUtils and then use `Seurat` for downstream analysis, retaining cells with gene counts \>200 and \<2,500 and a mitochondrial content \<10%.
- [Xu et al. (2021)](https://academic.oup.com/hmg/article/30/5/370/6131713?login=false#supplementary-data) "removed doublet cells using “DoubletFinder”. Then we use the function “decontX” of “celda” to correct the probable cross containment. Then the corrected expression matrix was processed by “Seurat V3.14”. For the quality control, cells with 0~30000 UMIs, 200~5000 genes, and less than 50% mitochondrial expression percentage were filtered out for the next analysis."

For single-nuclei RNA, [Hardwick et al. (2022)](https://www.nature.com/articles/s41587-022-01231-3) excluded nuclei with unique gene counts >7,500 or <200 or >4% mitochondrial gene expression. UMI numbers and mitochondrial gene expression % were regressed from each nucleus and the matrix was log-normalised and scaled to 10,000 reads per cell. Performed both tSNE and UMAP non-lnear reduction techniques... cell types assigned by canonical marker genes for each cluster... cell type annotation confirmed by aligning to <other data>. (Hardwick et al. 2022, Nature Biotechnology).

7. **Feature selection, dimensionality reduction and visualisation**. 

> Feature selection  
Here, the data is filtered for informative genes, such as *highly variable genes* (HVGs), usually between 1,000 and 5,000. Tools for selecting HVGs are provided by Seurat and Scanpy, which bin genes by their mean expression and select genes with the highest variance-to-mean ratio. Make sure to check what type of data your method expects (e.g. raw count data or log-transformated data). For your final analysis, HVGs should be selected *after* normalization and pre-processing. (Note some tools select HVGs as part of the pre-processing step; these are not the same as your final list of HVGs which give your clusters).

For dimensionaly reduction, UMAP is reported to be the optimal method over several other popular alternatives (e.g. t-SNE, PCA and MDS) [(Yang et al. 2021)](https://www.sciencedirect.com/science/article/pii/S2211124721008597)

> Dimensionality reduction  
**Clustering**. In the 2019 review by Leuken and Theis, the most popular method for clustering was **multi-resolution modularity optimization** algorithm as implemented in the Louvian algorithm. This algorithm detects groups of cells that have "more links between them than expected from the number of links the cells have in total". (Implemented by Seurat and Scanpy). 

![Figure 1 (Leuken and Theis, 2019)](https://github.com/CebolaLab/scRNA/blob/main/Figures/Leucken_Theis_Table1.png)
Table 1 from [(Leucken and Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746) summarises the types of data expected as input for downstream analysis.

**DATASET INTEGRATION** 

Integrating multiple scRNA-seq datasets presents an additional challenge, which may again be tackled with different methods. Note, the below are specifically for integrating multiple sequencing runs of different GEM Wells. For samples sequenced in the same GEM well, pass the multiple fastq files to `cellranger count` using the `--fastqs` argument. 

- [(Brüning et al. 2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) integrate expression matrices using Suerat, including: normalization with the `SCTransform` function, rank the features using the `SelectIntegrationFeatures` function, with the resulting features controlled using the function `PrepSCTIntegration`. Anchors were determined by `FindIntegrationAnchors` and afterwards used with the `IntegrateData` function. 
- With **CellRanger** using `cellranger aggr` which aggregates multiple runs of `cellranger count`, normalizes runs to the same effective sequencing depth (the pipeline equalizes the average read depth per cell between groups before merging), and then performs secondary analysis on the combined data (this produces one single-feature barcode matrix and a .cloupe file for visualizing with Loupe Browser). `cellranger aggr` uses Chemistry Batch Correction when aggregating resuts from a combination of 5' and 3', or 3' v2 and 3' v3 Gene Expression data, which improves the mixing of the batches in the t-SNE visualization and clustering results (note that residual batch effects may still be present). 

Seurat anchors?

- Seurat `FindClusters`. 
- CellRanger provides the `Loupe Browser` which can be used to explore data.

#### UMAP plots

[(Brüning et al. 2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) UMAP run with the first 20 PCs.
**Seurat**: `FindAllMarkers`
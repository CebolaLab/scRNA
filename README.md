# Analysis of 10X Genomics scRNA-seq data
This Github describes the pipeline used by the Cebola Lab to analyse single-cell RNA-seq data (scRNA-seq) from 10X Genomics. 

## 1. Introduction

This pipeline has been developed by carefully reviewing current tools and best practises used in the analysis of 10X Genomics scRNA-seq data, as of March 2022. This Github will first present an overview of various available tools, followed by the Cebola Lab pipeline. Resources used are shown in the [References](#references) at the bottom of this page.

#### Alignment, demultiplexing and quantification

There are several alignment algorithms to choose from, including CellRanger, STARsolo, Alevin, Alevin-fry and Kallisto; these are compared in [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741). This pipeline will use the 10X Genomics toolbox, **CellRanger**. If memory requirement is an issue, [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) suggest STARsolo as an alternative.

#### Secondary analysis

CellRanger includes the demultiplexing of reads from individual droplets, or 'GEMs' (Gel Beads in EMulsion), based on unique barcodes. These should correspond to unique single cells *in theory*. However, droplets can contain more or less than one cell, or contain damaged or low quality cells. Therefore, the secondary analysis will assess the quality of the dataset and filter the data to retain only droplets with high-quality single cells. 

The pipeline for secondary analysis, including references, is discussed below.

1. **Empty droplets and ambient gene expression**. A proportion of the barcodes in the count matrix will correspond to empty droplets. A simple method to identify such droplets may use a minimum threshold of counts. However, this has the disadvantage of removing potential informative clusters of cells with low gene-expression levels, such as quiescent cells. More sophisticated methods make use of *ambient gene expression*, which refers to the background expression detected in all droplets, including empty droplets, from contaminating cell-free RNA. Modelling ambient gene expression in empty droplets can also allow other droplets to be corrected for contaminating counts, thus avoiding problems in downstream analysis and clustering. Several sophisticated methods exist with which to remove empty droplets and correct the remaining data for ambient or 'background' gene expression. Tools include **SoupX** [(Young and Behjati, 2020)](https://doi.org/10.1093/gigascience/giaa151), **EmptyDrops** [(Lun et al. 2019)](https://doi.org/10.1186%2Fs13059-019-1662-y) and **DEIM** [(Alvarez et al. 2020](https://www.nature.com/articles/s41598-020-67513-5) (note DEIM was designed for single-nuclei data, although shown to work with scRNA as well). This pipeline will use SoupX. Briefly, SoupX quantifies cell-free RNA contamination from the profiles of empty droplets and outputs background-corrected counts for the remaining GEMs.

2. **Outliers**. Common metrics used to filter uninformative GEMs such as damaged cells include the number of read counts and genes, as well as the % of mitochondrial DNA expression. For example, very high mtDNA expression can indicate cells where the cytoplasmic RNA leaked from the cell membrane, while gene expressed in little to no cells (gene-level filtering), will remove ininformative data. Earlier methods set hard cut-off thresholds, however this may remove meaningful biology, such as quiesent cells with low gene counts or highly active cells with a high level of mitochondria gene expression. 

An improved approach may be to identify clusters of cells with 

If setting cut-offs It is advisable to explore your data to set 

Any threshold for filtering genes should be informed by your experimental design, including the number of cells in the dataset and the number of cells in the smallest cluster of interest [(Leucken and Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746).

Filtering typically remove barcodes based on three typical QC measures: (1) the number of read counts per barcode, (2) the number of genes detected per barcode, and (3) the proportion of mitochondrial DNA.

The thresholds will be specific to the distribution of reads in your data, with the aim of removing barcodes with outliers. For example, [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) retained cells with gene counts \>200 and \<2,500 and a mitochondrial content \<10%. [Xu et al. (2021)](https://academic.oup.com/hmg/article/30/5/370/6131713?login=false#supplementary-data) retained cells with >30000 UMIs, 200-5000 genes, and less than 50% mitochondrial expression percentage. [Hardwick et al. (2022)](https://www.nature.com/articles/s41587-022-01231-3) excluded nuclei with unique gene counts >7,500 or <200 or >4% mitochondrial gene expression (note this is single-**nuclei** RNA-seq). 

It is challenging to determine optimal QC thresholds *a priori*, particularly in samples with heterogenous cell-types; thus it may be advisable to first apply permisive threshold and then revist QC steps one or more times to optimise the distribution of QC metrics and clusters (*not* to alter or improve the results and outcome statistics, such as differential gene expression, which falls into the trap of 'cherry picking' or 'data peeking'!) [(Leucken and Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746). It is also advisable to carry out QC filtering seperately for independent samples, due to differences in sample quality [(Plasschaert et al. 2018; ](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6108322/)[Chen, Ning and Shi. 2019)](https://doi.org/10.3389/fgene.2019.00317).

4. **Doublets**. Doublets have been removed in the past by simply removing cells with particularly high gene counts (e.g. with significant deviation). More sophisticated methods typically identify doublets based on their similarity with simulated doublets. Currently available tools include: doubletCells (Lun et al., 2016), Scrublet [(Wolock et al., 2019)](https://www.sciencedirect.com/science/article/pii/S2405471220301952#bib32), cxds (Bais and Kostka, 2020), bcds (Bais and Kostka, 2020), hybrid (Bais and Kostka, 2020), Solo [(Bernstein et al. 2020)](https://doi.org/10.1016/j.cels.2020.05.010), DoubletDetection (Gayoso and Shor, 2018), DoubletFinder ([McGinnis et al., 2019a](https://www.sciencedirect.com/science/article/pii/S2405471220304592#bib44), [2019b](https://www.sciencedirect.com/science/article/pii/S2405471220304592#bib45)), and DoubletDecon (DePasquale et al., 2019). These methods were recently compared by [Xi and Li (2021)](https://www.sciencedirect.com/science/article/pii/S2405471220304592), who report **DoubletFinder** and **Solo** as the top two performing methods. Briefly, **DoubletFinder** uses a *k*-nearest neighbors (kNN) algorithm to identify doublets based on their clustering with simulated doublets in principal component space. **Solo** (included in the scvi-tools suite from the Yosef Lab at UC Berkeley) uses a semi-supervised deep learning approach and claims improvements over DoubletFinder by not assuming linear gene expression circuits (note [Xi and Li (2021)](https://www.sciencedirect.com/science/article/pii/S2405471220304592) reported DoubletFinder as the top method).

5. **Normalization (within-sample)**. 
Within-sample normalization aims to normalise counts across cells which can differ due to sequencing depth, RNA content, and efficiency of lysis and reverse transcription [(Saket Choudhary & Rahul Satija, 2022; ](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02584-9)[Grün, Kester and van Oudenaarden. 2014)](http://scholar.google.com/scholar_lookup?&title=Validation%20of%20noise%20models%20for%20single-cell%20transcriptomics&journal=Nat%20Methods&volume=11&issue=6&pages=637-40&publication_year=2014&author=Grün%2CD&author=Kester%2CL&author=van%20Oudenaarden%2CA). Some methods normalize based on the total expression detected per-cell, however this approach is more appropriate for bulk RNA-seq as it assumes that each cell started with the same number of RNA molecules (e.g **Seurat** `NormalizeData` which normalizes by the total expression, multiplied by a scale factor and log-transformed). More recent methods apply downsampling or statistical models. Tools include Linnorm [(Yip et al., 2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7019105/#B42), SCnorm [(Bacher et al., 2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7019105/#B3), scran [(Lun et al., 2016)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7019105/#B25), and more recently, Normalisr [(Wang, 2021)](https://www.nature.com/articles/s41467-021-26682-1), sctransform [(Hafemeister and Satija, 2019)](https://doi.org/10.1186%2Fs13059-019-1874-1), bayNorm [(Tang et al. 2020)](https://www.nature.com/articles/s41467-021-26682-1#ref-CR15), and Sanity [(Breda et al. 2019, *preprint*)](https://www.biorxiv.org/content/10.1101/2019.12.28.889956v1). Other normalization may be considered depending on the specific study design. For example, normalising for biological covariates (such as cell cycle stage) may be useful for trajectory inference. See discussion in [(Leucken and Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746).

**Variance stabilization**. The aim here is to correct for the relationship between gene expression and variation in expression (a well-known effect which is corrected for in bulk RNA-seq pipelines, for example by DESeq2).

5. **Normalization (between-sample)**. 

6. Imputation?

#### Downstream analysis

7. **Feature selection, dimensionality reduction and visualisation**. 

> Feature selection  
Here, the data is filtered for informative genes, such as *highly variable genes* (HVGs), usually between 1,000 and 5,000. Tools for selecting HVGs are provided by Seurat and Scanpy, which bin genes by their mean expression and select genes with the highest variance-to-mean ratio. Make sure to check what type of data your method expects (e.g. raw count data or log-transformated data). For your final analysis, HVGs should be selected *after* normalization and pre-processing. (Note some tools select HVGs as part of the pre-processing step; these are not the same as your final list of HVGs which give your clusters).

For dimensionaly reduction, UMAP is reported to be the optimal method over several other popular alternatives (e.g. t-SNE, PCA and MDS) [(Yang et al. 2021)](https://www.sciencedirect.com/science/article/pii/S2211124721008597)

> Dimensionality reduction  
**Clustering**. In the 2019 review by Leuken and Theis, the most popular method for clustering was **multi-resolution modularity optimization** algorithm as implemented in the Louvian algorithm. This algorithm detects groups of cells that have "more links between them than expected from the number of links the cells have in total". (Implemented by Seurat and Scanpy). 

![Figure 1 (Leuken and Theis, 2019)](https://github.com/CebolaLab/scRNA/blob/main/Figures/Leucken_Theis_Table1.png)
Table 1 from [(Leucken and Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746) summarises the types of data expected as input for downstream analysis.

#### Dataset integration

Integrating multiple scRNA-seq datasets presents an additional challenge, which may again be tackled with different methods. Note, the below are specifically for integrating multiple sequencing runs of different GEM Wells. For samples sequenced in the same GEM well, pass the multiple fastq files to `cellranger count` using the `--fastqs` argument. 

- [(Brüning et al. 2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) integrate expression matrices using Suerat, including: normalization with the `SCTransform` function, rank the features using the `SelectIntegrationFeatures` function, with the resulting features controlled using the function `PrepSCTIntegration`. Anchors were determined by `FindIntegrationAnchors` and afterwards used with the `IntegrateData` function. 
- With **CellRanger** using `cellranger aggr` which aggregates multiple runs of `cellranger count`, normalizes runs to the same effective sequencing depth (the pipeline equalizes the average read depth per cell between groups before merging), and then performs secondary analysis on the combined data (this produces one single-feature barcode matrix and a .cloupe file for visualizing with Loupe Browser). `cellranger aggr` uses Chemistry Batch Correction when aggregating resuts from a combination of 5' and 3', or 3' v2 and 3' v3 Gene Expression data, which improves the mixing of the batches in the t-SNE visualization and clustering results (note that residual batch effects may still be present). 

Seurat anchors?

#### Clustering and assigning cell types

- Seurat `FindClusters`. 
- CellRanger provides the `Loupe Browser` which can be used to explore data.

#### UMAP plots

[(Brüning et al. 2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) UMAP run with the first 20 PCs

#### Differential expression

**Seurat**: `FindAllMarkers`

## 2. Installations

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

Create a folder to save the new CellRanger code, `cd` and run the `wget` or `curl` command provided following the registration. 

```bash
tar -zxvf cellranger-6.1.2.tar.gz
```
Next add the path with the CellRanger executable to your PATH. **NOTE: you will have to run this command every time you want to use CellRanger**.
```bash
export PATH=/rds/general/user/hm1412/home/anaconda3/envs/scRNA2/bin/cellranger-6.1.2:$PATH
```

## 2. Generate a reference transcriptome

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

## 3. Cellranger count

To run `cellranger count`, make sure your files are in the `bcl2fastq` naming convention e.g. `SRR10009414_S1_L00X_R1_001.fastq.gz` (and the corresponding `I1` and `R2`). The below command should be run, where `<ID>` is the sample ID at the start of the filename (e.g. SRR10009414) and the `<PATH>` should direct to the reference directory created by the previous command. **Technical replicates can be combined here, or in the next stage in R.**

```bash
#Run cellranger count with the sampleID and cellranger reference directory
cellranger count --id <ID> --transcriptome <PATH>

#If working with public data i.e. pre-computed clusters:
cellranger count --nosecondary --id <ID> --transcriptome <PATH>

#Example for donor SAMN12614700, with four sets of fastq files from four runs in the SAMN12614700 directory:
#cellranger count --nosecondary --id SAMN12614700 --sample SRR10009414,SRR10009415,SRR10009416,SRR10009417 --transcriptome $GENOMEDIR/GRCh38 --fastqs SAMN12614700/
``` 

## 4. Secondary analysis

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

Read in the data filtered by CellRanger, here shown for the example sample SRR10009414 from [Ramachandran et al. (2019)](https://www.nature.com/articles/s41586-019-1631-3). All the cells will initially be included, in order to carry out normalization and preliminary clustering required for SoupX.

```R
#Read in data for Seurat
SAMN12614700.data=Read10X("SAMN12614700/outs/filtered_feature_bc_matrix/")

#Initialize the Seurat object with the raw (non-normalized data).
#No filtering at this stage
SAMN12614700 <- CreateSeuratObject(counts = SAMN12614700.data, project = "SAMN12614700", min.cells = 1, min.features = 1)

#To combine technical replicates:
#SAMN12614700.data2=Read10X("SAMN12614700_2/outs/filtered_feature_bc_matrix/")
#SAMN12614700.2 <- CreateSeuratObject(counts = SAMN12614700.data2, project = "SAMN12614700", min.cells = 1, min.features = 1)

#SAMN12614700 <- merge(SAMN12614700, y = SAMN12614700.2, add.cell.ids = c("rep1", "rep1", project = "SAMN12614700")
#SAMN12614700
```

The output shows that there is `#24175 features across 1392 samples within 1 assay`, meaning 24,175 expressed genes and 1,392 cells. 

Next, an initial round of pre-processing will include normalization, clustering and marker gene identification. The normalization will be carried out using [`SCTransform`](https://satijalab.org/seurat/articles/sctransform_vignette.html), which normalizes using a negative binominal in place of a scale factor. SCTransform has been reported to recover improved biological meaning compared to log-normalization based on a scale factor. 

```R
#SCTransform can be used in place of the NormalizeData, FindVariableFeatures, ScaleData workflow.
SRR10009414_control.SCT <- SCTransform(SRR10009414_control, conserve.memory=TRUE,return.only.var.genes=TRUE)
#Carry out dimensionality reduction and clustering:
SRR10009414_control.SCT <- RunPCA(object = SRR10009414_control.SCT, verbose = FALSE)
SRR10009414_control.SCT <- RunUMAP(object = SRR10009414_control.SCT, dims = 1:20, verbose = FALSE)
SRR10009414_control.SCT <- FindNeighbors(object = SRR10009414_control.SCT, dims = 1:20, verbose = FALSE)
SRR10009414_control.SCT <- FindClusters(object = SRR10009414_control.SCT, verbose = FALSE)
```

Identify cluster marker genes:
```R
SRR10009414.markers <- FindAllMarkers(SRR10009414_control.SCT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers=SRR10009414.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
```

These initial clusters can be visualised in a UMAP plot. Note that these clusters form the temporary input to the next stage of analysis (`SoupX`) and are *not* the final clusters.

```R
#UMAP plot
DimPlot(object = SRR10009414_control.SCT, label = TRUE, reduction = "umap") + NoLegend() + ggtitle("sctransform")
```

![UMAP plot 1](https://github.com/CebolaLab/scRNA/blob/main/Figures/UMAP1.png)

We will add the UMAP coordinates to the metaData so they can be read by SoupX:
```R
#Add the UMAP coordinates to the metaData for later compatability with SoupX
umapCoord <- as.data.frame(Embeddings(object = SRR10009414_control.SCT[["umap"]]))

#You can then add the PC columns to the meta.data as you please (cbind, sapply, paste(..,collapse=),...)
metaData <- SRR10009414_control.SCT@meta.data

metaData$UMI_id <- rownames(metaData)
umapCoord$UMI_id <- rownames(umapCoord)

metaData <- merge.data.frame(metaData,umapCoord,by = "UMI_id")
metaData$UMI_id <- NULL
rownames(metaData) = rownames(SRR10009414_control.SCT@meta.data)

#Add the new dataframe to the object as metaData
SRR10009414_control.SCT@meta.data <- metaData
```

### SoupX

`SoupX` will be used to correct ambient gene expression. SoupX will read the data from CellRanger, including both the gene count matrices for the raw data (all droplets) and the filtered droplets. The raw data is used to estimate the contaminating gene counts, the *ambient gene expression*. The 'soup' model is created as part of the `load10X` command. 

```R
#Read in data for SoupX
sc = load10X("SRR10009414_control/outs/")
#sc contains droplets estimated by cellRanger to contain cells.
#Add the previously defined clusters and UMAP coordinates to this new data object.
sc = setClusters(sc, setNames(metaData[rownames(sc$metaData),]$seurat_clusters, rownames(sc$metaData)))
sc = setDR(sc, metaData[colnames(sc$toc), c("UMAP_1", "UMAP_2")])

#Set "top" marker genes (here they are ordered by adjusted p-value)
top.markers=SRR10009414.markers[order(SRR10009414.markers$p_val_adj),]
sc = autoEstCont(sc,topMarkers = top.markers)
# Estimate rho
sc = autoEstCont(sc)
# Clean the data
out = adjustCounts(sc)
```

The output data includes the genes whose contaminating reads make a significant contribution to the "soup", i.e. contribute the *ambient gene expression*. These will most likely include several mitochondrial genes.

```R
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
plotMarkerDistribution(sc)
```

You can highlight cells within the UMAP plot where a gene is expressed *at all* vs expressed above the background. For example for the mitochondrial gene `MT-CO1`. 

```R
#Plot the expression of an example gene
dd$CO1 = sc$toc["MT-CO1", ] #FCN1 is kupffer, CFTR is Cholangiocyte
ggplot(dd, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = CO1 > 0))
#CLEC4G PECAM1 LYVE1
plotMarkerMap(sc, "MT-CO1",DR=sc$metaData[,c('UMAP_1','UMAP_2')])
```

Cells where MT-CO1 is expressed at >1 count:  
![MT-CO1_all](https://github.com/CebolaLab/scRNA/blob/main/Figures/MT-CO1.1.png)

Cells where MT-CO1 has significant expression, above that of the background "soup" contamination in all cells:  
![MT-CO1_significant](https://github.com/CebolaLab/scRNA/blob/main/Figures/MT-CO1.2.png)

Correct the count matrix for the ambient gene expresison and create a Seurat data object.

```R 
out = adjustCounts(sc)
```

## Quality-control and filtering

Using the new, background-corrected count matrix (`srat`), we will now explore the QC of the data and filter out low-quality cells and/or clusters. First, the initial pre-processing will be rerun:

```R
SRR10009414_control.noSoup <- CreateSeuratObject(counts = out, project = "SRR10009414_control", min.cells = 3, min.features = 100)
SRR10009414_control.noSoup
```

We now have `18033 features across 1279 samples within 1 assay`.

```R
#SCTransform can be used in place of the NormalizeData, FindVariableFeatures, ScaleData workflow.
SRR10009414_control.noSoup.SCT <- SCTransform(SRR10009414_control.noSoup, conserve.memory=TRUE,return.only.var.genes=TRUE)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
SRR10009414_control.noSoup.SCT[["percent.mt"]] <- PercentageFeatureSet(SRR10009414_control.noSoup.SCT, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(SRR10009414_control.noSoup.SCT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

![QC1](https://github.com/CebolaLab/scRNA/blob/main/Figures/QC1.png)

Here, the data can be filtered to remove outliers. This pipeline will preferentially filter the data using clustering later on, but here we will remove cells with %mtDNA > 50%.

```R
#Remove droplets with %mtDNA>50
SRR10009414_control.noSoup.SCT <- subset(SRR10009414_control.noSoup.SCT, subset = percent.mt < 50)
#Rerun dimensionality reduction
SRR10009414_control.noSoup.SCT <- RunPCA(object = SRR10009414_control.noSoup.SCT, verbose = FALSE)
SRR10009414_control.noSoup.SCT <- RunUMAP(object = SRR10009414_control.noSoup.SCT, dims = 1:20, verbose = FALSE)
SRR10009414_control.noSoup.SCT <- FindNeighbors(object = SRR10009414_control.noSoup.SCT, dims = 1:20, verbose = FALSE)
SRR10009414_control.noSoup.SCT <- FindClusters(object = SRR10009414_control.noSoup.SCT, verbose = FALSE)
```

## Remove doublets

Next, doublets will be identified and removed by `doubletFinder`. (Note, you can check the number of informative dimensions with an ElbowPlot `ElbowPlot(SRR10009414_control.noSoup.SCT)`).

```R
## pK Identification (no ground-truth) 
sweep.res.list_SRR10009414 <- paramSweep_v3(SRR10009414.NoSoup.no7.SCT, sct = TRUE, PCs = 1:20)

## Homotypic Doublet Proportion Estimate 
homotypic.prop <- modelHomotypic(SRR10009414_control.noSoup.SCT@meta.data$seurat_clusters)  ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(SRR10009414_control.noSoup.SCT@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies 
SRR10009414_control.noSoup.SCT <- doubletFinder_v3(SRR10009414_control.noSoup.SCT, PCs = 1:16, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

#How many doublets?
table(SRR10009414_control.noSoup.SCT@meta.data$DF.classifications_0.25_0.09_90)

#Where are the doublets?
DimPlot(SRR10009414_control.noSoup.SCT, group.by = "DF.classifications_0.25_0.09_90", pt.size = 0.01, cols = c("red", "azure3"))

# check the nUMI for doublet and singlet
VlnPlot(SRR10009414_control.noSoup.SCT,
        features = "nCount_RNA",
        pt.size = 0,
        group.by = "DF.classifications_0.25_0.09_90") + NoLegend()
```

This code from the scrublet vignette can show how many doublets were present per-cluster:
```R
#https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/scrublet-doublet-validation.html
df <- data.table(SRR10009414_control.noSoup.SCT@meta.data)
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

Remove doublets and repeat the clustering:

```R
#Remove the doublets, keeping only Singlets
SRR10009414_control.noSoup.Singlet.SCT=subset(x = SRR10009414_control.noSoup.SCT, subset = DF.classifications_0.25_0.09_90 == 'Singlet')
#Repeat the dimensionality reduction
SRR10009414_control.noSoup.Singlet.SCT <- RunPCA(object = SRR10009414_control.noSoup.Singlet.SCT, verbose = FALSE)
SRR10009414_control.noSoup.Singlet.SCT <- RunUMAP(object = SRR10009414_control.noSoup.Singlet.SCT, dims = 1:20, verbose = FALSE)
SRR10009414_control.noSoup.Singlet.SCT <- FindNeighbors(object = SRR10009414_control.noSoup.Singlet.SCT, dims = 1:20, verbose = FALSE)
SRR10009414_control.noSoup.Singlet.SCT <- FindClusters(object = SRR10009414_control.noSoup.Singlet.SCT, verbose = FALSE)
plot1 <- DimPlot(object = SRR10009414_control.noSoup.Singlet.SCT, label = TRUE, reduction = "umap") + NoLegend() + ggtitle("sctransform. noSoup. noDoublet.")
plot1
```


To check how many cells are in each cluster:
```R
#How many cells are in each cluster?
table(Idents(object = SRR10009414_control.noSoup.SCT))
```

* * * * * * * * * * * * * * * * * * * * * * * * * *

`DropletUtils`


```R
#Recommended pipeline by the author of DoubletFinder
## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
seu_kidney <- CreateSeuratObject(kidney.data)
seu_kidney <- NormalizeData(seu_kidney)
seu_kidney <- FindVariableFeatures(seu_kidney, selection.method = "vst", nfeatures = 2000)
seu_kidney <- ScaleData(seu_kidney)
seu_kidney <- RunPCA(seu_kidney)
seu_kidney <- RunUMAP(seu_kidney, dims = 1:10)

## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
seu_kidney <- CreateSeuratObject(kidney.data)
seu_kidney <- SCTransform(seu_kidney)
seu_kidney <- RunPCA(seu_kidney)
seu_kidney <- RunUMAP(seu_kidney, dims = 1:10)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
gt.calls <- seu_kidney@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"].   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_kidney@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
```

- [(Brüning et al. 2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) filter cells using the R packages DropletUtils and then use `Seurat` for downstream analysis, retaining cells with gene counts \>200 and \<2,500 and a mitochondrial content \<10%.
- [Xu et al. (2021)](https://academic.oup.com/hmg/article/30/5/370/6131713?login=false#supplementary-data) "removed doublet cells using “DoubletFinder”. Then we use the function “decontX” of “celda” to correct the probable cross containment. Then the corrected expression matrix was processed by “Seurat V3.14”. For the quality control, cells with 0~30000 UMIs, 200~5000 genes, and less than 50% mitochondrial expression percentage were filtered out for the next analysis."

For single-nuclei RNA, [Hardwick et al. (2022)](https://www.nature.com/articles/s41587-022-01231-3) excluded nuclei with unique gene counts >7,500 or <200 or >4% mitochondrial gene expression. UMI numbers and mitochondrial gene expression % were regressed from each nucleus and the matrix was log-normalised and scaled to 10,000 reads per cell. Performed both tSNE and UMAP non-lnear reduction techniques... cell types assigned by canonical marker genes for each cluster... cell type annotation confirmed by aligning to <other data>. (Hardwick et al. 2022, Nature Biotechnology).

# References

- [Leucken and Theis (2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746): A 2019 effort to compile current best practises in scRNA-seq. This paper is very useful for "newbies" and gives an excellent overview of the essential steps of scRNA-seq analysis (note that some specific tools mentioned are superceeded by more recent published tools).
- https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

e.g. [Xu et al. (2021)](https://academic.oup.com/hmg/article/30/5/370/6131713?login=false#supplementary-data) use `DoubletFinder` and then `decontX` of `celda` to correct for "cross containment".
3. **CellRanger**: cellranger also provides it's own tools for secondary analysis, the "`cellranger reanalyze` command reruns secondary analysis performed on the feature-barcode matrix (dimensionality reduction, clustering and visualization) using different parameter settings."


[`DropletUtils`](https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html) is an R package for filtering droplet-based scRNA-seq data. (See example use in [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741)). This package can be used to read in the count matrix to R, filter empty droplets and visualise QC plots for the distribution of counts across barcodes, as well as downsample the count matrix or raw reads etc.


 **Recommendations**: *do not* run `DoubletFinder` on aggregated scRNA-seq data representing multiple distinct samples and filter your data first to remove low-quality cell clusters. The authors recommends (1) manually threshold raw gene expression matrices according to RNA nUMIs, (2) pre-process data using [standard workflow](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html), including data normalization and scaling (3) identify clusters with low RNA UMIs, high % mitochondrial reads and/or uninformative marker genes, (4) remove clusters, pre-process again, and run DoubletFinder.



Note...sctransform was used by the Satija lab in their 2022 [publication](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02584-9#Sec8)

"In particular, two recent studies proposed to use generalized linear models (GLMs), where cellular sequencing depth was included as a covariate, as part of scRNA-seq preprocessing workflows. Our sctransform [9] approach utilizes the Pearson residuals from negative binomial regression as input to standard dimensional reduction techniques, while GLM-PCA [10] focuses on a generalized version of principal component analysis (PCA) for data with Poisson-distributed errors. More broadly, multiple techniques aim to learn a latent state that captures biologically relevant cellular heterogeneity using either matrix factorization or neural networks [11–13], alongside a defined error model that describes the variation that is not captured by the latent space.
# Indivudal dataset analysis

The following analysis will be carried out in R. I recommend to run a Jupyter notebook with a conda environment (see the recommended installations below).

An overview of the secondary analysis steps:

1. [Read data into R](#1-read-data-into-r)
2. [Pre-processing and preliminary clustering](#2-Pre-processing-and-preliminary-clustering)
3. [Correct for ambient gene expression using SoupX](#3-Correct-for-ambient-gene-expression-using-SoupX) (vignette available [here](https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html))
4. [QC and removal of outlying clusters](#4-QC-and-removal-of-outlying-clusters)
5. [Remove doublets with DoubletFinder](#5-Remove-doublets-with-DoubletFinder)
6. [Final round of processing and clustering](#6-Final-round-of-processing-and-clustering)

A conda environment for the analysis can be created as follows:

```bash
conda create -n scRNA -c conda-forge r-base r-essentials
conda install -n scRNA -c r r-irkernel
conda install -n scRNA -c bioconda r-seurat
conda install -n scRNA -c r r-devtools
conda install -n scRNA -c conda-forge scvi-tools
```

## 1. Read data into R
Once you have opened R (e.g. in a Jupyter notebook), load the required R libraries:

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
#Initialize the Seurat object with the filtered count data.
SAMN12614700.data=Read10X("SAMN12614700_male_healthy/SAMN12614700/outs/filtered_feature_bc_matrix/")
SAMN12614700 <- CreateSeuratObject(counts = SAMN12614700.data, project = "SAMN12614700", 
                                          min.cells = 1, min.features = 1)

#If technical replicates are being merged at this stage, the following code can be used to merge Seurat objects:
#merged <- merge(rep1, y = rep2, add.cell.ids = c("rep1", "rep2"),project = "SAMN12614700").
```

The output shows that there is `26276 features across 5560 samples within 1 assay`, meaning 26,276 expressed genes and 5,560 cells. 

## 2. Pre-processing and preliminary clustering

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

## 3. Correct for ambient gene expression using SoupX

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

## 4. QC and removal of outlying clusters

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

## 5. Remove doublets with DoubletFinder

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

Remove doublets:

```R
#Remove the doublets, keeping only Singlets
SAMN12614700.filtered=subset(x = SAMN12614700.filtered, subset = DF.classifications_0.25_0.09_90 == 'Singlet')
```

## 6. Final round of processing and clustering 

```bash
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

Next, see the [Donor integration tutorial](https://github.com/CebolaLab/scRNA/tree/main/4.%20Donor%20integration).
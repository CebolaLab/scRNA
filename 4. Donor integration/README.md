
## 8. Donor integration

Up to now, biological replicates have been processed independently. This provides an opportunity to assess replicate similarity including by calculating Pearson correlation for cell-type clusters across replicates, and to remove low-quality cells from each donor.

Following this [Seurat tutorial](https://satijalab.org/seurat/articles/integration_introduction.html), we will [integrate biological replicates which have been processed and normalised using SCTranform](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) with `SelectIntegrationFeatures` > `PrepSCTIntegration()` > `FindIntegrationAnchors()` > `IntegrateData()` (use the normalization.method = "SCT" option where appropriate). See also the [Fast integration using reciprocal PCA (RPCA) vignette](https://satijalab.org/seurat/articles/integration_rpca.html). The steps will follow:

- Create a list of Seurat objects to integrate
- Integrate data using `FindIntegrationAnchors`
- Repeat dimensionality reduction (PCA/UMAP) and clustering 

First, load (and if necessary, install) the required libraries:
```R
#Load libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
library(ape)
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
library(harmony)
#install.packages('ggbeeswarm')
library(ggbeeswarm)
#install.packages('ggthemes')
library(ggthemes)
#BiocManager::install("destiny")
library(destiny)
library(scater)
#library(SingleCellExperiment)
#BiocManager::install("scater")
library(cowplot)
```

```R
#Create a list with the Seurat objects for the donors
liver.list=list(SAMN12614699.renamed,SAMN12614700.renamed,SAMN12614710.clusters,SAMN12614716.filtered,SAMN12614717.filtered)

#Note you may want to increase the memory availability: options(future.globals.maxSize = 8000 * 1024^2)
liver.features <- SelectIntegrationFeatures(object.list = liver.list, nfeatures = 3000) #3000
#Prepare suerat object, ensures all appropriate residuals have been calculated
liver.list <- PrepSCTIntegration(object.list = liver.list, anchor.features = liver.features, 
    verbose = FALSE)

#Find anchors
liver.anchors <- FindIntegrationAnchors(object.list = liver.list, normalization.method = "SCT", 
    anchor.features = liver.features, verbose = FALSE)
#Intergate data
liver.integrated <- IntegrateData(anchorset = liver.anchors, normalization.method = "SCT", 
    verbose = FALSE)
```

Note the strength of integration can be increased by increasing the `k.anchor` parameter, in the `FindIntegrationAnchors` command.

An alternative workflow for **large datasets** employs [reference-based integration](https://satijalab.org/seurat/articles/integration_large_datasets.html). Here, one (or two, if you wish for one male and one female) dataset is used as reference to reduce the number of comparisons, rather than identifying anchors between al pairs of query datasets. This involves using the `reference` option (e.g. `reference = c(1, 2)` in the `FindIntegrationAnchors` command. 

As before, dimensionality reduction and clustering will be run on the now integrated data. 

```R
liver.integrated <- RunPCA(liver.integrated, verbose = FALSE)
liver.integrated <- RunUMAP(liver.integrated, dims = 1:20)
liver.integrated <- FindNeighbors(object = liver.integrated, dims = 1:20, verbose = FALSE)
liver.integrated <- FindClusters(object = liver.integrated, verbose = FALSE)
```

### Visualise UMAP and donor integration

You can now view the UMAP and colour the cells by original donor, to access how successful the integration has been. You may want to test different integration methods and select the one which shows the best mixing of the donors (i.e.the data is not clustering by donor, but by cell type).

```R
#Dimensionality reduction plot coloured by both original identiti (donor) and seurat clusters (new following donor integrati)
plots <- DimPlot(liver.integrated, group.by = c("orig.ident", "seurat_clusters")) 
plots & theme(legend.position = "top", ) & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
    override.aes = list(size = 3))) 
```

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/first_integrated_UMAP.png" height="600">

The next step is to read in a list of defined marker genes, here obtained from PangloaDB:

```R
#Upload cell types markers and create a merged dataframe 
markers=read.table('liver.markers',sep='\t')
LSEC.markers=read.table('LSEC.markers')
markers=subset(markers,markers[,1] %in% rownames(liver.integrated@assays$SCT@counts))
LSEC.markers=subset(LSEC.markers,LSEC.markers[,1] %in% rownames(liver.integrated@assays$SCT@counts))
markers=rbind(markers,cbind(V1=LSEC.markers,V2='LSEC'))

#Using PercentageFeatureSet, assign to each cell the % gene expression from each set of marker genes
#Assign to the metadata, for each cell, the % of total gene expression which is from the gene set. 
for(x in unique(markers[,2])){
    name=gsub(' ','.',x)
    features=as.character(subset(markers,markers[,2]==x)[,1])
    #class(features)
    liver.integrated[[name]]<-PercentageFeatureSet(liver.integrated,features = features, assay = 'RNA')
}

```

Next, see the [pseudobulk count and bigwig visualisation integration tutorial](https://github.com/CebolaLab/scRNA/tree/main/9.pseudobulk_counts_bigwigs).


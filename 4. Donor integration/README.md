
## 4. Donor integration

Up to now, biological replicates have been processed independently. This provides an opportunity to assess replicate similarity including by calculating Pearson correlation for cell-type clusters across replicates, and to remove low-quality cells from each donor.

Following this [Seurat tutorial](https://satijalab.org/seurat/articles/integration_introduction.html), we will [integrate biological replicates which have been processed and normalised using SCTranform](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) with `SelectIntegrationFeatures` > `PrepSCTIntegration()` > `FindIntegrationAnchors()` > `IntegrateData()` (use the normalization.method = "SCT" option where appropriate). See also the [Fast integration using reciprocal PCA (RPCA) vignette](https://satijalab.org/seurat/articles/integration_rpca.html). The steps will follow:

1. [Integrate donor data](#Step-one---Integrate-donor-data)
2. [Dimensionality reduction (PCA/UMAP) and clustering](#Step-2---dimensionality-reduction-and-clustering)
3. [Marker gene expression and cluster annotation](#Step-3---Marker-gene-expression-and-cluster-annotation)

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

### Step one - Integrate donor data

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

### Step 2 - Dimensionality reduction and clustering

```R
liver.integrated <- RunPCA(liver.integrated, verbose = FALSE)
liver.integrated <- RunUMAP(liver.integrated, dims = 1:20)
liver.integrated <- FindNeighbors(object = liver.integrated, dims = 1:20, verbose = FALSE)
liver.integrated <- FindClusters(object = liver.integrated, verbose = FALSE)
```

You can now view the UMAP and colour the cells by original donor, to access how successful the integration has been. You may want to test different integration methods and select the one which shows the best mixing of the donors (i.e.the data is not clustering by donor, but by cell type).

```R
#Dimensionality reduction plot coloured by both original identiti (donor) and seurat clusters (new following donor integrati)
plots <- DimPlot(liver.integrated, group.by = c("orig.ident", "seurat_clusters")) 
plots & theme(legend.position = "top", ) & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
    override.aes = list(size = 3))) 
```

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/first_integrated_UMAP.png" height="600">

The next step is to read in a list of defined marker genes, here obtained from PangloaDB:

### Step 3 - Marker gene expression and cluster annotation 

These following steps will:

1. Read in a list of marker genes (gene sets) for liver cell types
2. Calculate the % expression from the genes in each gene set per-cell
3. Identify marker genes for each Seurat cluster

```R
#Upload cell types markers and create a merged dataframe 
markers=read.table('liver.markers',sep='\t')
LSEC.markers=read.table('LSEC.markers')
markers=subset(markers,markers[,1] %in% rownames(liver.integrated@assays$SCT@counts))
LSEC.markers=subset(LSEC.markers,LSEC.markers[,1] %in% rownames(liver.integrated@assays$SCT@counts))
markers=rbind(markers,cbind(V1=LSEC.markers,V2='LSEC'))

#Using PercentageFeatureSet, assign to the metadata for each cell the % gene expression 
#from each set of marker genes
for(x in unique(markers[,2])){
    name=gsub(' ','.',x)
    features=as.character(subset(markers,markers[,2]==x)[,1])
    liver.integrated[[name]]<-PercentageFeatureSet(liver.integrated,features = features, assay = 'RNA')
}

#Calculte the marker genes for each cluster
#Find markers for every cluster compared to all remaining cells, report only the positive ones
liver.cluster.markers <- FindAllMarkers(liver.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers=liver.cluster.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
```

Check `cluster.markers` and if any cluster show mitochondrial marker genes, these clusters should be removed and [step two](#Step-2-dimensionality-reduction-and-clustering) repeated. 


Optional additional data exploration:
```R
liver.integrated=BuildClusterTree(liver.integrated)
Tool(object = liver.integrated, slot = 'BuildClusterTree')
PlotClusterTree(object = liver.integrated)
```

You can now colour the UMAP by the % expression for the gene sets, as assigned above. Adjust the `max.cutoff` parameter for the colour to be more or less intensive.
```R
FeaturePlot(object = liver.integrated, features = "Kupffer.cells",label=TRUE, max.cutoff=20) 
FeaturePlot(object = liver.integrated, features = "Endothelial.cells",label=TRUE)
FeaturePlot(object = liver.integrated, features = "Cholangiocytes",label=TRUE, max.cutoff = 5) 
FeaturePlot(object = liver.integrated, features = "Hepatic.stellate.cells", label=TRUE, max.cutoff=10)
```

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/cell_colour_UMAP.png" height="600">

Clusters can be renamed as shown:

```R
liver.integrated <- RenameIdents(liver.integrated, '5'='Stellate cells','9'='Stellate cells',
                                 '11'='Cholangiocyte','12'='Cholangiocyte',
                                 '0'='Endothelial','1'='Endothelial','2'='Endothelial','3'='Endothelial','8'='Endothelial',
                                 '10'='Endothelial','13'='Endothelial','14'='Endothelial','20'='Endothelial')
```

Another way of viewing which clusters score most highly in terms of gene expression from the gene sets is via heatmap:

```R
metaData=liver.integrated@meta.data
cell.types=gsub(' ','.',unique(markers[,2]))
cell_type.scores=aggregate(metaData[,cell.types], list(Idents(liver.integrated)), mean)
heatmap(as.matrix(cell_type.scores[,-1]),scale="column",margins=c(10,6),labRow=cell_type.scores[,1])
```

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/cell_colour_heatmap.png" height="600">

Once your clusters have been labelled, (1) save these identities in a file and (2) create pseudobulk count data which will reflect cell-type specific expression levels.

```R
identity=as.data.frame(Idents(liver.integrated))
write.table(identity,'cell_identity.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=FALSE)

pseudocount=AggregateExpression(object=liver.integrated,slot="counts",assays="RNA")
write.table(pseudocount$RNA,'pseudo_counts.txt',sep='\t',quote=FALSE,col.names=TRUE,row.names=TRUE)
```

Next, further clustering can be carried out. In this case, we will identify liver sinusoidal endothelial cells (LSECs) from the endothelial cluster, then further cluster periportal, pericentral and central venous LSECs using pseudotime and diffusion analysis. 


```R
Endothelial=subset(liver.integrated, idents = "Endothelial")
#SCTransform can be used in place of the NormalizeData, FindVariableFeatures, ScaleData workflow.
Endothelial <- SCTransform(Endothelial, conserve.memory=TRUE,return.only.var.genes=FALSE)
all.genes <- rownames(Endothelial)
Endothelial <- RunPCA(Endothelial, features = all.genes)
Endothelial <- RunUMAP(Endothelial, dims = 1:20)
Endothelial <- FindNeighbors(object = Endothelial, dims = 1:20, verbose = FALSE)
Endothelial <- FindClusters(object = Endothelial, verbose = FALSE)
```

Create a list of known marker genes for zonated LSECs:
```R
#Central venous LSECs
LSEC.central.venous.extended=c('FCGR2B','STAB2','LYVE1','CD14','ICAM1',"FCN3","FCN2") #ICAM1 = CD54,FCGR2B = CD32B,"THBD" 
#Periportal LSECs
#Note periportal express PECAM1, but so do non-LSECs
LSEC.periportal.extended=c('PECAM1','F8',"SPARCL1","CLEC14A") #"DLL4","MSR1","LTBP4","NTN4") #CD31 = PECAM1 PECAM1, MRS1
nonLSECs=c('VWF',"CD34","ENG","ACKR1",'PECAM1') #CD34?  ACKR1','CD34','PECAM1
#CD36
LSEC.mid=c('LYVE1','CTSL')
LSEC.pericentral='KIT'

endo.markers=rbind(cbind(V1=nonLSECs,V2='nonLSECs'),
                   cbind(V1=LSEC.central.venous.extended,V2="LSEC.central.venous.extended"),
                   cbind(V1=LSEC.periportal.extended, V2="LSEC.periportal.extended"),
                   cbind(V1=LSEC.mid, V2='LSEC.mid'),
                   cbind(V1=LSEC.pericentral,V2='LSEC.pericentral'))

#Again add the % gene expression to the metadata
#?PercentageFeatureSet
for(x in unique(endo.markers[,2])){
    name=gsub(' ','.',x)
    features=as.character(subset(endo.markers,endo.markers[,2]==x)[,1])
    #class(features)
    Endothelial[[name]]<-PercentageFeatureSet(Endothelial,features = features, assay = 'RNA')
}

FeaturePlot(object = Endothelial, features = c("LSEC","nonLSECs","LSEC.central.venous.extended",
                                               "LSEC.periportal.extended"),
            label=FALSE, max.cutoff = c(20,0.6,4,1)) 
```

Note that the peripotal LSEC markers have some overlapping expression with non-LSEC endothelial cells, so the periportal markers are likely to also "light up" the non-LSEC clusters. The cluster with periportal LSECs is likely to have the "nonLSECs" markers showing little to no blue but a positive score from the periportal markers. 

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/Endo_colour_UMAP.png" height="400">

Clusters which don't score positively with any of the endothelial sets and cluster seperately may indicate carry over of other cells types. In this case, repeat the steps to determine per-cluster marker genes. Here, it may be useful to search the [Human Protein Atlas](https://www.proteinatlas.org/ENSG00000164825-DEFB1/single+cell+type) (HPA) to explore the expression of the identified marker genes in published scRNA-seq data from liver. For example:

```R
#Find markers for every cluster compared to all remaining cells, report only the positive ones
Endo.markers <- FindAllMarkers(Endothelial, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers=Endo.markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
```

Here, clusters 9, 11 and 13 were not clearly assigned to any of the endothelial populations. The marker genes identified by `FindAllMarkers` are:

- Cluster 9: DEFB1, TM4SF4, SERPINA1, ANXA4, KRT8
- Cluster 11: TAGLN, ACTA2, MYL9, RGS5, TPM2
- Cluster 13: CCL21, TFF3, FABP4, NTS, FABP5

According to the HPA, the genes for cluster 9 are highly expressed in cholangiocytes, while the cluster 11 genes are highly expressed in stellate cells. The cluster 13 genes are less specific, but are consistently expressed in stellate cells.

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/Endo_UMAP.png" height="400">

Again, save the cell identities and pseudocount data before extracting the LSECs for further clustering:

```R
identity.endo=as.data.frame(Idents(Endothelial))
write.table(identity.endo,'cell_Endothelial_identity.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=FALSE)

#head(liver.integrated@assays$RNA@counts)
pseudocount.endo=AggregateExpression(object=Endothelial,slot="counts",assays="RNA")
write.table(pseudocount.endo$RNA,'pseudo_counts_Endo.txt',sep='\t',quote=FALSE,col.names=TRUE,row.names=TRUE)
```

## LSEC clustering

LSECs have been described as reflecting liver zonation in their function and gene expression. Previous publications (**get citations**) have assigned (**three**)? zones:

- Zone 1 = periportal LSECs
- Zone 3 


[Su et al.](https://www.sciencedirect.com/science/article/pii/S2352345X2030206X) **mouse study**:
Arterial-like ECs = Vwf 
Midzonal markers = Lyve1 and Ctsl
Periportal landmarks: CD36, Dll4 and Efnb2... Msr1, Ltbp4, Ntn4, and Adam23
Pericentral: Kit
Central venous: Thbd
Lymphatic ECs: IL7

[MacParland et al. 2018](https://www.nature.com/articles/s41467-018-06318-7). (Note citation of bulk RNA-seq studies on sorted LSECs (Shahani et al. 2014),

**Periportal ("most abundant")**: Recently, immunofluorescent staining was used to describe the zonation of human LSECs and a population of CD36hi CD32B− CD14− LYVE1− LSECs in Zone 1 of the hepatic acinus (the periportal area). Enriched expression of F8, PECAM1, with little expression of CD32B, LYVE-1, STAB2, and CD14. (Top DE genes: MGP, SPARCL1, TM4SF1, CLEC14A, ID1, IGFBP7, ADIRF, CTGF, VWF, CD9, C7, SRPX, ID3, CAV1, GNG11, AQP1, HSPG2, EMP1, SOX18, CLDN5). In line with previous work32, we propose that these endothelial cells are likely periportal LSECs (Zone 1). 

**Central venous ("second most abundant")**: enriched expression of CD32B, LYVE1, STAB2, with little expression of VWF. LYVE-1+, CD32Bhi, CD14+, CD54+, CD36mid-lo. (Top DE genes: CCL14, CLEC1B, FCN2, S100A13, FCN3, CRHBP, STAB1, GNG11, IFI27, CLEC4G, CLDN5, CCL23, OIT3, RAMP3, SGK1, DNASE1L3, LIFR, SPARC, ADGRL4, EGFL7, PCAT19, CDKN1C). (Note Strauss et al. histological examination). 

**Non-LSECs**: The least abundant hepatic endothelial cell population (Cluster 13) was characterized by low or no expression of LSEC markers (LYVE1, STAB2, CD32B). These cells are likely non-LSEC endothelial cells including central vein and portal arterial and venous endothelial cells based on the expression of ENG (protein alias CD105) and PECAM1 (protein alias CD31) as has been described in the human liver via immunohistochemistry32. (Top DE genes: RAMP3, INMT, DNASE1L3, LIFR, PTGDS, C7, CTGF, TIMP3, RNASE1, ID3, ENG, MGP, PCAT19, HSPG2, GPM6A, PTPRB, VWF, FAM167B, SRPX, LTC4S, IFI27)


Three endothelial cell populations which were less proliferative than immune cells and expressed CALCRL (Fig. 5c.i) and RAMP2, suggesting sensitivity to adrenomedullin signaling37. 

```R
LSECs=subset(Endothelial, idents = "LSECs")
#SCTransform can be used in place of the NormalizeData, FindVariableFeatures, ScaleData workflow.
LSECs <- SCTransform(LSECs, conserve.memory=TRUE,return.only.var.genes=FALSE)
all.genes <- rownames(LSECs)
LSECs <- RunPCA(LSECs, features = all.genes)
LSECs <- RunUMAP(LSECs, dims = 1:20)
LSECs <- FindNeighbors(object = LSECs, dims = 1:20, verbose = FALSE)
LSECs <- FindClusters(object = LSECs, verbose = FALSE)

#?PercentageFeatureSet
for(x in unique(endo.markers[,2])){
    name=gsub(' ','.',x)
    features=unique(as.character(subset(endo.markers,endo.markers[,2]==x)[,1]))
    LSECs[[name]]<-PercentageFeatureSet(LSECs,features = features, assay = 'RNA')
}

#LSEC.central.venous, LSEC.periportal
FeaturePlot(object = LSECs, features = c("LSEC.central.venous.extended",
            "LSEC.periportal.extended","LSEC.mid"), label=TRUE,
            max.cutoff = c(5,0.5,1))
```

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/LSEC_UMAP.png" height="400">


### Functional pseudotime analysis

Next, we will carry out functional pseudotime analysis, following this Broad Institute [tutorial](https://broadinstitute.github.io/2019_scWorkshop/functional-pseudotime-analysis.html).

```R
LSECs.sc=as.SingleCellExperiment(LSECs)
# Run PCA. Use the runPCA function from the SingleCellExperiment package.
LSECs.sc <- runPCA(LSECs.sc, ncomponents = 50)

# Use the reducedDim function to access the PCA and store the results. 
pca <- reducedDim(LSECs.sc, "PCA")
LSECs.sc$PC1 <- pca[, 1]
LSECs.sc$PC2 <- pca[, 2]

ggplot(as.data.frame(colData(LSECs.sc)), aes(x = PC1, y = PC2, color = seurat_clusters)) + geom_quasirandom(groupOnX = FALSE) +
    theme_classic() + #scale_color_tableau() + 
    xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")

LSECs.sc$pseudotime_PC1 <- rank(LSECs.sc$PC1)  # rank cells by their PC1 score
ggplot(as.data.frame(colData(LSECs.sc)), aes(x = pseudotime_PC1, y = seurat_clusters, 
    colour = seurat_clusters)) + geom_quasirandom(groupOnX = FALSE) +
    theme_classic() + xlab("PC1") + ylab("Timepoint") +
    ggtitle("Cells ordered by first principal component")
```

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/LSEC_PCA.png" height="400">

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/LSEC_PCAtime.png" height="400">

We can see the expression of marker genes according to principal component 1:

```R
#Central venous markers
plotExpression(LSECs.sc, c("CD14","LYVE1","ICAM1","STAB2","FCGR2B","THBD"), x = "PC1", colour_by = "seurat_clusters", show_violin = FALSE,show_smooth = TRUE)
#Periportal markers
plotExpression(LSECs.sc, c("CD36","PECAM1","F8","DLL4","LTBP4","NTN4"), x = "PC1", colour_by = "seurat_clusters", show_violin = FALSE,show_smooth = TRUE)
```

The data can also be explored using a diffusion map approach:
```R
#  Prepare a counts matrix with labeled rows and columns. 
LSECs.counts <- as.data.frame(logcounts(LSECs.sc))  # access log-transformed counts matrix
cellLabels <- LSECs.sc$seurat_clusters
colnames(LSECs.counts) <- cellLabels

# Make a diffusion map.
dm <- DiffusionMap(t(LSECs.counts),n_pcs = 50)
# Plot diffusion component 1 vs diffusion component 2 (DC1 vs DC2). 
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  Timepoint = LSECs.sc$seurat_clusters)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
    geom_point() + #scale_color_tableau() + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic()
```

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/LSEC_diffusion.png" height="400">

To rename clusters in the Single Cell Experiment object:

```R
LSECs.sc@colData$seurat_clusters=gsub('1','Periportal',LSECs.sc@colData$seurat_clusters)
```

Next, see the [bigwig visualisation tutorial](https://github.com/CebolaLab/scRNA/tree/main/5.%20Bigwig%20visualisation).



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
    Endothelial2[[name]]<-PercentageFeatureSet(Endothelial2,features = features, assay = 'RNA')
}

FeaturePlot(object = Endothelial2, features = c("LSEC","nonLSECs","LSEC.central.venous.extended",
                                               "LSEC.periportal.extended"),
            label=FALSE, max.cutoff = c(20,0.6,4,1)) 
```


<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/Endo_colour_UMAP.png" height="600">

Next, see the [bigwig visualisation tutorial](https://github.com/CebolaLab/scRNA/tree/main/5.%20Bigwig%20visualisation).


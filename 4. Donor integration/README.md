
## 8. Donor integration

Up to now, biological replicates have been processed independently. This provides an opportunity to assess replicate similarity including by calculating Pearson correlation for cell-type clusters across replicates.

Following this [Seurat tutorial](https://satijalab.org/seurat/articles/integration_introduction.html), we will [integrate biological replicates which have been processed and normalised using SCTranform](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) with `SelectIntegrationFeatures` > `PrepSCTIntegration()` > `FindIntegrationAnchors()` > `IntegrateData()` (use the normalization.method = "SCT" option where appropriate). See also the [Fast integration using reciprocal PCA (RPCA) vignette](https://satijalab.org/seurat/articles/integration_rpca.html). The steps will follow:

- Create a list of Seurat objects to integrate
- Integrate data using `FindIntegrationAnchors`
- Repeat dimensionality reduction (PCA/UMAP) and clustering 

```bash
#Create a list with the Seurat objects for the donors
liver.list=list(SAMN12614699.renamed,SAMN12614700.renamed,SAMN12614710.clusters,SAMN12614716.filtered,SAMN12614717.filtered)

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

Next, see the [pseudobulk count and bigwig visualisation integration tutorial](https://github.com/CebolaLab/scRNA/tree/main/9.pseudobulk_counts_bigwigs).
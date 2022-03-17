# Analysis of 10X Genomics scRNA-seq data
This Github describes the pipeline used by the Cebola Lab to analyse single-cell RNA-seq data (scRNA-seq) from 10X Genomics. 

## 1. Introduction

This pipeline has been developed by carefully reviewing the literature for current tools and pipelines used in the analysis of 10X Genomics scRNA-seq data. This introduction presents an overview of various approaches used to analysis scRNA data, which is followed by the Cebola Lab pipeline.

Resources used are shown in the References at the bottom of this page.

#### Alignment, demultiplexing and quantification

There are several alignment algorithms to choose from, including CellRanger, STARsolo, Alevin, Alevin-fry and Kallisto; these are compared in [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741). This pipeline will use **CellRanger**, which is the analysis toolbox from 10X Genomics. [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) suggest STARsolo as an avisable alternative if memory requirement is an issue.

CellRanger includes the demultiplexing of reads based on unique barcodes, which should correspond to unique cells following secondary filtering and QC. 

#### Secondary analysis

The above tools generate a 'count matrix', which contains counts of reads for each gene, per-cell. The count matrices then undergo quality checks and filtering. The aim is to remove instances of droplets or "GEMs" (Gel Beads in EMulsion) with more or less than one cell, damaged cells, or other uninformative data. 

The pipeline for secondary analysis, including references, is discussed below.

1. **Empty droplets and ambient gene expression**. A proportion of the barcodes in the count matrix will correspond to empty droplets. A simple method to identify such droplets may use a minimum threshold of counts. However, this has the disadvantage of removing potential informative clusters of cells with low gene-expression levels, such as quiescent cells. More sophisticated methods make use of *ambient gene expression*, which refers to the background expression detected in all droplets, including empty droplets, from contaminating cell-free RNA. Modelling ambient gene expression in empty droplets can also allow other droplets to be corrected for contaminating counts, thus avoiding problems in downstream analysis and clustering. Several sophisticated methods exist with which to remove empty droplets and correct the remaining data for ambient or 'background' gene expression. Tools include **SoupX** [(Young and Behjati, 2020)](https://doi.org/10.1093/gigascience/giaa151), **EmptyDrops** [(Lun et al. 2019)](https://doi.org/10.1186%2Fs13059-019-1662-y) and **DEIM** [(Alvarez et al. 2020](https://www.nature.com/articles/s41598-020-67513-5) (note DEIM was designed for single-nuclei data, although shown to work with scRNA as well). Briefly, SoupX quantifies cell-free RNA contamination from the profiles of empty droplets and outputs background-corrected counts for the remaining GEMs. ... DEIM... *tbc*

2. **Outliers**. Common metrics used to filter uninformative GEMs such as damaged cells include the number of read counts and genes, as well as the % of mitochondrial DNA expression. For example, very high mtDNA expression can indicate cells where the cytoplasmic RNA leaked from the cell membrane, while gene expressed in little to no cells (gene-level filtering), will remove ininformative data. Earlier methods set hard cut-off thresholds, however this may remove meaningful biology, such as quiesent cells with low gene counts or highly active cells with a high level of mitochondria gene expression. 

An improved approach may be to identify clusters of cells with 

If setting cut-offs It is advisable to explore your data to set 

Any threshold for filtering genes should be informed by your experimental design, including the number of cells in the dataset and the number of cells in the smallest cluster of interest [(Leucken and Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746).

Filtering typically remove barcodes based on three typical QC measures: (1) the number of read counts per barcode, (2) the number of genes detected per barcode, and (3) the proportion of mitochondrial DNA.

The thresholds will be specific to the distribution of reads in your data, with the aim of removing barcodes with outliers. For example, [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) retained cells with gene counts \>200 and \<2,500 and a mitochondrial content \<10%. [Xu et al. (2021)](https://academic.oup.com/hmg/article/30/5/370/6131713?login=false#supplementary-data) retained cells with >30000 UMIs, 200-5000 genes, and less than 50% mitochondrial expression percentage. [Hardwick et al. (2022)](https://www.nature.com/articles/s41587-022-01231-3) excluded nuclei with unique gene counts >7,500 or <200 or >4% mitochondrial gene expression (note this is single-**nuclei** RNA-seq). 

It is challenging to determine optimal QC thresholds *a priori*, particularly in samples with heterogenous cell-types; thus it may be advisable to first apply permisive threshold and then revist QC steps one or more times to optimise the distribution of QC metrics and clusters (*not* to alter or improve the results and outcome statistics, such as differential gene expression, which falls into the trap of 'cherry picking' or 'data peeking'!) [(Leucken and Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746). It is also advisable to carry out QC filtering seperately for independent samples, due to differences in sample quality [(Plasschaert et al. 2018; ](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6108322/)[Chen, Ning and Shi. 2019)](https://doi.org/10.3389/fgene.2019.00317).

4. **Doublets**. Doublets have been removed in the past by simply removing cells with particularly high gene counts (e.g. with significant deviation). More sophisticated methods typically identify doublets based on their similarity with simulated doublets. Currently available tools include: doubletCells (Lun et al., 2016), Scrublet [(Wolock et al., 2019)](https://www.sciencedirect.com/science/article/pii/S2405471220301952#bib32), cxds (Bais and Kostka, 2020), bcds (Bais and Kostka, 2020), hybrid (Bais and Kostka, 2020), Solo [(Bernstein et al. 2020)](https://doi.org/10.1016/j.cels.2020.05.010), DoubletDetection (Gayoso and Shor, 2018), DoubletFinder ([McGinnis et al., 2019a](https://www.sciencedirect.com/science/article/pii/S2405471220304592#bib44), [2019b](https://www.sciencedirect.com/science/article/pii/S2405471220304592#bib45)), and DoubletDecon (DePasquale et al., 2019). These methods were recently compared by [Xi and Li (2021)](https://www.sciencedirect.com/science/article/pii/S2405471220304592), who report **DoubletFinder** and **Solo** as the top two performing methods. Briefly, **DoubletFinder** uses a *k*-nearest neighbors (kNN) algorithm to identify doublets based on their clustering with simulated doublets in principal component space. **Solo** (included in the scvi-tools suite from the Yosef Lab at UC Berkeley) uses a semi-supervised deep learning approach and claims improvements over DoubletFinder by not assuming linear gene expression circuits (note [Xi and Li (2021)](https://www.sciencedirect.com/science/article/pii/S2405471220304592) reported DoubletFinder as the top method).

5. **Normalization (within-sample)**. 
Within-sample normalization aims to normalise counts across cells which can differ due to sequencing depth, RNA content, and efficiency of lysis and reverse transcription [(Saket Choudhary & Rahul Satija, 2022; ](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02584-9)[Grün, Kester and van Oudenaarden. 2014)](http://scholar.google.com/scholar_lookup?&title=Validation%20of%20noise%20models%20for%20single-cell%20transcriptomics&journal=Nat%20Methods&volume=11&issue=6&pages=637-40&publication_year=2014&author=Grün%2CD&author=Kester%2CL&author=van%20Oudenaarden%2CA). Some methods normalize based on the total expression detected per-cell, however this approach is more appropriate for bulk RNA-seq as it assumes that each cell started with the same number of RNA molecules (e.g the popular toolkit **Seurat** `NormalizeData` command normalizes gene expression by the total expression, multiplied by a scale factor and log-transformed). More recent methods apply downsampling or statistical models. Tools include [Linnorm][(Yip et al., 2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7019105/#B42), SCnorm [(Bacher et al., 2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7019105/#B3), scran [(Lun et al., 2016)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7019105/#B25), and more recently, Normalisr [(Wang, 2021)](https://www.nature.com/articles/s41467-021-26682-1), sctransform[(Hafemeister and Satija, 2019)](https://doi.org/10.1186%2Fs13059-019-1874-1), bayNorm[(Tang et al. 2020)](https://www.nature.com/articles/s41467-021-26682-1#ref-CR15), and [Sanity](https://www.nature.com/articles/s41467-021-26682-1#ref-CR16).

Other methods of normalization may be considered depending on the specific study design. For example, normalising for biological covariates (such as cell cycle stage) may be useful for trajectory inference. See discussion in [(Leucken and Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746).


Note...sctransform was used by the Satija lab in their 2022 [publication](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02584-9#Sec8)

"In particular, two recent studies proposed to use generalized linear models (GLMs), where cellular sequencing depth was included as a covariate, as part of scRNA-seq preprocessing workflows. Our sctransform [9] approach utilizes the Pearson residuals from negative binomial regression as input to standard dimensional reduction techniques, while GLM-PCA [10] focuses on a generalized version of principal component analysis (PCA) for data with Poisson-distributed errors. More broadly, multiple techniques aim to learn a latent state that captures biologically relevant cellular heterogeneity using either matrix factorization or neural networks [11–13], alongside a defined error model that describes the variation that is not captured by the latent space.

**Variance stabilization**. The aim here is to correct for the relationship between gene expression and variation in expression (a well-known effect which is corrected for in bulk RNA-seq pipelines, for example by DESeq2).

5. **Normalization (between-sample)**. 


6. Imputation?

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

To run `cellranger count`, make sure your files are in the `bcl2fastq` naming convention e.g. `SRR10009414_S1_L00X_R1_001.fastq.gz` (and the corresponding `I1` and `R2`). The below command should be run, where `<ID>` is the sample ID at the start of the filename (e.g. SRR10009414) and the `<PATH>` should direct to the reference directory created by the previous command.

```bash
#Run cellranger count with the sampleID and cellranger reference directory
cellranger count --id <ID> --transcriptome <PATH>

#If working with public data i.e. pre-computed clusters:
cellranger count --nosecondary --id <ID> --transcriptome <PATH>
``` 

## 4. Secondary analysis 

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

e.g. [Xu et al. (2021)](https://academic.oup.com/hmg/article/30/5/370/6131713?login=false#supplementary-data) use `DoubletFinder` and then `decontX` of `celda` to correct for "cross containment".
3. **CellRanger**: cellranger also provides it's own tools for secondary analysis, the "`cellranger reanalyze` command reruns secondary analysis performed on the feature-barcode matrix (dimensionality reduction, clustering and visualization) using different parameter settings."


[`DropletUtils`](https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html) is an R package for filtering droplet-based scRNA-seq data. (See example use in [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741)). This package can be used to read in the count matrix to R, filter empty droplets and visualise QC plots for the distribution of counts across barcodes, as well as downsample the count matrix or raw reads etc.


 **Recommendations**: *do not* run `DoubletFinder` on aggregated scRNA-seq data representing multiple distinct samples and filter your data first to remove low-quality cell clusters. The authors recommends (1) manually threshold raw gene expression matrices according to RNA nUMIs, (2) pre-process data using [standard workflow](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html), including data normalization and scaling (3) identify clusters with low RNA UMIs, high % mitochondrial reads and/or uninformative marker genes, (4) remove clusters, pre-process again, and run DoubletFinder.
# Analysis of Liver scRNA-seq (10X genomics)
This Github describes the pipeline used by the Cebola Lab to analyse single-cell RNA-seq data (scRNA-seq) from 10X Genomics. 

This scRNA-seq analysis tutorial is split into several chapters and will focus on analysis of liver samples.

## Table of contents

1. Introduction to single-cell transcriptomics and the liver
2. [Data preparation using cellranger](https://github.com/CebolaLab/scRNA/tree/main/2.%20Data%20preparation%20with%20CellRanger)
3. [Individual dataset analysis](https://github.com/CebolaLab/scRNA/tree/main/3.%20Individual%20dataset%20analysis)
4. [Donor integration, cell type identification and pseudocount data](https://github.com/CebolaLab/scRNA/tree/main/4.%20Donor%20integration)
5. [Visualising bigwigs](https://github.com/CebolaLab/scRNA/tree/main/5.%20Pseudobulk%20counts%20and%20bigwig%20visualisation)

In this first Chapter 1, "Introduction to single-cell transcriptomics and the liver", the following topics will be discussed:

1. [Pipeline overview](#1-Pipeline-overview)  
2. [Background](#2-Background): 
    - [Hepatic cell types](#Hepatic-cell-types) and marker genes
    - [Published liver scRNA-seq studies](#Published-Liver-scRNA-seq-studies)
    - [Methods for the analysis of scRNA-seq data](#Methods-for-the-analysis-of-scRNA-seq-data)
    
## 1. Pipeline overview

1. **Generate count matrix**: `CellRanger count` is used to generate count matrices with some initial filtering to remove empty droplets. 
2. **Pre-processing**: each **biological replicate is processing seperately**. the count matrix is read into R using `Seurat`.  Initial pre-processing is carried out to prepare the data for the next steps, including normalization with a negative binominal model (`SCTransform`), merging of technical replicates and initial dimensionality reduction and clustering using `RunPCA`, `RunUMAP`, `FindNeighbors` and `FindClusters`.
3. **Correct for ambient gene expression** `SoupX` is used to correct for ambient gene expression. 
4. **QC filtering**: repeat the normalization and clustering with the corrected data. Identify and remove clusters of low-quality cells. Several rounds of pre-processing, clustering and filtering may be required.
5. **DoubletFinder**: identify and remove droplets with doublets i.e. two (or more) cells using `doubletFinder`. 
6. **Final clustering**: the cleaned data is processed for a final time and clusters are labelled using known marker genes. (Supervised and/or unsupervised clustering may be carried out). 
7. **Integrate biological replicates**
8. **Pseudo-bulk gene expression**

This pipeline has been developed by carefully reviewing current tools and best practises used in the analysis of 10X Genomics scRNA-seq data, as of March 2022. This Github will first present an overview of various available tools, followed by the Cebola Lab pipeline. Resources used are shown in the [References](#references) at the bottom of this page.

## 2. Background

The background section will be split into (1) liver cell-types, (2) a summary of previous scRNA-seq studies of the liver and (3) background on the methods used to analyse scRNA-seq data.

### Hepatic cell types

Briefly, the cell types found in human liver include:

Parenchymal cells:
- **Hepatocytes**

Non-parenchymal cells (NPCs):
- **Endothelial cells (sinusoidal)**: liver sinusoidal endothelial cells, or 'LSECs'. Highly permeable and specialised.
- **Endothelial cells (vascular)**
- **Endothelial cells (arterial)**
- **Cholangiocytes**: epithelial cells which line the bile ducts
- **Kupffer cells**: liver resident macrophages
- **Stellate cells**: secrete vitamin A and collagen in healthy and unhealthy livers, respectively

Bloodborne cells:
- T-cells
- NK-like cells
- B cells 
- Plasma cells 

The cells have important functions which contribute to the hepatic niche. Importantly, all of these cell types may be further sub-divided depending on their location, or "zonation" in the liver, with distinct cell differences between periportal and central venous regions.  
  
These are summarised in Figure 10 from [MacParland et al. (2018)](https://www.nature.com/articles/s41467-018-06318-7):  

<img src="https://github.com/CebolaLab/scRNA/blob/main/Figures/MacParland_liver_schematic.png" width="800">

> Marker genes  

Marker genes to distinguish each cell type are provided with this Github. The `xxxx` spreadsheet includes detailed information such as the sources and reported level of specificity. 

### Published Liver scRNA-seq studies

The studies listed below are published scRNA-seq studies of *human* liver:

- [MacParland et al. (2018), *Nature Comms*.](https://www.nature.com/articles/s41467-018-06318-7) - the first scRNA-seq atlas of human liver (5 livers, 8,444 cells)
- [Aizarani et al. (2019), *Nature*.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6687507/) - Nine human donors, healthy livers (~10,000 cells).
- [Ramachandran et al. (2019), *Nature*.](https://www.nature.com/articles/s41586-019-1631-3) - scRNA-seq of both healthy and cirrhotic liver, including NAFLD.
- [Zhang et al. (2020), *Front. Oncol.*](https://www.frontiersin.org/articles/10.3389/fonc.2020.596318/full) - scRNA of HCC livers with cirrhosis
- [Payen et al. (2021)](https://www.jhep-reports.eu/article/S2589-5559(21)00054-9/fulltext) - >25,000 human liver cells
- [Wang et al. (2021), *Scientific Reports*](https://www.nature.com/articles/s41598-021-98806-y) - 17,810 non-parenchymal cells from six healthy human livers. Integrated with public bulk/scRNA-seq of NASH livers.

Additional related studies:
- [Rocque et al. (2021), *Front. Immunol.*](https://doi.org/10.3389/fimmu.2021.679521) - integrated published studies to create a "RNASeq Meta-Atlas"
- [Brancale and Vilarinho]() not sure?

Reviews on the topic include: 

- [Xiong et al. (2020), *Hepatology*](10.1002/hep.31149) "A single-cell perspective of the mammalian liver in health and disease"
- [Stamataki and Swadling (2020)](https://onlinelibrary.wiley.com/doi/full/10.1111/imm.13193)
- [Ramachandran et al. (2020), *Nature*](https://www.nature.com/articles/s41575-020-0304-x) "Single-cell technologies in hepatology: new insights into liver biology and disease pathogenesis"
- [Saviano et al. (2020)](https://www.sciencedirect.com/science/article/pii/S016882782030372X) "Single-cell genomics and spatial transcriptomics: Discovery of novel cell states and cellular interactions in liver physiology and disease biology"
- [He et al. (2021), *Ann Transl Med*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8576673/) "Application of single-cell RNA sequencing technology in liver diseases: a narrative review"

The published scRNA-seq studies of the liver have uncovered many novel cell types and states in both healthy and diseased human liver. For example, [MacParland et al. (2018)](https://www.nature.com/articles/s41467-018-06318-7) reported **two distinct populations of Kupffer cells (liver macrophages), one pro-inflammatory and one immunoregulatory**. A year later and with a larger sample size, [Aizarani et al. (2019)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6687507/) reported novel cellular sub-types and zonated profiles of endothelial cells, Kupffer cells and hepatocytes, and heterogeneity of EPCAM+ hepatocyte-biased and cholangiocyte cells, as well as a population of progenitor cells with "strong potential to form bipotent liver organoids". 

Briefly, studies in mice have also uncovered interesting mechanisms, including zone-Specific alterations of LSECs in cirrhotic mouse liver [(Su et al. 2021)](https://www.sciencedirect.com/science/article/pii/S2352345X2030206X).

Look at Xiong et al. 2019 (Landscape of intercellular crosstalk in healthy and NASH liver revealed by single-cell secretome gene analysis) and Terkelsen et al. (GSE145086).

### Methods for the analysis of scRNA-seq data

#### Alignment, demultiplexing and quantification

There are several alignment algorithms to choose from, including CellRanger, STARsolo, Alevin, Alevin-fry and Kallisto; these are compared in [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741). This pipeline will use the 10X Genomics toolbox, **CellRanger**. If memory requirement is an issue, [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) suggest STARsolo as an alternative.

#### Secondary analysis

CellRanger includes the demultiplexing of reads from individual droplets, or 'GEMs' (Gel Beads in EMulsion), based on unique barcodes. These should *in theory* correspond to unique single cells. However, droplets can contain more or less than one cell, or contain damaged or low quality cells. Therefore, the secondary analysis will assess the quality of the dataset and filter the data to retain only droplets with high-quality single cells. 

The pipeline for secondary analysis, including references, is discussed below.

> **Normalization (within-sample)**  

Within-sample normalization aims to normalise counts across cells which can differ due to sequencing depth, RNA content, and efficiency of lysis and reverse transcription [(Saket Choudhary & Rahul Satija, 2022; ](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02584-9)[Grün, Kester and van Oudenaarden. 2014)](http://scholar.google.com/scholar_lookup?&title=Validation%20of%20noise%20models%20for%20single-cell%20transcriptomics&journal=Nat%20Methods&volume=11&issue=6&pages=637-40&publication_year=2014&author=Grün%2CD&author=Kester%2CL&author=van%20Oudenaarden%2CA). Some methods normalize based on the total expression detected per-cell, however this approach is more appropriate for bulk RNA-seq as it assumes that each cell started with the same number of RNA molecules (e.g **Seurat** `NormalizeData` which normalizes by the total expression, multiplied by a scale factor and log-transformed). More recent methods apply downsampling or statistical models. Tools include Linnorm [(Yip et al., 2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7019105/#B42), SCnorm [(Bacher et al., 2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7019105/#B3), scran [(Lun et al., 2016)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7019105/#B25), and more recently, Normalisr [(Wang, 2021)](https://www.nature.com/articles/s41467-021-26682-1), sctransform [(Hafemeister and Satija, 2019)](https://doi.org/10.1186%2Fs13059-019-1874-1), bayNorm [(Tang et al. 2020)](https://www.nature.com/articles/s41467-021-26682-1#ref-CR15), and Sanity [(Breda et al. 2019, *preprint*)](https://www.biorxiv.org/content/10.1101/2019.12.28.889956v1). Other normalization may be considered depending on the specific study design. For example, normalising for biological covariates (such as cell cycle stage) may be useful for trajectory inference. See discussion in [(Leucken and Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746).

See the SCTranform [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) and [vignette](https://satijalab.org/seurat/articles/sctransform_vignette.html). 


Others?? *tbc* **Variance stabilization**. The aim here is to correct for the relationship between gene expression and variation in expression (a well-known effect which is corrected for in bulk RNA-seq pipelines, for example by DESeq2. Between-sample normalization, imputation

> **Empty droplets and ambient gene expression**  

A proportion of the barcodes in the count matrix will correspond to empty droplets. A simple method to identify such droplets may use a minimum threshold of counts. However, this risks removing cells with low levels of gene expression, such as quiescent cells. More sophisticated methods identify empty droplets through *ambient gene expression*, which refers to the background expression detected in all droplets from contaminating cell-free RNA. Modelling ambient gene expression in empty droplets can also allow other droplets to be corrected for contaminating counts, thus avoiding problems in downstream analysis and clustering. Several sophisticated methods exist with which to remove empty droplets and correct the remaining data for ambient or 'background' gene expression. These include SoupX [(Young and Behjati, 2020)](https://doi.org/10.1093/gigascience/giaa151), EmptyDrops [(Lun et al. 2019)](https://doi.org/10.1186%2Fs13059-019-1662-y) and DEIM [(Alvarez et al. 2020](https://www.nature.com/articles/s41598-020-67513-5) (note DEIM was designed for single-nuclei data, although shown to work with scRNA as well). This pipeline will use **SoupX**. 

> **Filtering low-quality and uninformative barcodes**  

GEMs containing damaged cells can be identified using commonly-used metrics: **(1) the number of read counts per barcode, (2) the number of genes detected per barcode & (3) the proportion of mitochondrial DNA.** For example, very high mtDNA expression can indicate cells where the cytoplasmic RNA leaked from the cell membrane. Hard cut-offs may be used to filter the data, however this may remove meaningful biology, such as quiesent cells with low gene counts or highly active cells with a high level of mitochondria gene expression. Prior knowledge of the expected cell types may inform the use of thresholds. 

**QC thresholds should be specific to your data**. For example, [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) retained cells with gene counts \>200 and \<2,500 and a mitochondrial content \<10%. [Xu et al. (2021)](https://academic.oup.com/hmg/article/30/5/370/6131713?login=false#supplementary-data) retained cells with >30,000 UMIs, 200-5,000 genes, and less than 50% mitochondrial expression. It is challenging to determine optimal QC thresholds *a priori*, particularly in samples with heterogenous cell-types; thus it is advisable to first apply permisive threshold and then revist QC steps one or more times to optimise the distribution of QC metrics and clusters (*not* to alter or improve the results and outcome statistics, such as differential gene expression, which falls into the trap of 'cherry picking' or 'data peeking'!) [(Leucken and Theis, 2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746). It is also advisable to carry out QC filtering seperately for independent samples, due to differences in sample quality [(Plasschaert et al. 2018; ](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6108322/)[Chen, Ning and Shi. 2019)](https://doi.org/10.3389/fgene.2019.00317).

> **Doublets**   

Doublets have been removed in the past by simply removing cells with particularly high gene counts (e.g. with significant deviation). However, this can remove active cell types with high levels of gene expression. In addition, hepatocytes and other liver cells types can be bi-nucleated, leading to single-cell droplets with increased counts. More sophisticated methods typically identify doublets based on their similarity with simulated doublets. Currently available tools include: doubletCells (Lun et al., 2016), Scrublet [(Wolock et al., 2019)](https://www.sciencedirect.com/science/article/pii/S2405471220301952#bib32), cxds (Bais and Kostka, 2020), bcds (Bais and Kostka, 2020), hybrid (Bais and Kostka, 2020), Solo [(Bernstein et al. 2020)](https://doi.org/10.1016/j.cels.2020.05.010), DoubletDetection (Gayoso and Shor, 2018), DoubletFinder ([McGinnis et al., 2019a](https://www.sciencedirect.com/science/article/pii/S2405471220304592#bib44), [2019b](https://www.sciencedirect.com/science/article/pii/S2405471220304592#bib45)), and DoubletDecon (DePasquale et al., 2019). These methods were recently compared by [Xi and Li (2021)](https://www.sciencedirect.com/science/article/pii/S2405471220304592), who report **DoubletFinder** and **Solo** as the top two performing methods. Briefly, **DoubletFinder** uses a *k*-nearest neighbors (kNN) algorithm to identify doublets based on their clustering with simulated doublets in principal component space. **Solo** (included in the scvi-tools suite from the Yosef Lab at UC Berkeley) uses a semi-supervised deep learning approach and claims improvements over DoubletFinder by not assuming linear gene expression circuits (note [Xi and Li (2021)](https://www.sciencedirect.com/science/article/pii/S2405471220304592) reported DoubletFinder as the top method). ([DEIM](https://www.nature.com/articles/s41598-020-67513-5) may be useful for single-nuclei experiments).

See more information in the [OSCA vignette](http://bioconductor.org/books/3.13/OSCA.advanced/doublet-detection.html). "This site contains the advanced analysis chapters for the “Orchestrating Single-Cell Analysis with Bioconductor” book." (2021-05)

#### Downstream analysis

> **Feature selection, dimensionality reduction and visualisation**    

Including *clustering and cell identity*
*TBC*

Batch effect correction: methods include [ComBat-seq](https://academic.oup.com/nargab/article/2/3/lqaa078/5909519)
Batch correction using [Harmony](https://portals.broadinstitute.org/harmony/articles/quickstart.html), BBKNN, and Seurat v3 with processed data and highly variable genes. Harmony better than Seurat v3 according to [https://www.nature.com/articles/s41597-021-00809-x#Sec2]

Integrate datasets which have been normalised using SCTransform [vignette](https://satijalab.org/seurat/archive/v3.1/integration.html#sctransform),

#### Dataset integration

Here, we opt to process biological replicates independently and then integrate them using the [recommended Seurat pipeline](https://satijalab.org/seurat/articles/integration_introduction.html) on integrating datasets (see also [this link](https://satijalab.org/seurat/archive/v3.1/integration.html#sctransform)). Specifically, we follow steps to integrate datasets normalized with SCTransform (`PrepSCTIntegration()` > `FindIntegrationAnchors()` > `IntegrateData()`, with the normalization.method parameter to the value SCT).

(Note, [Harmony](https://portals.broadinstitute.org/harmony/articles/quickstart.html) may be an alternative.)

#### Differential expression and pseudobulk expression  
*TBC*

Next, see the [donor integration tutorial](https://github.com/CebolaLab/scRNA/tree/main/8.donor_integration).



# References and resources

- [Leucken and Theis (2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746): A 2019 effort to compile current best practises in scRNA-seq. This paper is very useful for "newbies" and gives an excellent overview of the essential steps of scRNA-seq analysis (note that some specific tools mentioned are superceeded by more recent published tools).
- https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
- [svci-tools User Guide](https://docs.scvi-tools.org/en/stable/user_guide/index.html)
- [SoupX vignette](https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html)
- Wellcome single-cell course [using Seurat](https://www.singlecellcourse.org/single-cell-rna-seq-analysis-using-seurat.html)
- [Advanced chapters from the “Orchestrating Single-Cell Analysis with Bioconductor” book](http://bioconductor.org/books/3.13/OSCA.advanced/)

e.g. [Xu et al. (2021)](https://academic.oup.com/hmg/article/30/5/370/6131713?login=false#supplementary-data) use `DoubletFinder` and then `decontX` of `celda` to correct for "cross containment".
3. **CellRanger**: cellranger also provides it's own tools for secondary analysis, the "`cellranger reanalyze` command reruns secondary analysis performed on the feature-barcode matrix (dimensionality reduction, clustering and visualization) using different parameter settings."

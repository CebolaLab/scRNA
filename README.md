# Analysis of 10X Genomics single-cell RNA-seq
This Github describes the pipeline used by the Cebola Lab to analyse single-cell RNA-seq data (scRNA-seq) from 10X Genomics. 

## 1. Introduction

This pipeline has been developed by carefully reviewing the literature for current tools and pipelines used in the analysis of 10X Genomics scRNA-seq data. This introduction presents an overview of various approaches used to analysis scRNA data, which is followed by the Cebola Lab pipeline.

Resources used are shown in the References at the bottom of this page.

#### Alignment, demultiplexing and quantification

There are several alignment algorithms to choose from, including CellRanger, STARsolo, Alevin, Alevin-fry and Kallisto; these are compared in [Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741). This pipeline will use **CellRanger**, which is the analysis toolbox from 10X Genomics ([Brüning et al. (2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) suggest STARsolo as an avisable alternative if memory requirement is an issue).

CellRanger includes the demultiplexing of reads based on unique barcodes, which should correspond to unique cells following secondary filtering and QC. 

#### Secondary analysis

The above tools generate a 'count matrix', which contains counts of reads for each gene, per-cell. The count matrices then undergo secondary analysis, including QC and filtering. The aim is to remove instances of droplets or "GEMs" (Gel Beads in EMulsion) with more or less than one cell. Filtering typically remove barcodes based on three typical QC measures: (1) the number of read counts per barcode, (2) the number of genes detected per barcode, and (3) the proportion of mitochondrial DNA.

The thresholds will be specific to the distribution of reads in your data, with the aim of removing barcodes with outliers. For example, [(Brüning et al. 2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) retained cells with gene counts \>200 and \<2,500 and a mitochondrial content \<10%. [Xu et al. (2021)](https://academic.oup.com/hmg/article/30/5/370/6131713?login=false#supplementary-data) retained cells with >30000 UMIs, 200-5000 genes, and less than 50% mitochondrial expression percentage. [Hardwick et al. (2022)](https://www.nature.com/articles/s41587-022-01231-3) excluded nuclei with unique gene counts >7,500 or <200 or >4% mitochondrial gene expression (note this is single-**nuclei** RNA-seq). 

Some examples of tools available for secondary analysis are discussed below.

- **Remove doublets** e.g. [Xu et al. (2021)](https://academic.oup.com/hmg/article/30/5/370/6131713?login=false#supplementary-data) use `DoubletFinder` and then `decontX` of `celda` to correct for "cross containment".
-  [`DropletUtils`](https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html) is an R package for filtering droplet-based scRNA-seq data. (See example use in [(Brüning et al. 2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741)). Functions include "load in the matrix of unique molecule identifier (UMI) counts as well as the raw molecule information... downsampling the UMI count matrix or the raw reads; for distinguishing cells from empty droplets, based on the UMI counts; and to eliminate the effects of barcode swapping on Illumina 4000 sequencing machine."
- **CellRanger**: cellranger also provides it's own tools for secondary analysis, the "`cellranger reanalyze` command reruns secondary analysis performed on the feature-barcode matrix (dimensionality reduction, clustering and visualization) using different parameter settings."

#### Dataset integration

Integrating multiple scRNA-seq datasets presents an additional challenge, which may again be tackled with different methods. Note, the below are specifically for integrating multiple sequencing runs of different GEM Wells. For samples sequenced in the same GEM well, pass the multiple fastq files to `cellranger count` using the `--fastqs` argument. 

- [(Brüning et al. 2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) integrate expression matrices using Suerat, including: normalization with the `SCTransform` function, rank the features using the `SelectIntegrationFeatures` function, with the resulting features controlled using the function `PrepSCTIntegration`. Anchors were determined by `FindIntegrationAnchors` and afterwards used with the `IntegrateData` function. 
- With **CellRanger** using `cellranger aggr` which aggregates multiple runs of `cellranger count`, normalizes runs to the same effective sequencing depth (the pipeline equalizes the average read depth per cell between groups before merging), and then performs secondary analysis on the combined data (this produces one single-feature barcode matrix and a .cloupe file for visualizing with Loupe Browser). `cellranger aggr` uses Chemistry Batch Correction when aggregating resuts from a combination of 5' and 3', or 3' v2 and 3' v3 Gene Expression data, which improves the mixing of the batches in the t-SNE visualization and clustering results (note that residual batch effects may still be present). 

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



- [(Brüning et al. 2022)](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741) filter cells using the R packages DropletUtils and then use `Seurat` for downstream analysis, retaining cells with gene counts \>200 and \<2,500 and a mitochondrial content \<10%.
- [Xu et al. (2021)](https://academic.oup.com/hmg/article/30/5/370/6131713?login=false#supplementary-data) "removed doublet cells using “DoubletFinder”. Then we use the function “decontX” of “celda” to correct the probable cross containment. Then the corrected expression matrix was processed by “Seurat V3.14”. For the quality control, cells with 0~30000 UMIs, 200~5000 genes, and less than 50% mitochondrial expression percentage were filtered out for the next analysis."

For single-nuclei RNA, [Hardwick et al. (2022)](https://www.nature.com/articles/s41587-022-01231-3) excluded nuclei with unique gene counts >7,500 or <200 or >4% mitochondrial gene expression. UMI numbers and mitochondrial gene expression % were regressed from each nucleus and the matrix was log-normalised and scaled to 10,000 reads per cell. Performed both tSNE and UMAP non-lnear reduction techniques... cell types assigned by canonical marker genes for each cluster... cell type annotation confirmed by aligning to <other data>. (Hardwick et al. 2022, Nature Biotechnology).

# References

- [Leucken and Theis (2019)](https://www.embopress.org/doi/full/10.15252/msb.20188746): A 2019 effort to compile current best practises in scRNA-seq. This paper is very useful for "newbies" and gives an excellent overview of the essential steps of scRNA-seq analysis (note that some specific tools mentioned are superceeded by more recent published tools).
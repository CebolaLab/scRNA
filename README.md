### Create a conda environment

This pipeline will process data from Chromium Single Cell A Chip Kit (10X Genomics) using Seurat v4 and Cell Ranger v6.1.

A conda environment for the analysis can be created as follows:
```bash
conda create -n scRNA -c conda-forge r-base r-essentials
conda install -n scRNA -c r r-irkernel
conda install -n scRNA -c bioconda r-seurat
conda install -n scRNA -c hcc cellranger
```
The first steps will follow the recommended pipeline from [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ov).
First, trim align fastq files to the reference genome (GRCh38).


```R

```

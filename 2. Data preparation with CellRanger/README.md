# Data preparation with CellRanger

This pipeline will process data from Chromium Single Cell A Chip Kit (10X Genomics) using Seurat v4 and Cell Ranger v6.1.
First, the raw sequence data from 10X Genomics should be processed using Cell Ranger to generate raw and filtered counts.

You should install the latest version of Cell Ranger, by registering at [this link](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest). 

Create a folder to save the new CellRanger code, `cd` and run the `wget` or `curl` command provided following the registration. Then:

```bash
tar -zxvf cellranger-6.1.2.tar.gz
```

Next add the path with the CellRanger executable to your PATH. **NOTE: you will have to run this command every time you want to use CellRanger, or include it in your CellRanger scripts**.

```bash
export PATH=/rds/general/user/hm1412/home/anaconda3/envs/scRNA2/bin/cellranger-6.1.2:$PATH
```

The next steps will be to:

1. [Generate a reference transcriptome](#2-Generate-a-reference-transcriptome)
2. [Cellranger count](#5-Cellranger-count)

## 1. Generate a reference transcriptome

The first steps will follow the recommended pipeline from [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ov).

When run on the raw sequence data, `cellranger count` will generate a count matrix reflecting the number of reads for each gene in the reference transcriptome. This will require you to either download the prebuilt reference, or generate you own. Here, here we will compile our own transcriptome using our preferred versions of the reference genome and gene annotation:

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

## 2. Cellranger count

To run `cellranger count`, make sure your files are in the `bcl2fastq` naming convention e.g. `SRR10009414_S1_L00X_R1_001.fastq.gz` (and the corresponding `I1` and `R2`). The below command should be run, where `<ID>` is the sample ID at the start of the filename (e.g. SRR10009414) and the `<PATH>` should direct to the reference directory created by the previous command. **Technical replicates should be combined here**.

```bash
#Run cellranger count with the sampleID and cellranger reference directory
cellranger count --id <ID> --transcriptome <PATH>

#If working with public data i.e. pre-computed clusters:
cellranger count --nosecondary --id <ID> --transcriptome <PATH>

#Example for donor SAMN12614700, with four sets of fastq files from four runs in the SAMN12614700 directory:
#cellranger count --nosecondary --id SAMN12614700 --sample SRR10009414,SRR10009415,SRR10009416,SRR10009417 --transcriptome $GENOMEDIR/GRCh38 --fastqs SAMN12614700/
``` 

Next, check out the (individual dataset analysis pipeline)[https://github.com/CebolaLab/scRNA/tree/main/3.%20Individual%20dataset%20analysis] tutorial.
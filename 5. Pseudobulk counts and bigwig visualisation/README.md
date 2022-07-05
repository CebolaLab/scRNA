
## 9. Pseudobulk RNA counts

```bash
#Save the cluster identity for each cell
identity=as.data.frame(Idents(LSECs.no9))
write.table(identity,'cell_identity.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=FALSE)

#Generate the pseudobulk data and save
pseudocount.LSEC=AggregateExpression(object=LSECs,slot="counts",assays="RNA")
write.table(pseudocount.LSEC$RNA,'pseudo_counts_LSEC.txt',sep='\t',quote=FALSE,col.names=TRUE,row.names=TRUE)
```

## 10. Generate bigwigs and visualise data

Here, we will generate bigwig files for each cell type. This will use `sinto` to subset bam files for each donor according to the cell type ID. Then, corresponding bam files will be merged and bigwigs generated using `deeptools` and `macs2`.

1. Split bam files according to cluster
2. Merge cluster-specific bam files across donors
3. Generate coverge bigwigs
4. Generate signal bigwigs

**Split bam files by cluster**

First, you will need the `cell_identity.txt` file created above. This will contain the cell barcode, appended with the sample number (added during the merging step), and the assigned cluster (cell type). It will will look something like:

```bash
AAACCTGCAAAGGTGC-1_1	Kupffer cell
AAACCTGGTCTCGTTC-1_1	monocyte
AAAGATGAGGGATACC-1_1	NK cell
AAAGTAGCAAGGACTG-1_1	Kupffer cell
AAATGCCCAGCATGAG-1_1	monocyte
```

You will need to know the order in which the donors were merged, e.g:

_1 = donor 1 \\
_2 = donor 2 \\
_3 = donor 3 \\

Here, we assume that each donor has a diretory, such that the following two files are visible:

```bash
donor_1/outs/possorted_genome_bam.bam
donor_1/outs/possorted_genome_bam.bam.bai
donor_2/outs/possorted_genome_bam.bam
donor_2/outs/possorted_genome_bam.bam.bai
donor_3/outs/possorted_genome_bam.bam
donor_3/outs/possorted_genome_bam.bam.bai
```

Using `sinto filterbarcodes`, the bam files for each donor will be split according to the cluster. An *array job* can be used to split the position sorted bam file for each donor:

```bash
tid=$PBS_ARRAY_INDEX
declare -a donors=(donor1 donor2 donor3)
declare -a id=(1 2 3)

#Setting an array to run three jobs (i.e. 0-2) will loop through the arrays to analyse donor1, 2 and 3 independently. 
#This next step will extract the corresponding cell IDs for the donor, e.g. lines which have _1, _2 or _3 depending on which job in the array is running.
grep "_${id[$tid]}" cell_identity.txt > cellIDs
#cellIDs is a file which will only have the cellIDs with _1. We now need to remove _1 as this was added in the merging step and is missing from the cell IDs in the original bam file.
sed -i "s/_${id[$tid]}//g" cellIDs 
#Make sure to use the double quotes as this will allow for a variable within the command.


sinto filterbarcodes -c cellIDs -b ${donors[$tid]}/outs/possorted_genome_bam.bam -p 6 --outdir ${donors[$tid]}
#-p is the number of processors
```

Next, we will use sinto 


**Merging bam files**
```bash
#Create an array with the names of the clusters
#Option 1 - use a file, here cell_types.txt with each cluster name on a new line
cells=$(cat $DIR/make_bigwigs/cell_types.txt)
arr=($cells)
#Option 2 - define an array and specify each cluster name
declare -a array cells=(7 B Central Cholangiocyte Kupffer macrophage monocyte NK otherEC Pericentral)

#With an array job, each job will loop through the list 
#e.g. on a PBS system, tid=$PBS_ARRAY_INDEX
cell=${arr[$tid]}

#Use samtools merge to merge the bam files. The output file will have the cluster/cell type name and will merge the corresponding bam files in the donor subdirectories. 
samtools merge -o "$cell".bam $DIR/*/"$cell".bam
samtools index "$cell".bam

#The coverage bigwig is generated using deeptools bamCoverage
#The effective genome size will differ depending on the read length and reference genome used
bamCoverage -b "$cell".bam -o "$cell".bw --normalizeUsing BPM --effectiveGenomeSize 2913022398
```


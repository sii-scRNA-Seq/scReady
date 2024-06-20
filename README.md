# scRNAseq-Standardised-pipeline
A Standardised basic pipeline, to quickly go from Cellranger mapping to the Seurat object

From fastqfiles it generates SeuratObject, QC and save it as RDS (quickly!)

At moment you can run the script on MARS with

> sbatch cellranger_seurat.sh LIBINFO.csv

LIBINFO.csv should look like:

> fastqs,sample,count introns
> 
> /path/to/fastq/files,FASTQID1,true
>
> /path/to/fastq/files,FASTQID2,true

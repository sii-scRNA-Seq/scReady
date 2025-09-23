# scRNAseq-Standardised-pipeline
A Standardised basic pipeline, to quickly go from Cellranger mapping to the Seurat object

From fastqfiles it generates SeuratObject, QC and save it as RDS (quickly!)

At moment you can run the script on MARS with

> sbatch cellranger_seurat.sh lib.csv <optional_metadata_file.csv>

lib.csv should look like:

> fastqs,sample,count introns,generate bam,output folder
> 
> /path/to/fastq/files/,FASTQID1,true,false,/path/to/output/folder1
> 
> /path/to/fastq/files/,FASTQID2,false,true,/path/to/output/folder1
> 
> /path/to/fastq/files/,FASTQID3,true,true,/path/to/output/AnotherFolder

You can also run just the Seurat part with:

> sbatch seurat_run_only.sh /path/to/output/

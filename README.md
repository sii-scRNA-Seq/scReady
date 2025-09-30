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


## Quick start
Requirements:
* Docker for local runs, Apptainer/Singularity on HPC
* An input directory `/path/to/output/` containing Cellranger outputs

### Running locally

1) Pull Docker image `docker pull ghcr.io/sii-scrna-seq/scrnaseq-standardised-pipeline:latest`. If you want a specific version, substitute this for `latest`.
2) Run Seurat part of pipeline:
```
docker run --rm -v "$PWD":/work \
    ghcr.io/sii-scrna-seq/scrnaseq-standardised-pipeline:latest \
    /path/to/output/
```

### Running on MARS/HPC
Run `sbatch run_slurm.sh /path/to/output/`

## Contributing
If you need to update R dependencies, edit `deps.R` and regenerate the `renv.lock`:
1) (optional) Run `docker build -f Dockerfile.locksmith -t renv-locksmith .` if the `renv-locksmith` container hasn't been built before.
2) Run `docker run --rm -v "$PWD":/work renv-locksmith` to save updated `renv.lock` - make sure to commit this alongside any code changes!

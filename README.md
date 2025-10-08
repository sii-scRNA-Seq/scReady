# scRready
A pipeline, that automates and standardizes preprocessing for scRNA-seq data

## Quick start
Requirements:
* Docker for local runs, Apptainer/Singularity on HPC
* An input directory `/path/to/output/` containing Cellranger outputs

### Running locally
1) Pull Docker image `docker pull ghcr.io/sii-scrna-seq/scready:latest`. If you want a specific version, substitute this for `latest`.
2) (optional) Save a custom config in current directory which will override the default config. Either:

```
docker run ghcr.io/sii-scrna-seq/scready --print-default-config > scReady.config
``` 
Or:

```
docker run ghcr.io/sii-scrna-seq/scready init-config
```
Then edit `scReady.config` with desired parameters. 

3) Run the pipeline:
```
docker run --rm -v "$PWD":/work \
    ghcr.io/sii-scrna-seq/scready \
    /path/to/output/
```

### Running on MARS/HPC
Run `sbatch run_slurm.sh /path/to/output/`

## Contributing
If you need to update R dependencies, edit `deps.R` and regenerate the `renv.lock`:
1) (optional) Run `docker build -f Dockerfile.locksmith -t renv-locksmith .` if the `renv-locksmith` container hasn't been built before.
2) Run `docker run --rm -v "$PWD":/work renv-locksmith` to save updated `renv.lock` - make sure to commit this alongside any code changes!

[![DOI](https://zenodo.org/badge/817825035.svg)](https://doi.org/10.5281/zenodo.17789689)
[![DeepWiki](https://img.shields.io/badge/DeepWiki-sii--scRNA--Seq%2FscReady-blue.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACwAAAAyCAYAAAAnWDnqAAAAAXNSR0IArs4c6QAAA05JREFUaEPtmUtyEzEQhtWTQyQLHNak2AB7ZnyXZMEjXMGeK/AIi+QuHrMnbChYY7MIh8g01fJoopFb0uhhEqqcbWTp06/uv1saEDv4O3n3dV60RfP947Mm9/SQc0ICFQgzfc4CYZoTPAswgSJCCUJUnAAoRHOAUOcATwbmVLWdGoH//PB8mnKqScAhsD0kYP3j/Yt5LPQe2KvcXmGvRHcDnpxfL2zOYJ1mFwrryWTz0advv1Ut4CJgf5uhDuDj5eUcAUoahrdY/56ebRWeraTjMt/00Sh3UDtjgHtQNHwcRGOC98BJEAEymycmYcWwOprTgcB6VZ5JK5TAJ+fXGLBm3FDAmn6oPPjR4rKCAoJCal2eAiQp2x0vxTPB3ALO2CRkwmDy5WohzBDwSEFKRwPbknEggCPB/imwrycgxX2NzoMCHhPkDwqYMr9tRcP5qNrMZHkVnOjRMWwLCcr8ohBVb1OMjxLwGCvjTikrsBOiA6fNyCrm8V1rP93iVPpwaE+gO0SsWmPiXB+jikdf6SizrT5qKasx5j8ABbHpFTx+vFXp9EnYQmLx02h1QTTrl6eDqxLnGjporxl3NL3agEvXdT0WmEost648sQOYAeJS9Q7bfUVoMGnjo4AZdUMQku50McDcMWcBPvr0SzbTAFDfvJqwLzgxwATnCgnp4wDl6Aa+Ax283gghmj+vj7feE2KBBRMW3FzOpLOADl0Isb5587h/U4gGvkt5v60Z1VLG8BhYjbzRwyQZemwAd6cCR5/XFWLYZRIMpX39AR0tjaGGiGzLVyhse5C9RKC6ai42ppWPKiBagOvaYk8lO7DajerabOZP46Lby5wKjw1HCRx7p9sVMOWGzb/vA1hwiWc6jm3MvQDTogQkiqIhJV0nBQBTU+3okKCFDy9WwferkHjtxib7t3xIUQtHxnIwtx4mpg26/HfwVNVDb4oI9RHmx5WGelRVlrtiw43zboCLaxv46AZeB3IlTkwouebTr1y2NjSpHz68WNFjHvupy3q8TFn3Hos2IAk4Ju5dCo8B3wP7VPr/FGaKiG+T+v+TQqIrOqMTL1VdWV1DdmcbO8KXBz6esmYWYKPwDL5b5FA1a0hwapHiom0r/cKaoqr+27/XcrS5UwSMbQAAAABJRU5ErkJggg==)](https://deepwiki.com/sii-scRNA-Seq/scReady)

# scRready
A pipeline, that automates and standardizes preprocessing for scRNA-seq data

## Ask DeepWiki

You can ask [DeepWiki](https://deepwiki.com/sii-scRNA-Seq/scReady) about scReady documentation

## Quick start
Requirements:
- Docker for local runs, Apptainer/Singularity on HPC
- An input directory (`/path/to/output/`) containing **one subfolder per sample** (most mapping tools will prepare this by default), where each subfolder includes:
  - Cell Ranger outputs (e.g., `outs/filtered_feature_bc_matrix/`) or `.h5` files.
  - Example structure (all these options are accepted):
    ```
    /path/to/output/
    ├── sample1/
    │   └── outs/
    │       ├── filtered_feature_bc_matrix/
    │       ├── raw_feature_bc_matrix/
    │       ├── analysis
    │       └── ...
    ├── sample2/
    │   ├── filtered_feature_bc_matrix/
    │   ├── raw_feature_bc_matrix/
    │   ├── analysis
    │   └── ...
    ├── sample3/
    │   ├── filtered_feature_bc_matrix.h5
    │   ├── raw_feature_bc_matrix.h5
    │   └── ...
    └── ...
    ```
- A metadata file (`/path/to/metadata.csv`) with sample annotations (e.g., `sample_id`, `protocol`). Example:
  ```csv
  sample_id,protocol,condition
  sample1,Protocol1,Control
  sample2,Protocol1,Treated
  sample3,Protocol2,Control

### Running locally

#### Docker

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
docker run -v "$PWD":/work ghcr.io/sii-scrna-seq/scready /path/to/output/ /path/to/metadata.csv
```

don't forget to add nohup/& to run it in background
```
nohup docker run -v "$PWD":/work ghcr.io/sii-scrna-seq/scready /path/to/output/ /path/to/metadata.csv &
```


#### Apptainer

1) Pull scReady image (if you want a specific version, substitute this for `latest`).
```
apptainer build scready_latest.sif docker://ghcr.io/sii-scrna-seq/scready:latest
```

2) (optional) Save a custom config in current directory which will override the default config:
```
apptainer run docker://ghcr.io/sii-scrna-seq/scready --print-default-config > scReady.config
```
Then edit `scReady.config` with desired parameters.

Run scReady (in background):
```
nohup apptainer run scready_latest.sif /path/to/output/ /path/to/metadata.csv > scready.log 2>&1 &
```

You might need to adjust the --bind flag if your data is not in the current directory.

### Running on HPC using slurm
1) Run the following to collect the latest `run_slurm.sh` script
```
wget -qO run_slurm.sh https://raw.githubusercontent.com/sii-scrna-seq/scready/main/dist/run_slurm.sh && chmod +x run_slurm.sh
wget -qO scReady.config https://github.com/sii-scRNA-Seq/scReady/raw/refs/heads/main/config/scReady.config
```
Modify it a accordingly to your HPC requests (ask your admin/facility)

2) Run `sbatch run_slurm.sh /path/to/output/ /path/to/metadata.csv`

## Contributing
If you need to update R dependencies, edit `deps.R` and regenerate the `renv.lock`:
1) (optional) Run `docker build -f Dockerfile.locksmith -t renv-locksmith .` if the `renv-locksmith` container hasn't been built before.
2) Run `docker run --rm -v "$PWD":/work renv-locksmith` to save updated `renv.lock` - make sure to commit this alongside any code changes!

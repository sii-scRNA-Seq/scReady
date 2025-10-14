# scRready
A pipeline, that automates and standardizes preprocessing for scRNA-seq data

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
Or
```
wget https://github.com/sii-scRNA-Seq/scReady/raw/refs/heads/main/config/scReady.config
```
Then edit `scReady.config` with desired parameters.

Run scReady (in background):
```
nohup apptainer run scready_latest.sif /path/to/output/ /path/to/metadata.csv > scready.log 2>&1 &
```

You might need to adjust the --bind flag if your data is not in the current directory.

### Running on MARS/HPC
1) Run the following to collect the latest `run_slurm.sh` script
```
wget -qO run_slurm.sh https://raw.githubusercontent.com/sii-scrna-seq/scready/main/dist/run_slurm.sh && chmod +x run_slurm.sh
```
Modify it a accordingly to your HPC requests (ask your admin/facility)

2) Run `sbatch run_slurm.sh /path/to/output/ /path/to/metadata.csv`

## Contributing
If you need to update R dependencies, edit `deps.R` and regenerate the `renv.lock`:
1) (optional) Run `docker build -f Dockerfile.locksmith -t renv-locksmith .` if the `renv-locksmith` container hasn't been built before.
2) Run `docker run --rm -v "$PWD":/work renv-locksmith` to save updated `renv.lock` - make sure to commit this alongside any code changes!

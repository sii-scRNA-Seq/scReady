FROM bioconductor/bioconductor_docker:RELEASE_3_20

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libglpk-dev \
        libcairo2-dev libxt-dev xorg-dev libpng-dev libjpeg-dev libtiff5-dev \
        libmagick++-dev \
    && rm -rf /var/lib/apt/lists/*

# App dir; separate from the working dir you’ll mount
WORKDIR /opt/app

# Install renv and restore from lockfile
COPY renv.lock /opt/app/renv.lock

RUN R -q -e "install.packages(c('renv','BiocManager'), repos='https://cloud.r-project.org')"

ENV RENV_DOWNLOAD_FILE_METHOD=curl \
    R_DEFAULT_INTERNET_TIMEOUT=1800

RUN R -q -e "options(repos = BiocManager::repositories(), \
                    download.file.method = 'libcurl', \
                    timeout = 1800); \
             renv::restore(lockfile='/opt/app/renv.lock', prompt=FALSE)"

COPY R /opt/app/R
COPY run.sh /usr/local/bin/run.sh
RUN chmod +x /usr/local/bin/run.sh

# Working directory for user data. You will bind-mount your project here.
WORKDIR /work

# Nice-to-haves for deterministic behavior
# ENV RENV_CONFIG_SANDBOX_ENABLED=false \
#     RENV_CONFIG_CACHE_SYMLINKS=true \
#     R_DEFAULT_INTERNET_TIMEOUT=300

ENTRYPOINT ["/usr/local/bin/run.sh"]

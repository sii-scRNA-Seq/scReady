install.packages("BiocManager", repos="https://cloud.r-project.org")
options(repos = BiocManager::repositories())

renv::install("bioc::BiocVersion")

# Add runtime packages here (pin if you like with @version)
renv::install(c(
  "Seurat",
  "SoupX",
  "bioc::scDblFinder",
  # "bioc::Rhtslib" tarball is corrupted, using GitHub source instead.
  "git::https://git.bioconductor.org/packages/Rhtslib@RELEASE_3_20",
  "BimberLab/cellhashR",
  "harmony",
  "hdf5r"
))

#============================================================================================
#Create user-defined directory to install R packages
#Install single-cell RNA-seq analysis packages
#============================================================================================

# Create a user-writable library folder if it doesn't exist
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)

# Tell R to install packages into this user folder first
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# Install SeuratData and SeuratWrappers from GitHub
BiocManager::install("satijalab/seurat-data", update = FALSE)
BiocManager::install("satijalab/seurat-wrappers", update = FALSE)

# Install Azimuth (specify ref = "master" as required by the repo)
remotes::install_github("satijalab/azimuth", ref = "master")




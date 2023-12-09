# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c(
    "basilisk",
    "dplyr",
    "Matrix",
    "MatrixGenerics",
    "progress",
    "Seurat",
    "stringr",
    "withr"
  )
)

options(clustermq.scheduler = "multicore")

tar_source()
# source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_download(
    midgut.counts,
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120537/suppl/GSE120537%5Fcounts.csv.gz',
    'GSE120537_counts.csv.gz'
  ),
  tar_download(
    midgut.metadata,
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120537/suppl/GSE120537%5Fmetadata.csv.gz',
    'GSE120537_metadata.csv.gz'
  ),
  tar_download(
    midgut.gtf,
    'https://ftp.flybase.org/releases/FB2019_02/dmel_r6.27/gtf/dmel-all-r6.27.gtf.gz',
    'dmel-all-r6.27.gtf.gz'
  ),
  tar_target(
    OptimalSPCA,
    git_clone_OptimalSPCA(),
    format='file',
    # Git clone and checkout a specific commit hash - never needs to be
    # refreshed
    cue=tar_cue('never')
  ),
  tar_target(
    OptimalSPCADepot,
    julia_pkg_OptimalSPCADepot(),
    format='file',
    cue=tar_cue('never')
  ),
  tar_target(
    indrop,
    midgut_seurat_for_technology(midgut.counts, midgut.metadata, midgut.gtf, 'inDrop')
  ),
  tar_target(
    indrop.2,
    RunSparsePCA(indrop, "covar", varnum=8, npcs=60, eigen_gap=0.05, search_cap=500000)
  )
)

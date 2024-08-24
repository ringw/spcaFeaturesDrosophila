# Sets up the spcaFeatures function(s). User will need to setwd to the
# spcaFeaturesDrosophila directory when sourcing this script, and when running
# the functions including RunSparsePCA().
library(basilisk)
library(dplyr)
library(MatrixGenerics)
library(progress)
library(Seurat)
library(withr)

source("_targets.R")

# Does not affect the current R process Sets up the "OptimalSPCA" directory
# instead.
targets::tar_make(OptimalSPCA | OptimalSPCADepot)

# Packages are now loaded and installation created in the
# "spcaFeaturesDrosophila" working directory. The current R process is ready to
# call RunSparsePCA().

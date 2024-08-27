## Sparse PCA For Interpreting Single-Cell Clusters And Transcriptional Programs

The project is the data resource for
Sparse PCA in Genomics
(add final title here).
It contains the complete analysis (from single-cell genomics UMI counts from
URLs to figure panels and Seurat objects), set up using the
[targets](https://docs.ropensci.org/targets/) package.

### Running Sparse PCA

To run Sparse PCA on a Seurat object using our instrumentation, the system
requirements are: Linux or macOS, for x86-64 only (ARM as well as Windows OS are
not supported); "git" command line. The former requirements are because our package does not have an external dependency
on the Julia language, which it uses, but in keeping with other R packages,
installs the dependency fresh using the Basilisk package -> Miniconda tool. The
restriction of the Anaconda "julia" package to Linux/macOS and x86-64 only is an
Anaconda known issue.

Optimal SPCA takes 4-8 GB RAM (in a helper process) in excess of the R memory
consumption of the Seurat object in question. With our presented parameters
(FindVariableFeatures nfeatures = 1000, SparsePCA npcs = 50, varnum = 8), and a
dedicated desktop PC, the wall time was 24-36 hours. The Julia dependency can
utilize several CPU cores (generally in matrix-vector multiplications), and is
somewhat leaky (memory usage will creep up, but rather than swapping to disk or
causing memory pressure, the memory usage seems to eventually stabilize). Our
Docker recipe (Dockerfile) will help with setting fixed RAM for a long-running
Sparse PCA analysis. It will also be used from a cloud instance, where it is
affordable to run our own Sparse PCA analysis.

Our instrumentation (which the end user can evaluate using a helper script,
[setupRunSparsePCA.R](./setupRunSparsePCA.R)), has the following R packages as
requirements:

* seurat **< 5.0**. We developed SPCA using Seurat and SeuratObject v4. You
  would need to load the counts matrix into Seurat **<= 4.4.0**, running
  NormalizeData, FindVariableFeatures, ScaleData, and RunSparsePCA on this
  object. Then, the object `object[["spca"]]` may be saved (using `saveRDS`) and
  may be able to be applied to a Seurat v5 object in another R installation.
* basilisk
* dplyr
* MatrixGenerics
* progress
* targets
* tarchetypes
* tibble
* withr

Next, we need to install a (nearly empty) directory named spcaFeaturesDrosophila as one of the user's R installed R packages. This is a necessary step when we create a local installation of Julia here, as other R packages do.

`tar czf spcaFeaturesDrosophila.tar.gz spcaFeaturesDrosophila/_targets.R spcaFeaturesDrosophila/DESCRIPTION spcaFeaturesDrosophila/NAMESPACE spcaFeaturesDrosophila/R/git.R && R CMD INSTALL spcaFeaturesDrosophila.tar.gz && rm -v spcaFeaturesDrosophila.tar.gz`

With these packages and spcaFeaturesDrosophila installed, evaluating the script
[setupRunSparsePCA.R](./setupRunSparsePCA.R)
should take about 10 minutes on the first run, cloning the Julia project,
invoking Miniconda several times, and invoking the Julia package manager.
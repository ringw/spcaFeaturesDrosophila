# Reproducible environment used by the Sparse PCA in Genomics project!
FROM debian:bookworm

WORKDIR /root
COPY _targets.R setupRunSparsePCA.R DESCRIPTION NAMESPACE /root/spcaFeaturesDrosophila/
COPY R/*.R /root/spcaFeaturesDrosophila/R/
COPY setupRunSparsePCA.R /root/.Rprofile
RUN apt-get update && \
  apt-get install --no-install-recommends -y \
  build-essential \
  git \
  r-cran-seurat=4.3.0-1 \
  r-cran-base64url \
  r-bioc-basilisk \
  r-cran-callr \
  r-cran-dplyr=1.0.10-1 \
  r-bioc-matrixgenerics=1.10.0-1 \
  r-cran-processx \
  r-cran-progress=1.2.2-2 \
  r-cran-ps \
  r-cran-remotes \
  r-cran-tibble \
  r-cran-withr && \
  R -e 'remotes::install_url("https://cran.r-project.org/src/contrib/Archive/targets/targets_1.4.1.tar.gz", upgrade="never")' && \
  R -e 'remotes::install_url("https://cran.r-project.org/src/contrib/Archive/tarchetypes/tarchetypes_0.7.11.tar.gz", upgrade="never")' && \
  tar czf spcaFeaturesDrosophila.tar.gz spcaFeaturesDrosophila/DESCRIPTION spcaFeaturesDrosophila/NAMESPACE spcaFeaturesDrosophila/R/git.R && \
  R CMD INSTALL spcaFeaturesDrosophila.tar.gz && \
  cd /root/spcaFeaturesDrosophila && \
  R -e 'source("setupRunSparsePCA.R", ec=T)'
WORKDIR /root/spcaFeaturesDrosophila
CMD ["/usr/bin/R"]

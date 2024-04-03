# System Requirements:
# Debian 12 (bookworm)
# Debian Stable packages should have a stable package version archived. We
# installed the packages at this point in time:
# https://snapshot.debian.org/archive/debian/20240325T091329Z/

sudo apt update
sudo apt install --no-install-recommends r-base-core=4.2.2.20221110-2 \
  r-bioc-annotationdbi=1.60.0-1 \
  r-bioc-basilisk=1.10.2+ds-4 \
  r-bioc-biocparallel=1.32.5-1 \
  r-bioc-edger=3.40.2+dfsg-1 \
  r-bioc-glmgampoi=1.10.2+dfsg-1 \
  r-bioc-limma=3.54.1+dfsg-1 \
  r-bioc-org.hs.eg.db=3.16.0-1 \
  r-bioc-scran=1.26.2+dfsg-1 \
  r-bioc-scuttle=1.8.4+dfsg-1 \
  r-cran-ape=5.7-1 \
  r-cran-argparse=2.2.2+dfsg-1 \
  r-cran-base64url=1.4-2+b1 \
  r-cran-biocmanager=1.30.20+dfsg-1 \
  r-cran-cairo=1.6-0-1 \
  r-cran-coda=0.19-4-1 \
  r-cran-cpp11=0.4.3-1 \
  r-cran-cyclocomp=1.1.0-2 \
  r-cran-emdbook=1.3.12+ds-2 \
  r-cran-fastcluster=1.2.3-2 \
  r-cran-findpython=1.0.7-1 \
  r-cran-ggforce=0.4.1-1 \
  r-cran-ggplot2=3.4.1+dfsg-1 \
  r-cran-ggpubr=0.6.0-1 \
  r-cran-ggrastr=1.0.1-2 \
  r-cran-lintr=3.0.2-1 \
  r-cran-rcppeigen=0.3.3.9.3-1 \
  r-cran-rcppparallel=5.1.6+dfsg-1 \
  r-cran-remotes=2.4.2-1 \
  r-cran-rjags=1:4-13-1 \
  r-cran-seurat=4.3.0-1 \
  r-cran-xmlparsedata=1.0.5-2

R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/targets/targets_1.4.1.tar.gz", repos=NULL, type="source")'
R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/tarchetypes/tarchetypes_0.7.11.tar.gz", repos=NULL, type="source")'
# TODO: Add other dependencies here such as HiddenMarkov
R -e 'BiocManager::install(c("apeglm", "infercnv", "RcppNumerical"), dep=F, upd=F)'
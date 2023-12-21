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
    "apeglm",
    "basilisk",
    "dplyr",
    "forcats",
    "glmGamPoi",
    "Matrix",
    "MatrixGenerics",
    "progress",
    "scran",
    "scuttle",
    "Seurat",
    "stringr",
    "tibble",
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
  tar_download(
    midgut.flybase,
    'https://ftp.flybase.org/releases/FB2019_02/precomputed_files/genes/fbgn_annotation_ID_fb_2019_02.tsv.gz',
    'fbgn_annotation_ID_2019_02.tsv.gz'
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
    midgut.metafeatures,
    build_meta_features(midgut.counts, midgut.flybase, midgut.gtf)
  ),
  tar_target(
    midgut.augment.metadata,
    add_pooled_size_factors(midgut.counts, midgut.metadata),
    format='file'
  ),
  tar_target(
    midgut.pooledSizeFactors,
    midgut_pooled_size_factors(midgut.counts, midgut.metadata)
  ),
  tar_target(
    indrop.pca,
    # esg may be included as an spca gene as it helps separate ISC from dEC, and
    # also separates EB (esg hi, reason unknown) from ISC (esg mid).
    midgut_seurat_for_technology(
        midgut.counts, midgut.metadata, midgut.metafeatures, 'inDrop',
        midgut.pooledSizeFactors, spca_genes='esg')
    # Now call the pca clusters.
    %>% FindNeighbors(dims=1:9) %>% FindClusters(res=1.04, random.seed=1)
    %>% AddMetaData(
      Idents(.) %>% fct_recode(
        ISC='6', EB='7',
        `EC-like`='0', `EC-like2`='2',
        `EC-like3`='20',
        `EC.anterior1`='1', `EC.anterior2`='5',
        `EC.anterior3`='8', EC.anterior4='9',
        EC.anterior5='25',
        `EE.Ast`='3',
        cardia1='18', cardia2='23',
        dEC1='4', dEC2='11',
        `copper/iron`='10',
        EC.meso='12',
        LFC1='13', LFC2='19',
        others.1='14', others.2='18', others.3='22',
        EC.posterior1='15', EC.posterior2='16',
        EC.posterior3='17', EC.posterior4='21',
        EC.posterior5='24'
      ) %>% fct_relevel('ISC', 'EB', 'EC.anterior1'),
      'pca_clusters'
    ) %>% midgut_classify_cell_types('pca_clusters')
  ),
  tar_target(
    indrop.spca.dimreduc,
    RunSparsePCA(
      indrop.pca, 'covar', varnum=8, npcs=60, eigen_gap=0.05, search_cap=500000
    )[['spca']],
    # We are not going to change the construction of the 'covar' matrix, so
    # don't rerun this expensive step.
    cue=tar_cue('never')
  ),
  tar_target(
    indrop.sct,
    SCTransform(
      indrop.pca %>% subset(
        # rRNA contamination found using the inDrop technology. In SCTransform,
        # the rRNA genomic features might 
        features = rownames(.[['RNA']]) %>% subset(!grepl('rRNA', .))
      ),
      do.correct.umi = F
    )[['SCT']]
  ),
  tar_target(
    indrop.sct.pca,
    RunPCA(indrop.sct, verb=F, assay='SCT')
  ),
  tar_target(
   indrop,
   spca_with_centered_umap(indrop.pca, indrop.spca.dimreduc, dims=1:50)
   %>% FindNeighbors(dims=1:36, red='spca')
   %>% FindClusters(res=1.28, random.seed=1)
   %>% AddMetaData(
     Idents(.) %>% fct_recode(
       ISC='14', EB='5'
     ) %>% fct_relevel('ISC', 'EB', '1'),
     'spca_clusters'
   ) # %>% midgut_classify_cell_types('spca_clusters')
  ),
  tar_target(
    indrop.glm,
    build_glms(indrop, midgut.pooledSizeFactors, c('pca_clusters', 'spca_clusters'))
  ),
  tar_target(
    indrop.present.genes,
    build_present_gene_list(indrop, c('pca_clusters', 'spca_clusters'))
  ),
  tar_target(
    indrop.deg,
    build_de_data(indrop.glm, indrop.present.genes)
  )
)

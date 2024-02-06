# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(dplyr)
library(stringr)
library(targets)
library(tarchetypes)
library(tibble)

# Set target options:
tar_option_set(
  packages = c(
    "AnnotationDbi",
    "apeglm",
    "basilisk",
    "colorspace",
    "dplyr",
    "forcats",
    "ggforce",
    "ggnewscale",
    "ggplot2",
    "ggpubr",
    "ggrastr",
    "glmGamPoi",
    "infercnv",
    "Matrix",
    "MatrixGenerics",
    "org.Hs.eg.db",
    "pillar",
    "processx",
    "progress",
    "reshape2",
    "scales",
    "scran",
    "scuttle",
    "Seurat",
    "stringr",
    "tibble",
    "viridis",
    "withr"
  )
)

options(clustermq.scheduler = "multicore")

tar_source()
# source("other_functions.R") # Source other scripts as needed.

midgut_figures = list(
  tar_target(
    fig.indrop.pca,
    save_figure(
      "figure/Midgut/UMAP-1-PCA.pdf",
      plot_indrop_pca(indrop),
      width = 5, height = 4
    ),
    format = "file"
  ),
  tar_target(
    fig.indrop.spca,
    save_figure(
      "figure/Midgut/UMAP-1-SPCA.pdf",
      plot_indrop_spca(indrop),
      width = 5, height = 4
    ),
    format = "file"
  )
)

# Patient identifiers.
acc_individual <- paste0("ACC", c(2,5,7,15,19,21,22))

acc_mutually_exclusive <- tribble(
  ~feature1, ~feature2,
  "NOTCH3", "DLL1",
  "N3+", "N1N2",
  "targets", "ligands",
  "SPARSE_14", "SPARSE_1",
  "SPARSE_26", "SPARSE_3",
  "SPARSE_11", "SPARSE_12"
)
acc_mutually_exclusive$name <- acc_mutually_exclusive %>% with(
  paste(str_replace(feature1, "\\+", ""), str_replace(feature2, "\\+", ""), sep="_")
)

acc_feature_loadings = bind_rows(
  acc_mutually_exclusive %>%
    subset(select=c(feature1, feature2)) %>%
    apply(2, \(v) v %>% grep("^SPARSE_", ., val=T) %>% data.frame(feature = .), simplify = FALSE) %>%
    setNames(c("magma", "viridis")),
  .id = "viridis_option"
) %>% within(viridis_option <- as.character(viridis_option))
acc_feature_loadings = bind_rows(
  acc_mutually_exclusive %>%
    subset(select=c(feature1, feature2)) %>%
    apply(2, \(v) v %>% grep("^SPARSE_", ., val=T) %>% data.frame(feature = .), simplify = FALSE)
) %>% within(viridis_option <- "viridis")

# New features: Plot a single feature, use the bounding box of ACC and CAF
acc_features <- tribble(
  ~feature, ~max_scale, ~annotations,
  # "DLL1", 5, quote(list(c("tile", 0, 1, width=4, height=3), c("tile", 0, -7, width=4, height=5.5))),
  "DLL1", 8, quote(list(list(geom="tile", x=0, y=1, width=4, height=3))),
  "NOTCH3", 8, NULL,
  "SPARSE_11", 8, NULL,
  "SPARSE_26", NA, NULL,
  "SPARSE_12", NA, NULL,
  "SPARSE_3", NA, NULL,
  "FN1", NA, NULL
)
acc_features_max_scale <- 8
acc_query <- tribble(
  ~ident, ~feature,
  "ACC", "NOTCH3",
  "ACC", "DLL1",
  "CAF", "NOTCH3",
  "CAF", "DLL1"
) # %>%
  # left_join(
    # acc_features %>%
      # subset(select=c(feature, max_scale)),
    # by = "feature"
  # )

acc_figures = list(
  tar_target(
    fig.acc.annotated,
    save_figure("figure/ACC/ACC-Annotation.pdf",
    acc_annotated_figure(acc, "aneuploidy", guide = "Copy Number", limits = c(0, 0.2), oob_squish = TRUE), width=4, height=3),
    format = "file"
  ),
  tar_target(
    fig.acc.arranged,
    save_figure("figure/ACC/ACC-Preview.pdf", acc_arrange_figure(acc), width=20, height=4),
    format = "file"
  ),
  tar_map(
    data.frame(individual = acc_individual),
    tar_target(
      fig.acc.individual,
      save_figure(
        paste0(
          "figure/ACC/ACC-Preview-",
          individual,
          ".pdf"
        ),
        acc_arrange_figure(acc, individual), width=20, height=4
      ),
      format = "file"
    )
  ) %>% append(
    list(
      tar_combine(
        fig.acc.allindividuals,
        .,
        command = {
          filename <- "figure/ACC/ACC-Preview-Individuals.pdf"
          run(
            "pdfunite",
            c(
              !!!.x,
              filename
            )
          )
          filename
        },
        format = "file"
      )
    )
  ),
  tar_map(
    acc_mutually_exclusive,
    names = name,
    tar_target(
      fig.acc.blend,
      save_figure(
        paste0(
          "figure/ACC/ACC-Features-", feature1, "-", feature2, ".pdf"
        ),
        acc_blend_features_inset(acc, feature1, feature2, show.legend = F),
        width = 2, height = 4
      )
    ),
    tar_target(
      fig.acc.blend.legend,
      save_figure(
        paste0(
          "figure/ACC/ACC-Features-", feature1, "-", feature2, "-Legend.pdf"
        ),
        acc_blend_features_inset(acc, feature1, feature2) %>% get_legend,
        width = 0.75, height = 3.1
      ),
      format = "file"
    )
  ),
  tar_map(
    acc_feature_loadings,
    names = feature,
    tar_target(
      fig.acc.loadings,
      save_figure(
        paste0("figure/ACC/ACC-Feature-Loadings-", feature, ".pdf"),
        (
          spc_tile_plot(
            acc.spca.dimreduc, feature, option=viridis_option, end=0.5,
            fontface = "plain"
          )
          + theme(legend.position = "bottom", legend.margin = margin(r = 5))
          + guides(
            fill = guide_colorbar(barwidth = 3, barheight = 0.75, title = "")
          )
        ),
        width = 0.9,
        height = 3.5
      ),
      format = "file"
    )
  ),
  tar_map(
    acc_features,
    names = feature,
    tar_target(
      fig.acc.feature,
      save_figure(
        paste0("figure/ACC/ACC-Feature-Inset-", feature, ".pdf"),
        nonneg_feature_plot_annotate(
          acc, feature, max_scale = 8, annotations = NULL
        ) + theme(legend.position = "none"),
        width = 2,
        height = 4
      ),
      format = "file"
    ),
    tar_target(
      fig.acc.feature.caf,
      save_figure(
        paste0("figure/ACC/ACC-Feature-Inset-CAFMarked-", feature, ".pdf"),
        nonneg_feature_plot_annotate(
          acc, feature, max_scale = 8, annotations = NULL
        ) + theme(
          legend.position = "none",
          plot.margin = margin(r = 25.5, l = 5.5, t = 5.5, b = 5.5)
        ) + annotate(
          "rect", xmin=-2.25, ymin=-9.75, xmax=6.25, ymax=-4.25, fill = "transparent",
          color = "magenta"
        ) + annotate(
          "segment", c(6.25, 6.25), c(-9.75, -4.25), xend=c(9, 9), yend=c(-11, 8.5),
          color = "magenta"
        ),
        width = 2.28,
        height = 4
      ),
      format = "file"
    ),
    tar_target(
      fig.acc.feature.umap.caf,
      save_figure(
        paste0("figure/ACC/ACC-Feature-CAF-UMAP-", feature, ".pdf"),
        nonneg_feature_plot_caf(
          caf, feature, max_scale = 8
        ) + theme(
          legend.position = "none",
        ),
        width = 2,
        height = 3
      ),
      format = "file"
    )
  ),
  tar_map(
    acc_query,
    names = all_of(c("feature", "ident")),
    tar_target(
      fig.acc.feature.query,
      save_figure(
        paste0("figure/ACC/ACC-Feature-Query-", feature, "-", ident, ".pdf"),
        nonneg_feature_plot_query(acc %>% acc_add_components_from_umap, feature, max_scale = 8, ident)
          + theme(legend.position = "none"),
        width = 2,
        height = 4
      ),
      format = "file"
    )
  ),
  tar_target(
    fig.acc.feature.legend,
    save_figure(
      "figure/ACC/ACC-Feature-Legend.pdf",
      (
        ggplot(data.frame(x=0, y=0, color=c(0, 8)), aes(x, y, color=color))
          + geom_point()
          + scale_color_gradientn(
            colors = acc_nonneg_feature_plot_gradient
          )
          + labs(color = "LogNormalize")
      ) %>% get_legend,
      width = 1,
      height = 2
    )
  ),
  tar_target(
    fig.acc.individual.discrete,
    save_figure(
      "figure/ACC/ACC-Individuals.pdf",
      acc_dim_plot_individual(acc, "individual"),
      width = 4,
      height = 6
    )
  ),
  tar_target(
    fig.acc.dotplot,
    save_figure(
      "figure/ACC/ACC-Dot-Plot.pdf",
      acc_dot_plot(
        acc.glm,
        c("MYB", "DLL1", "CLDN3", "NOTCH3", "FN1", "MMP11", "CXCL12", "CXCL2", "ACTA2", "PAX7")
      ),
      width = 4,
      height = 3
    )
  ),
  tar_target(
    fig.acc.cafident,
    save_figure(
      "figure/ACC/ACC-CAF-UMAP.pdf",
      idents_plot_caf(caf),
      width = 2,
      height = 2.5
    )
  )
)

# Replace the target list below with your own:
list(
  tar_download(
    midgut.counts,
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120537/suppl/GSE120537%5Fcounts.csv.gz',
    'GSE120537_counts.csv.gz',
    cue=tar_cue('never')
  ),
  tar_download(
    midgut.metadata,
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120537/suppl/GSE120537%5Fmetadata.csv.gz',
    'GSE120537_metadata.csv.gz',
    cue=tar_cue('never')
  ),
  tar_download(
    midgut.gtf,
    'https://ftp.flybase.org/releases/FB2019_02/dmel_r6.27/gtf/dmel-all-r6.27.gtf.gz',
    'dmel-all-r6.27.gtf.gz',
    cue=tar_cue('never')

  ),
  tar_download(
    midgut.flybase,
    'https://ftp.flybase.org/releases/FB2019_02/precomputed_files/genes/fbgn_annotation_ID_fb_2019_02.tsv.gz',
    'fbgn_annotation_ID_2019_02.tsv.gz',
    cue=tar_cue('never')

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
        midgut.pooledSizeFactors, spca_genes='esg') %>%
      # Now call the pca clusters.
      FindNeighbors(dims=1:9) %>%
      FindClusters(res=1.04, random.seed=1) %>%
      AddMetaData(
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
      ) %>%
      midgut_classify_cell_types('pca_clusters')
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
        # features are not unit-scaled, and the rRNA features might have an
        # unwanted amount of variance (ideally little to none).
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
    spca_with_centered_umap(
      indrop.pca, indrop.spca.dimreduc, dims=1:50,
      # Flip the UMAP horizontally and vertically.
      umap_transform=-diag(2)
    ) %>%
      FindNeighbors(dims=1:36, red='spca') %>%
      FindClusters(res=1.28, random.seed=1) %>%
      AddMetaData(
        Idents(.) %>% fct_recode(
          ISC='14', EB='5',
          `EC-like`='0',
          EC.anterior='1',
          `copper/iron`='2',
          EC2='3',
          dEC1='4',
          EC3='6',
          EC4='7',
          EE1='8',
          LFC1='9',
          dEC2='10',
          EE2='11',
          EC5='12',
          others.1='13',
          cardia='15',
          EC6='16',
          EC7='17',
          `EC-like2`='18',
          LFC2='19',
          others.2='20',
          EC8='21',
          `EC-like3`='22',
          others.3='23',
          `EC-like4`='24',
          EC9='25'
        ) %>% fct_relevel('ISC', 'EB'),
        'spca_clusters'
      ) %>%
      midgut_classify_cell_types('spca_clusters') %>%
      AddMetaData(
        compute_sct_clusters(indrop.sct.pca),
        'sct_clusters'
      ) %>%
      midgut_classify_cell_types('sct_clusters')
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
  ),
  tar_target(
    indrop.misc,
    build_midgut_misc_stats(indrop, indrop.sct.pca)
  ),

  # Adenoid Cystic Carcinoma sample
  tar_download(
    acc.counts,
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210171/suppl/GSE210171%5Facc%5Fscrnaseq%5Fcounts.txt.gz',
    'GSE210171_acc_scrnaseq_counts.txt.gz',
    cue=tar_cue('never')
  ),
  tar_target(acc.rownames, make_rownames.acc(acc.counts)),
  tar_target(
    acc.gene.order,
    write_gene_order.acc(
      acc.counts, acc.rownames, "GSE210171_acc_gene_order.tsv"
    ),
    format = "file"
  ),
  tar_target(acc.rna, preprocess_acc(acc.counts, acc.rownames)),
  tar_target(
    acc.spca.dimreduc,
    RunSparsePCA(
      acc.rna, 'covar', varnum=10, npcs=50, eigen_gap=0.05, search_cap=100000
    )[['spca']],
    # We are not going to change the construction of the 'covar' matrix, so
    # don't rerun this expensive step.
    cue=tar_cue('never')
  ),
  tar_target(
    acc.spca,
    process_acc_spca(acc.rna, acc.spca.dimreduc)
  ),
  tar_target(
    acc.annotations.pca,
    write_seurat_column(acc.spca, "pca_coarse", "GSE210171_acc_pca_coarse.tsv"),
    format = "file"
  ),
  tar_target(
    acc.infercnv.pca,
    {
      output_path <- "acc_infercnv_pca"
      if (file.exists(output_path))
        file.rename(output_path, paste0(output_path, "~"))
      with_options(
        list(scipen=100),
        infercnv::run(
          infercnv::CreateInfercnvObject(
            acc.rna[['RNA']]@counts,
            acc.gene.order,
            acc.annotations.pca,
            # FindMarkers: cluster '0' is TP63+, while cluster '1' is not.
            ref_group_names='1'
          ),
          out_dir = output_path,
          HMM = TRUE
        )
      )
      output_path
    },
    format = "file",
    cue = tar_cue("never")
  ),
  tar_target(
    acc,
    acc.spca %>%
      infercnv::add_to_seurat(infercnv_output_path = acc.infercnv.pca)
  ),
  acc_figures,
  tar_combine(acc.figures, acc_figures)
)

# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(dplyr)
library(grDevices)
library(stringr)
library(targets)
library(tarchetypes)
library(tibble)

# Set font for figures.
pdfFonts(sans = pdfFonts()$Helvetica)
postscriptFonts(sans = postscriptFonts()$Helvetica)

# Set target options:
tar_option_set(
  packages = c(
    "AnnotationDbi",
    "basilisk",
    "colorspace",
    "cowplot",
    "dplyr",
    "forcats",
    "ggforce",
    "ggnewscale",
    "ggplot2",
    "ggpubr",
    "ggrastr",
    "grDevices",
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
    "tidyselect",
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
      "figure/Midgut/UMAP-5-SPCA.pdf",
      plot_indrop_spca(indrop),
      width = 5, height = 4
    ),
    format = "file"
  ),
  tar_target(
    fig.indrop.pca.new,
    save_figure(
      "figure/Midgut/Feature-PCA-Clustering.pdf",
      plot_indrop_pca(indrop),
      width = 5, height = 4
    ),
    format = "file"
  ),
  tar_target(
    fig.indrop.spca.new,
    save_figure(
      "figure/Midgut/Feature-SPCA-Clustering.pdf",
      plot_indrop_spca(indrop),
      width = 5, height = 4
    ),
    format = "file"
  ),
  tar_target(
    fig.indrop.legend,
    save_figure(
      "figure/Midgut/Feature-Legend-Clustering.pdf",
      plot_midgut_legend(indrop),
      width = 4, height = 0.6
    ),
    format = "file"
  ),
  tar_target(
    fig.indrop.pca.ec,
    save_figure(
      "figure/Midgut/UMAP-Feature-PCA-betaTry.pdf",
      plot_midgut_feature(indrop, "PCA", "umap", "betaTry"),
      width = 5, height = 4
    ),
    format = "file"
  ),
  tar_target(
    fig.indrop.spca.ec,
    save_figure(
      "figure/Midgut/UMAP-Feature-SPCA-betaTry.pdf",
      plot_midgut_feature(indrop, "SPCA", "umap.spca", "betaTry"),
      width = 5, height = 4
    ),
    format = "file"
  ),
  tar_target(
    indrop.feature.list,
    c(
      # For calling ISC and EB cluster labels:
      "esg", "pros", "E(spl)mbeta-HLH", "N", "Dl",
      # aEC markers:
      "betaTry", "Amy-p", "Jon65Aiv",
      # mEC
      "CG13315", "Vha16-1",
      # pEC
      "LManVI", "Gs2", "Mur29B",
      # Novel markers (DE-CI)
      "His3.3A",
      # Others (not called yet as dEC markers?)
      "fmt", "trol",
      "SPARSE_1", "SPARSE_26"
    )
  ),
  tar_target(indrop.feature.printing.limits, list(`E(spl)mbeta-HLH`=c(0,4.5), SPARSE_26=c(0,4.5))),
  tar_target(
    indrop.pca.stemlike.inset,
    annotate("rect", xmin=-3.8, ymin=7.7, xmax=-0.2, ymax=2.5, fill="transparent", color="black", linewidth=0.5)
  ),
  tar_target(
    indrop.spca.stemlike.inset,
    annotate("rect",xmin=-6.5, ymin=4.9, xmax=-2.5, ymax=2.3, fill="transparent", color="black", linewidth=0.5)
  ),
  tar_target(
    indrop.feature.printing.annot,
    list(
      PCA=list(`E(spl)mbeta-HLH`=indrop.pca.stemlike.inset, SPARSE_26=indrop.pca.stemlike.inset),
      SPCA=list(`E(spl)mbeta-HLH`=indrop.spca.stemlike.inset, SPARSE_26=indrop.spca.stemlike.inset)
    )
  ),
  tar_target(
    indrop.feature.printing,
    tribble(~name, ~embedding, "PCA", "umap", "SPCA", "umap.spca")
  ),
  tar_target(
    fig.indrop.feature,
    save_figure(
      paste0("figure/Midgut/Feature-", indrop.feature.printing %>% pull(name), "-", indrop.feature.list, ".pdf"),
      plot_midgut_feature(
        indrop,
        bg_color=indrop.feature.printing %>% pull(name),
        indrop.feature.printing %>% pull(embedding),
        indrop.feature.list,
        limits=indrop.feature.printing.limits[[indrop.feature.list]]
      ) + indrop.feature.printing.annot[[
        indrop.feature.printing %>% pull(name)
      ]][[
        indrop.feature.list
      ]],
      width = 5, height = 4
    ),
    format = "file",
    pattern = cross(indrop.feature.list, indrop.feature.printing)
  ),
  tar_target(
    indrop.inset.theme,
    theme(
      aspect.ratio = 0.75,
      axis.text = element_text(size = 12)
    )
  ),
  tar_target(
    fig.indrop.pca.espl.inset,
    save_figure(
      "figure/Midgut/Feature-PCA-E(spl)mbeta-HLH-Inset.pdf",
      plot_midgut_feature(
        indrop,
        bg_color="PCA",
        "umap",
        "E(spl)mbeta-HLH",
        limits=indrop.feature.printing.limits[[1]]
      ) + scale_x_continuous(
        name=NULL, breaks=pretty_breaks(3)
      ) + scale_y_continuous(
        name=NULL, breaks=pretty_breaks(3)
      ) + coord_cartesian(
        c(indrop.pca.stemlike.inset$data$xmin, indrop.pca.stemlike.inset$data$xmax),
        c(indrop.pca.stemlike.inset$data$ymax, indrop.pca.stemlike.inset$data$ymin)
      ) + indrop.inset.theme,
      2,
      1.5
    )
  ),
  tar_target(
    fig.indrop.spca.espl.inset,
    save_figure(
      "figure/Midgut/Feature-SPCA-E(spl)mbeta-HLH-Inset.pdf",
      plot_midgut_feature(
        indrop,
        bg_color="SPCA",
        "umap.spca",
        "E(spl)mbeta-HLH",
        limits=indrop.feature.printing.limits[[1]]
      ) + scale_x_continuous(
        name=NULL, breaks=pretty_breaks(3)
      ) + scale_y_continuous(
        name=NULL, breaks=c(3,4)
      ) + coord_cartesian(
        c(indrop.spca.stemlike.inset$data$xmin, indrop.spca.stemlike.inset$data$xmax),
        c(indrop.spca.stemlike.inset$data$ymax, indrop.spca.stemlike.inset$data$ymin)
      ) + indrop.inset.theme,
      2,
      1.5
    )
  ),
  tar_file(
    fig.indrop.spca.spc26.inset,
    save_figure(
      "figure/Midgut/Feature-SPCA-SPARSE_26-Inset.pdf",
      plot_midgut_feature(
        indrop,
        bg_color="SPCA",
        "umap.spca",
        "SPARSE_26",
        limits=indrop.feature.printing.limits[[1]]
      ) + scale_x_continuous(
        name=NULL, breaks=pretty_breaks(3)
      ) + scale_y_continuous(
        name=NULL, breaks=c(3,4)
      ) + coord_cartesian(
        c(indrop.spca.stemlike.inset$data$xmin, indrop.spca.stemlike.inset$data$xmax),
        c(indrop.spca.stemlike.inset$data$ymax, indrop.spca.stemlike.inset$data$ymin)
      ) + indrop.inset.theme,
      2,
      1.5
    )
  ),
  tar_file(
    fig.indrop.legend.betaTry,
    save_figure(
      "figure/Midgut/Feature-Legend-betaTry.pdf",
      plot_midgut_feature_legend(indrop, "betaTry"),
      width = 2.65, height = 0.6
    )
  ),
  tar_file(
    fig.indrop.legend.thumbnail,
    save_figure(
      "figure/Midgut/Feature-Legend-Thumbnail.pdf",
      plot_midgut_feature_legend(
        indrop, limits=c(0,4.5),
        legend.direction="vertical",
        legend.name=NULL
      ),
      width = 0.5, height = 1.5
    )
  ),
  tar_file(
    fig.indrop.model.background.legend,
    save_figure(
      "figure/Midgut/Model-Background-Legend.pdf",
      plot_midgut_model_background_legend(),
      width = 2, height = 0.6
    )
  ),
  tar_target(
    midgut.cdf.cluster.list, c("ISC", "EB", "dEC", "EC", "EC-like", "EE")
  ),
  tar_target(
    indrop.violin.pcs,
    tribble(
      ~name, ~value, ~group_by,
      "SPC1", "SPARSE_1", "spca_classif",
      "SPC26", "SPARSE_26", "spca_classif",
      "PC3", "PC_3", "pca_classif",
      "PC7", "PC_7", "pca_classif"
    )
  ),
  tar_target(
    fig.indrop.violin,
    save_figure(
      paste0("figure/Midgut/Indrop-Violin-", indrop.violin.pcs$name, ".pdf"),
      tiny_violin_plot(
        indrop, indrop.violin.pcs$group_by,
        midgut.cdf.cluster.list, indrop.violin.pcs$value
      ),
      width=4.1, height=2.5
    ),
    pattern = map(indrop.violin.pcs),
    format = "file"
  ),
  tar_target(
    fig.indrop.thumbnail.violin,
    save_figure(
      paste0("figure/Midgut/Indrop-Thumbnail-Violin-", indrop.violin.pcs$name, ".pdf"),
      tiny_violin_plot(
        indrop %>%
          AddMetaData(.$pca_classif %>% recode(`EC-like`=""), "pca_classif") %>%
          AddMetaData(.$spca_classif %>% recode(`EC-like`=""), "spca_classif"),
        indrop.violin.pcs$group_by,
        midgut.cdf.cluster.list %>% replace(. == "EC-like", ""), indrop.violin.pcs$value
      ) + theme(
        plot.margin = margin(2.5, 1, 1, 1), aspect.ratio = 0.618,
        axis.text = element_text(size = 20)
      ),
      width=3.5, height=2.3
    ),
    pattern = map(indrop.violin.pcs),
    format = "file"
  ),
  tar_map(
    tribble(~name, ~value, "SPC1", "SPARSE_1", "SPC26", "SPARSE_26"),
    names = name,
    tar_file(
      fig.indrop.cdf,
      save_figure(
        paste0("figure/Midgut/Indrop-Violin-CDF-", name, ".pdf"),
        plot_spca_cdf(
          indrop,
          "spca_classif",
          midgut.cdf.cluster.list,
          value
        ),
        width = 4.1, height = 2.5
      )
    ),
    tar_file(
      fig.indrop.thumbnail.cdf,
      save_figure(
        paste0("figure/Midgut/Indrop-Thumbnail-Violin-CDF-", name, ".pdf"),
        plot_spca_cdf(
          indrop,
          "spca_classif",
          midgut.cdf.cluster.list,
          value
        ) + theme(
          plot.margin = margin(5.5, 1, 1, 1), aspect.ratio = 0.7,
          axis.text = element_text(size = 20)
        ),
        width = 3.5, height = 2.3
      )
    )
  ),
  tar_target(
    fig.indrop.deg.pca,
    save_figure(
      "figure/Midgut/Indrop-DEG-PCA.pdf",
      plot_arrange_deg_model_color_panels(indrop.deg$pca_clusters, limits=c(-4.5, 4.5))
      + annotate(
        "text",
        -5.05, 13.5,
        label="A",
        fontface="bold",
        size=2
      )
      + theme(
        aspect.ratio = 0.4, legend.position = "none",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.title.x = element_text(margin = margin(0, 0, 0, 0)),
        plot.margin = margin(0, 0, -10, 0)
      ),
      4.5, 2.12,
      device = CairoPDF
    ),
    format = "file",
    packages = tar_option_get("packages") %>% c("Cairo")
  ),
  tar_target(
    fig.indrop.deg.spca,
    save_figure(
      "figure/Midgut/Indrop-DEG-SPCA.pdf",
      plot_arrange_deg_model_color_panels(indrop.deg$spca_clusters, limits=c(-2.6, 2.6))
      + annotate(
        "text",
        -2.92, 13.5,
        label="B",
        fontface="bold",
        size=2
      )
      + theme(
        aspect.ratio = 0.4,
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.title.x = element_text(margin = margin(0, 0, 0, 0)),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(8, "pt"),
        plot.margin = margin(0, 0, -5, 0),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-6, -6, -6, -6)
      ),
      4.5, 2.36,
      device = CairoPDF
    ),
    format = "file",
    packages = tar_option_get("packages") %>% c("Cairo")
  ),
  tar_target(
    fig.indrop.deg.legend,
    save_figure(
      "figure/Midgut/DEG-Legend.pdf",
      cowplot::get_legend(
        plot_midgut_de_panel(
          tibble(
            rowname="a", displayName="a", log2FoldChange=2, lfcSE=1, q_val=1
          )
        )
        + guides(color = guide_legend(nrow = 2))
        + theme(
          legend.position = "bottom",
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.75, "lines")
        )
      ),
      5, 1
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
  "SPARSE_14", 8, NULL,
  "SPARSE_26", NA, NULL,
  "SPARSE_12", NA, NULL,
  "SPARSE_3", NA, NULL,
  "FN1", NA, NULL,
  "SPARSE_29", 8, NULL
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
    acc_annotated_figure(acc, "aneuploidy", guide = "CNV", limits = c(0, 0.2), oob_squish = TRUE)
      + theme(text = element_text(size=16), aspect.ratio = 0.978),
    width=4, height=3),
    format = "file"
  ),
  tar_target(
    fig.acc.annotated.narrow,
    save_figure(
      "figure/ACC/ACC-Annotation-Narrow.pdf",
      acc_annotated_figure(acc, "aneuploidy", guide = "CNV Score", limits = c(0, 0.2), oob_squish = TRUE)
      + theme(text = element_text(size=14), legend.position = "none"),
      width=2.85, height=3
    ),
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
        ) + scale_color_viridis_c(option="magma", end = 0.9, limits=c(0,8), oob=scales::squish)
        + theme(legend.position = "none"),
        width = 2,
        height = 3
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
    cue=tar_cue('never'),
    packages="withr"
  ),
  tar_target(
    OptimalSPCADepot,
    julia_pkg_OptimalSPCADepot(),
    format='file',
    cue=tar_cue('never'),
    packages=c("basilisk", "withr")
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
    midgut_seurat_for_technology(
        midgut.counts, midgut.metadata, midgut.metafeatures, 'inDrop',
        midgut.pooledSizeFactors) %>%
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
          EC.meso='12', LFC1='13', LFC2='19',
          others.1='14', others.2='22',
          EC.posterior1='15', EC.posterior2='16',
          EC.posterior3='17', EC.posterior4='21',
          EC.posterior5='24'
        ) %>% fct_relevel('ISC', 'EB', 'EC.anterior1'),
        'pca_clusters'
      ) %>%
      midgut_classify_cell_types('pca_clusters') %>%
      # We found that PC3 has a heavy right tail containing most of the EC
      # cells. It has a long tail of other EC cells, but that tail is low in
      # density and we will put that tail in the negative direction.
      replace_pca_embedding_feature(
        "PC_3", \(v) -v * sign(v)[which.max(abs(v))]
      ) %>%
      # EB (the extreme values of PC7) should be positive.
      replace_pca_embedding_feature(
        "PC_7", \(v) v * sign(v)[which.max(abs(v))]
      )
  ),
  tar_target(
    indrop.spca.param.k.input,
    tibble(k = 1:20)
  ),
  tar_target(
    indrop.spca.param.k,
    indrop.spca.param.k.input %>%
      rowwise %>%
      mutate(
        dimreduc = indrop.pca %>%
          seurat_spca("covar", varnum=k, npcs=2, eigen_gap=0.01, search_cap=500000) %>%
          list
      ),
    pattern = map(indrop.spca.param.k.input)
  ),
  tar_target(
    indrop.spca.k4,
    indrop.pca %>%
      seurat_spca("covar", varnum=4, npcs=10, eigen_gap=0.01, search_cap=200000)
  ),
  tar_target(
    indrop.spca.k5,
    indrop.pca %>%
      seurat_spca("covar", varnum=5, npcs=10, eigen_gap=0.01, search_cap=200000)
  ),
  tar_target(
    indrop.pca.quickCluster,
    indrop.pca %>%
      FindNeighbors(dims = 1:findPC(indrop.pca[['pca']]@stdev)) %>%
      FindClusters(res = 0.1) %>%
      Idents,
    packages = "findPC"
  ),
  tar_target(
    indrop.spca.modsel.features,
    rowAnys(
      indrop.spca.param.k %>%
        pull(dimreduc) %>%
        sapply(\(dimreduc) dimreduc@feature.loadings[, "SPARSE_1"]) %>%
        `!=`(0),
      useNames = TRUE
    ) %>%
      which %>%
      names
  ),
  tar_target(
    indrop.spca.select.k,
    indrop.spca.param.k %>%
      rowwise %>%
      mutate(
        logLik = score_univariate_mixture_model(
          indrop.pca[["RNA"]],
          indrop.spca.modsel.features,
          dimreduc@feature.loadings[, "SPARSE_1"],
          indrop.pca.quickCluster
        ),
        .keep = "unused"
      ),
    packages = "mvtnorm"
  ),
  tar_target(
    indrop.spca.models,
    replicate_spca_model_parallel(
      indrop.pca, memory_cgroups, seed=0:3,
      varnum=8, npcs=50, eigen_gap=0.001, search_cap=500000,
      do.correct.elbow = TRUE
    ),
    packages = "future",
    # 'covar' matrix has not changed in the 'indrop.pca' object, so do not rerun
    # spca models target which has an ETA of 12-24 hours.
    cue = tar_cue("never")
  ),
  tar_target(
    indrop.pca.validation.model.pattern,
    as.list(1:4)
  ),
  tar_target(
    indrop.pca.validation.models,
    irlba(
      indrop.pca[["RNA"]]@scale.data, nv=25,
      # Ensure that we get a random start to the irlba algorithm.
      v = rnorm(n = ncol(indrop.pca))
    ),
    pattern = map(indrop.pca.validation.model.pattern),
    iteration = "list",
    packages = "irlba"
  ),
  tar_target(
    indrop.pca.replicates.identical,
    full_join(
      tibble(pca1 = 1:4, model1 = indrop.pca.validation.models),
      tibble(pca2 = 1:4, model2 = indrop.pca.validation.models),
      by=character()
    ) %>%
      rowwise %>%
      mutate(
        all.equal.v = all.equal(
          diag(nrow=25),
          # Use abs to correct for svd being identical up to the sign.
          abs(t(model1$v) %*% model2$v)
        ),
        all.equal.d = all.equal(model1$d, model2$d)
      ) %>%
      subset(select=-c(model1, model2))
  ),
  tar_target(
    indrop.validation.neighbors,
    bind_rows(
      pca=tibble(replicate = factor(1:4)) %>%
        rowwise %>%
        summarise(
          replicate,
          # FindNeighbors returns nn and snn. Keep only the list element snn in
          # a new list-type column.
          snn = FindNeighbors(
            indrop.pca[["pca"]]@cell.embeddings[, 1:9] %*%
              diag(ifelse(runif(n = 9) < 0.5, -1, 1)),
            nn.method = "rann",
            verb = FALSE
          )["snn"]
        ),
      spca=tibble(replicate = factor(1:4), obj = indrop.spca.models) %>%
        rowwise %>%
        summarise(
          replicate,
          snn = FindNeighbors(
            obj@cell.embeddings[, seq(obj@stdev > sqrt(2.5))],
            nn.method = "rann",
            verb = FALSE
          )["snn"]
        ),
      .id = "model"
    )
  ),
  tar_target(
    indrop.validation.clusters,
    full_join(
      indrop.validation.neighbors,
      tibble(random.seed = 0:9),
      by=character()
    ) %>%
      rowwise %>%
      summarise(
        model,
        replicate,
        random.seed,
        resolution_search_lower = find_isc_clustering_res(indrop.pca, snn, random.seed = random.seed, by = 0.01) %>%
          binsearchChooseX %>%
          `/`(100),
        resolution_search_upper = find_isc_clustering_res(
          indrop.pca, snn, random.seed = random.seed, res_range = 100:800, by = 0.01, intended_count=3
        ) %>%
          binsearchChooseX %>%
          `/`(100),
        # best_resolution = (resolution_search_lower + resolution_search_upper) / 2,
        clustering = list(
          clusterAndLabelIscEb(indrop.pca, snn, res = resolution_search_lower, random.seed = random.seed)
        )
      )
  ),
  tar_map(
    tibble(
      evaluation = c("adj.rand.index", "cor"),
      FN = rlang::syms(c("guessRandIndex", "corIscEb"))
    ),
    names = evaluation,
    tar_target(
      indrop.eval.clusters,
      indrop.validation.clusters %>%
        group_by(model, replicate) %>%
        summarise(
          model,
          replicate,
          random.seed = unique(random.seed),
          score = full_join(
            tibble(random.seed, fct1 = clustering),
            tibble(random.seed2 = random.seed, fct2 = clustering),
            by = character()
          ) %>%
            filter(random.seed != random.seed2) %>%
            rowwise %>%
            mutate(score = FN(fct1, fct2)) %>%
            group_by(random.seed) %>%
            summarise(score = mean(score)) %>%
            pull(score)
        ),
      packages = "pdfCluster"
    ),
    tar_target(
      indrop.eval.clusters.replicated,
      indrop.validation.clusters %>%
        group_by(model) %>%
        summarise(
          model,
          replicate,
          score = full_join(
            tibble(replicate, random.seed, fct1 = clustering),
            tibble(fct2 = clustering),
            by=character()
          ) %>%
            rowwise %>%
            mutate(score = FN(fct1, fct2)) %>%
            group_by(replicate, random.seed) %>%
            summarise(score = mean(score)) %>%
            pull(score)
        ),
      packages = "pdfCluster"
    )
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
      indrop.pca, indrop.spca.models[[4]], dims=1:50,
      umap_transform=matrix(c(0,-1,-1,0), nrow=2)
    ) %>%
      FindNeighbors(dims=1:34, red="spca", nn.method = "rann") %>%
      FindClusters(res=1.0, random.seed=1) %>%
      AddMetaData(
        Idents(.) %>% fct_recode(
          ISC="9", EB="11",
          # A massive cluster of SPC1-high cells.
          `EC-like`="0",
          EC.anterior1="1",
          EC.anterior2="2",
          EC.posterior1="3",
          LFC1="4",
          `copper/iron`="5",
          EE.Ast1="6",
          EC.anterior3="7",
          EC.posterior2="8",
          dEC1="10",
          EC.meso="12",
          EE.Ast2="13",
          others.1="14",
          `EC-like2`="15",
          cardia="16",
          dEC2="17",
          EC.anterior4="18",
          others.2="19",
          LFC2="20",
          EC.posterior3="21",
          `EC-like3`="22",
          others.3="23",
          `EC-like3`="24"
        ) %>% fct_relevel("ISC", "EB"),
        "spca_clusters"
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
    build_glms(indrop, midgut.pooledSizeFactors, c('pca_clusters', 'spca_clusters')),
    packages = "glmGamPoi"
  ),
  tar_target(
    indrop.present.genes,
    build_present_gene_list(
      indrop,
      c('pca_clusters', 'spca_clusters'),
      # rRNA does not usually show up when using many RNA-seq technologies and
      # we are not interested in quantifying it, so remove it from the GLM
      # analysis. Mitochondrial rRNA in Drosophila does show up using different
      # RNA-seq technologies, so to avoid removing genes that are actually
      # abundant, we will include those rRNAs in the analysis.
      select=matches("^mt:") | !matches("rRNA")
    )
  ),
  tar_target(
    indrop.deg,
    build_de_data(indrop.glm, indrop.present.genes),
    packages = c("apeglm", "glmGamPoi", "SummarizedExperiment")
  ),
  tar_target(
    indrop.misc,
    build_midgut_misc_stats(indrop, indrop.sct.pca)
  ),

  # Midgut 10X Genomics batch recapitulation
  tar_target(
    tenx.pca,
    # esg may be included as an spca gene as it helps separate ISC from dEC, and
    # also separates EB (esg hi, reason unknown) from ISC (esg mid).
    midgut_seurat_for_technology(
        midgut.counts, midgut.metadata, midgut.metafeatures, '10x',
        midgut.pooledSizeFactors
    )
  ),
  tar_target(
    tenx.pca.clusters,
    tenx.pca %>%
      FindNeighbors(dims = 1:30) %>%
      FindClusters(res = 1.5) %>%
      Idents %>%
      fct_recode(
        ISC='17', EB='18',
        pEC1='0',
        aEC1='1',
        dEC='2',
        aEC2='3',
        # We created a feature that separates EC-like from EC on the inDrop
        # batch: UMI_per_Feature <- nCount_RNA / nFeature_RNA
        # The inDrop EC-like UMI_per_feature was much closer to that of the
        # dEC feature, although the cells are much closer in terms of cosine
        # similarity to aEC cells. Using the new feature, we called cluster '4'
        # as corresponding to the largest inDrop 'EC-like' cluster.
        `EC-like1`='4',
        EE1='5',
        pEC2='6',
        `copper/iron`='7',
        aEC3='8',
        `copper/iron2`='9',
        LFC='10',
        mEC='11',
        EE2='12',
        pEC3='13',
        aEC4='14',
        pEC4='15',
        aEC5='16',
        # Again, UMI_per_Feature suggests that cluster '19' matches 'EC-like'
        # and not 'aEC'.
        `EC-like2`='19',
        `copper/iron3`='20',
        pEC5='21',
        EE3='22',

        # EbpIII+
        others='23',
        # Pgant4+
        cardia='24'
      ) %>%
      fct_relevel(c("ISC", "EB"))
  ),
  tar_target(
    tenx.spca.dimreduc,
    RunSparsePCA(
      tenx.pca, 'covar', varnum=12, npcs=50, eigen_gap=0.001, search_cap=500000,
      cgroup=memory_cgroups[2], do.correct.elbow=TRUE
    )[['spca']],
    # 'covar' matrix has not changed in the 'tenx.pca' object, so do not rerun
    # dimreduc target which has an ETA of 12-24 hours.
    cue=tar_cue('never')
  ),
  tar_target(
    tenx,
    spca_with_centered_umap(
      tenx.pca, tenx.spca.dimreduc, dims=1:50
    ) %>%
      FindNeighbors(dims = 1:22, red = "spca", nn.method = "rann") %>%
      FindClusters(res = 2, random.seed = 1) %>%
      AddMetaData(
        Idents(.) %>% fct_recode(
          ISC='19', EB='20',

          pEC1='0',
          aEC1='1',
          pEC2='2',
          EE1='3',
          `copper/iron1`='4',
          aEC2='5',
          `EC-like1`='6',
          pEC3='7',
          LFC='8',
          `EC-like2`='9',
          aEC3='10',
          EE2='11',
          `copper/iron2`='12',
          mEC='13',
          pEC4='14',
          dEC='15',
          pEC5='16',
          aEC4='17',
          `EC-like3`='18',
          `copper/iron3`='21',
          cardia='22',
          aEC5='23',
          others='24',
          aEC6='25'
        ) %>%
          fct_relevel(c("ISC", "EB")),
        "spca_clusters"
      ) %>%
      AddMetaData(tenx.pca.clusters, "pca_clusters") %>%
      midgut_classify_cell_types("pca_clusters") %>%
      midgut_classify_cell_types("spca_clusters")
  ),
  tar_target(
    tenx.glm,
    build_glms(tenx, midgut.pooledSizeFactors, c('pca_clusters', 'spca_clusters')),
    packages = "glmGamPoi"
  ),
  tar_target(
    tenx.present.genes,
    build_present_gene_list(
      tenx, c('pca_clusters', 'spca_clusters'),
      select=matches("^mt:") | !matches("rRNA")
    )
  ),
  tar_target(
    tenx.deg,
    build_de_data(tenx.glm, tenx.present.genes),
    packages = c("apeglm", "glmGamPoi", "SummarizedExperiment")
  ),
  tar_target(
    tenx.spca.param.k.input,
    tibble(k = 1:12)
  ),
  tar_target(
    tenx.spca.param.k,
    tenx.spca.param.k.input %>%
      rowwise %>%
      mutate(
        dimreduc = tenx.pca %>%
          seurat_spca("covar", varnum=k, npcs=1, eigen_gap=0.01, search_cap=500000) %>%
          list
      ),
    pattern = map(tenx.spca.param.k.input)
  ),
  tar_target(
    tenx.pca.quickCluster,
    tenx.pca %>%
      FindNeighbors(dims = 1:findPC(tenx.pca[['pca']]@stdev)) %>%
      FindClusters(res = 0.1) %>%
      Idents,
    packages = "findPC"
  ),
  tar_target(
    tenx.spca.modsel.features,
    rowAnys(
      tenx.spca.param.k %>%
        pull(dimreduc) %>%
        sapply(\(dimreduc) dimreduc@feature.loadings[, "SPARSE_1"]) %>%
        `!=`(0),
      useNames = TRUE
    ) %>%
      which %>%
      names
  ),

  # Midgut (Both inDrop and 10X) figures
  midgut_figures,
  midgut_figures_2,
  # tar_combine(midgut.figures, list(midgut_figures, midgut_figures_2)),

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
      acc.rna, 'covar', varnum=10, npcs=50, eigen_gap=0.001, search_cap=100000,
      cgroup=memory_cgroups[3], do.correct.elbow=TRUE
    )[['spca']],
    # 'covar' matrix has not changed in the 'acc.rna' object, so do not rerun
    # dimreduc target which has an ETA of 12-24 hours.
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
          HMM = TRUE,
          no_plot = TRUE
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
  tar_target(
    caf,
    acc %>% acc_to_caf_spca
  ),
  tar_target(
    acc_colData,
    acc_call_idents(acc, caf) %>%
      FetchData("ident") %>%
      as.data.frame %>%
      mutate(size_factor = pooledSizeFactors(acc[['RNA']]@counts, clusters = ident))
  ),
  tar_target(
    acc.present.genes,
    (acc[['RNA']]@counts != 0) %>%
      rowSums %>%
      `>=`(100) %>%
      which %>%
      names %>%
      union(pan_caf_genes$pCAF) %>%
      union(pan_caf_genes$parikh_pCAF)
  ),
  tar_target(
    acc.glm,
    acc_glm(acc[['RNA']]@counts[acc.present.genes,], acc_colData),
    packages = "glmGamPoi"
  ),
  tar_target(
    acc.spca.param.k.input,
    tibble(k = 1:12)
  ),
  tar_target(
    acc.spca.param.k,
    acc.spca.param.k.input %>%
      rowwise %>%
      mutate(
        dimreduc = acc %>%
          seurat_spca("covar", varnum=k, npcs=1, eigen_gap=0.01, search_cap=500000) %>%
          list
      ),
    pattern = map(acc.spca.param.k.input)
  ),
  tar_target(
    acc.pca.quickCluster,
    acc %>%
      FindNeighbors(dims = 1:findPC(acc[['pca']]@stdev)) %>%
      FindClusters(res = 0.1) %>%
      Idents,
    packages = "findPC"
  ),
  tar_target(
    acc.spca.modsel.features,
    rowAnys(
      acc.spca.param.k %>%
        pull(dimreduc) %>%
        sapply(\(dimreduc) dimreduc@feature.loadings[, "SPARSE_1"]) %>%
        `!=`(0),
      useNames = TRUE
    ) %>%
      which %>%
      names
  ),

  tar_target(
    pan_caf_genes,
    list(
      pCAF=c("CDC45", "CDC25C", "CDK1", "TOP2A", "BIRC5"),
      parikh_pCAF=c("PAX7", "MYF5"),
      nCAF=c("CXCR4", "TPD52", "TPD52L1", "APOC1"),
      iCAF2=c("TNFAIP3", "ICAM1", "CXCL2", "CLU", "BDKRB1"),
      iCAF=c("CXCL12", "CXCL14", "C3", "CFD"),
      dCAF=c("STC1", "MMP11", "MMP1", "COL10A1", "COL3A1", "COL1A1"),
      myCAF=c("MYLK", "MCAM", "TAGLN", "MYH11", "ACTA2")
    )
  ),
  acc_figures,
  # tar_combine(acc.figures, acc_figures),
  midgut_supplement_sce,
  score_models_supplement
)

# Supplemental analysis of midgut - dot plot data.
midgut_supplement_sce <- list(
  # Dot plot data generation for batch crossed with PCA/SPCA workflow.
  tar_map(
    tribble(
      ~model, ~batch, ~seurat, ~glm_fits,
      "pca_subclassif", "indrop", rlang::sym("indrop"), rlang::sym("indrop.glm"),
      "spca_subclassif", "indrop", rlang::sym("indrop"), rlang::sym("indrop.glm"),
      "pca_subclassif", "tenx", rlang::sym("tenx"), rlang::sym("tenx.glm"),
      "spca_subclassif", "tenx", rlang::sym("tenx"), rlang::sym("tenx.glm")
    ) %>%
      rowwise %>%
      mutate(
        glm_model = model %>%
          str_replace("pca_subclassif", "pca_clusters") %>%
          rlang::sym() %>%
          call("$", glm_fits, .) %>%
          list
      ),
    names = c("model", "batch"),
    tar_target(
      dispersion.trend,
      with(
        glm_model,
        overdispersion_shrinkage_list$dispersion_trend %>%
          setNames(names(overdispersions))
      )
    ),
    tar_target(
      glm.abundances,
      build_glm_cpm(seurat, model, dispersion.trend),
      packages = "glmGamPoi"
    ),
    tar_target(
      cpm,
      predict_glm_cpm(glm.abundances)
    ),
    tar_target(
      pct.expressed,
      seurat[['RNA']] %>%
        GetAssayData %>%
        t %>%
        as.data.frame %>%
        split(FetchData(seurat, model)[, 1]) %>%
        sapply(\(df) colMeans(as.matrix(df) != 0))
    )
  ),
  # Figures that are plotting expression of the following genes:
  tar_target(
    midgut.dot.plot.features,
    c(
      "esg", "Dl", "E(spl)mbeta-HLH", "Rel",
      # dEC-enriched:
      "fmt", "trol",
      # aEC markers (also EC-like):
      "betaTry",
      "Amy-p",
      "Vha16-1",
      "Gs2",
      "MtnC",
      "mag",
      "Pgant4",
      "pros",
      "EbpIII"
    )
  ),
  tar_target(
    fig.indrop.dots.pca,
    save_figure(
      "figure/Midgut/Dot-Plot-Indrop-PCA.pdf",
      midgut_dot_plot(
        cpm_pca_subclassif_indrop[midgut.dot.plot.features, ],
        pct.expressed_pca_subclassif_indrop[midgut.dot.plot.features, ],
        midgut.model.colors.bg["PCA"]
      ),
      width = 6.4, height = 4.8, device = cairo_pdf
    )
  ),
  tar_target(
    fig.indrop.dots.spca,
    save_figure(
      "figure/Midgut/Dot-Plot-Indrop-SPCA.pdf",
      midgut_dot_plot(
        cpm_spca_subclassif_indrop[midgut.dot.plot.features, ],
        pct.expressed_spca_subclassif_indrop[midgut.dot.plot.features, ],
        midgut.model.colors.bg["SPCA"]
      ),
      width = 6.4, height = 4.8, device = cairo_pdf
    )
  )
)
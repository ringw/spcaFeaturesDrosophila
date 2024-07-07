midgut_figures_2 <- list(
  tar_target(
    fig.indrop.violin.features,
    save_figure(
      "figure/Midgut/Fig3-Violin.pdf",
      ggarrange(
        tiny_violin_plot(
          indrop, "pca_classif", c("EB", "ISC", "EC", "EC-like", "EE"), "PC_3"
        ) + labs(
          tag = "A"
        ),
        tiny_violin_plot(
          indrop, "spca_classif", c("EB", "ISC", "EC", "EC-like", "EE"), "SPARSE_1"
        ) + labs(
          tag = "B"
        ),
        plot_spca_cdf(
          indrop, "spca_classif", c("EB", "ISC", "EC", "EC-like", "EE"), "SPARSE_1"
        ) + labs(
          tag = "C"
        ),
        tiny_violin_plot(
          indrop, "pca_classif", c("EB", "ISC", "EC", "EC-like", "EE"), "PC_7"
        ) + labs(
          tag = "D"
        ),
        tiny_violin_plot(
          indrop, "spca_classif", c("EB", "ISC", "EC", "EC-like", "EE"), "SPARSE_26"
        ) + labs(
          tag = "E"
        ),
        plot_spca_cdf(
          indrop, "spca_classif", c("EB", "ISC", "EC", "EC-like", "EE"), "SPARSE_26"
        ) + labs(
          tag = "F"
        )
      ),
      12,
      6
    )
  ),
  tar_target(
    indrop.expl.var.stats.logumi,
    midgut_expl_var_stem_like(indrop, indrop.sct.pca)
  ),
  tar_target(
    indrop.expl.var.stemlike.inputs,
    c("pca", "spca2pca", "pca.sct")
  ),
  tar_file(
    fig.indrop.expl.var.stemlike,
    save_figure(
      paste0("figure/Midgut/ExplVar-", indrop.expl.var.stemlike.inputs, ".pdf"),
      plot_midgut_expl_var(indrop.expl.var.stats.logumi, indrop.expl.var.stemlike.inputs),
      4,
      1.25
    ),
    pattern = map(indrop.expl.var.stemlike.inputs)
  ),
  tar_file(
    fig.indrop.expl.var.stacked,
    save_figure(
      "figure/Midgut/ExplVar-ISC-EB.pdf",
      midgut_expl_var_stem_between(indrop),
      3, 3
    ),
    packages = c(tar_option_get("packages"), "ggpattern")
  ),

  tar_target(
    indrop.spca.dimreduc.sklearn,
    seurat_spca_from_feature_loadings(
      seurat_joint_spca(
        indrop.pca,
        alpha = 0.265,
        ridge_alpha = 0.01,
        n_components = 50
      ) %>%
        t,
      indrop.pca,
      "RNA",
      do.correct.elbow = TRUE
    )
  ),
  tar_target(
    indrop.spca.dimreduc.elasticnet,
    seurat_spca_from_feature_loadings(
      seurat_elasticnet_spca(
        indrop.pca,
        varnum = 8,
        n_components = 50
      ),
      indrop.pca,
      "RNA",
      do.correct.elbow = TRUE,
      do.rename.features = FALSE
    ),
    cue = tar_cue("never")
  ),
  tar_file(
    fig.indrop.feature.loadings,
    tribble(
      ~name, ~dimreduc, ~k,
      "SKLearn", indrop.spca.dimreduc.sklearn, 20,
      "LASSO", indrop.spca.dimreduc.elasticnet, 8,
      "OptimalSPCA", indrop.spca.models[[4]], 8
    ) %>%
      rowwise %>%
      mutate(
        gg = make_dimreduc_img(dimreduc, k=k) %>% list,
        filename = str_glue("figure/Midgut/Feature-Loadings-{name}.pdf"),
        ggsave = ggsave(filename, gg, width=4, height=6)
      ) %>%
      pull(filename)
  ),
  tar_map(
    tibble(feature = c("SPARSE_1", "SPARSE_26")) %>%
      full_join(
        tribble(~font_size, ~suffix_font_size, 4, "", 8, "-Thumbnail"),
        character(0)
      ),
    names = feature | font_size,
    tar_file(
      fig.indrop.loadings,
      save_figure(
        paste0("figure/Midgut/Feature-Loadings-", feature, suffix_font_size, ".pdf"),
        (
          spc_tile_plot(indrop.spca.models[[4]], feature, fontsize=font_size)
          + guides(fill = guide_none())
        ),
        width = if (font_size > 4) 1.75 else 1.2,
        height = 3.2
      ),
      packages = tar_option_get("packages") %>% c("Cairo")
    )
  ),
  tar_file(
    fig.indrop.loadings.thumbnail.legend,
    save_figure(
      "figure/Midgut/Feature-Loadings-Thumbnail-Legend.pdf",
      get_legend(
        ggplot(data.frame(x=c(0, 0.5)), aes(x, y=0, color=x))
        + geom_tile()
        + scale_color_viridis_c(
          breaks = c(0, 0.25, 0.5),
          labels = percent, end = 0.5,
          guide = guide_colorbar(title = NULL, barwidth = 2.5, barheight = 0.75)
        )
        + theme(legend.position = "bottom")
      ),
      width = 1.25,
      height = 0.5
    )
  )
)
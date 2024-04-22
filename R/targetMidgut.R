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
  )
)
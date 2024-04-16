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
  )
)
midgut_figures_2 <- list(
  tar_target(
    fig.indrop.umap,
    save_figure(
      "figure/Midgut/Fig2-Features.pdf",
      plot_indrop_fig2(indrop),
      6,
      4.75
    ),
    format = "file",
    packages = tar_option_get("packages") %>% c("egg")
  ),
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
      plot_midgut_expl_var(indrop.expl.var.stats.logumi, indrop.expl.var.stemlike.inputs)
      + theme(axis.text = element_text(size = rel(1.5))),
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
          spc_tile_plot(
            indrop.spca.models[[4]], feature, fontsize=font_size,
            fontface = c(`4`="bold", `8`="plain")[as.character(font_size)]
          )
          + guides(fill = guide_none())
        ),
        width = if (font_size > 4) 1.6 else 1.2,
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
        + geom_tile() +
        scale_color_gradientn(
          limits=c(-0.5, 0.5000001),
          colors = c(
            seq_gradient_pal("#e0524d", viridis(10)[1])(
              seq(0, 1, length.out=50)[-50]
            ),
            viridis(101)[1:51]
          ),
          guide = guide_colorbar(title = NULL, barwidth = 7, barheight = 0.75),
          labels = percent
        )
        + theme(legend.position = "bottom")
      ),
      width = 1.75,
      height = 0.5
    )
  ),
  tar_file(
    fig.indrop.loadings.pc3,
    save_figure(
      "figure/Midgut/Feature-Loadings-PC3-Thumbnail.pdf",
      (
        pc_tile_plot(indrop[["pca"]], "PC_3", fontsize=8, fontface="plain", begin=-0.5)
        + guides(fill = guide_none())
      ),
      width = 1.6,
      height = 3.6
    ),
    packages = tar_option_get("packages") %>% c("Cairo")
  ),
  tar_file(
    fig.indrop.loadings.pc7,
    save_figure(
      "figure/Midgut/Feature-Loadings-PC7-Thumbnail.pdf",
      (
        pc_tile_plot(indrop[["pca"]], "PC_7", fontsize=8, fontface="plain", begin=-0.5)
        + guides(fill = guide_none())
      ),
      width = 1.6,
      height = 3.6
    ),
    packages = tar_option_get("packages") %>% c("Cairo")
  ),
  tar_target(
    fig.indrop.thumbnail.density.pc3,
    save_figure(
      "figure/Midgut/Indrop-Thumbnail-Density-Schematic-PC3.pdf",
      tiny_density_plot(
        indrop,
        "pca_classif",
        "EC",
        "PC_3"
      ),
      width=3, height=1.5
    ),
    format = "file"
  ),
  tar_target(
    fig.indrop.thumbnail.density.spc26,
    save_figure(
      "figure/Midgut/Indrop-Thumbnail-Density-Schematic-SPC26.pdf",
      tiny_density_plot(
        indrop,
        "spca_classif",
        c("ISC", "EB"),
        "SPARSE_26"
      ),
      width=3, height=1.5
    ),
    format = "file"
  ),
  tar_target(
    indrop.heatmap.genes,
    c(
      # "betaTry", "CG30025", "CG12374", "alphaTry", "Bace",
      # "yip7",
      # "Amy-p", "Mal-A1",
      # "gammaTry", "CG30031", "deltaTry", "Jon65Aiv", "Jon65Aiii",
      # "Hsp23", "IA-2"
      "gammaTry", "CG30031", "deltaTry", "CG30025",
      "yip7",
      "Jon99Ciii", "Jon99Cii", "Jon65Aiv", "Jon65Aiii",
      "Pebp1",
      "N",
      "Hsp23"
    )
  ),
  tar_target(
    indrop.heatmap.colorscale,
    tibble(
      x = seq(-0.25, 1, by=0.05),
      y = c(seq_gradient_pal("#e0524d", viridis(10)[1])(seq(0, 1, length.out=6))[-6], viridis(21))
    )
  ),
  tar_target(
    indrop.heatmaps,
    {
      Heatmap <- indrop@misc$covar[indrop.heatmap.genes, indrop.heatmap.genes]
      eigs <- eigen(Heatmap, sym=TRUE)
      PC1 <- eigs$values[1] * tcrossprod(eigs$vectors[, 1])
      PC2 <- eigs$values[2] * tcrossprod(eigs$vectors[, 2])
      eigs_SPC1 <- eigen(
        diag(c(rep(1, 4), rep(0, 8)))
        %*%
        Heatmap
        %*%
        diag(c(rep(1, 4), rep(0, 8))),
        sym=T
      )
      SPC1 <- eigs_SPC1$values[1] * tcrossprod(eigs_SPC1$vectors[, 1])
      eigs_SPC2 <- eigen(
        diag(c(rep(0, 5), rep(1, 4), rep(0, 3)))
        %*%
        Heatmap
        %*%
        diag(c(rep(0, 5), rep(1, 4), rep(0, 3)))
        ,
        sym=T
      )
      SPC2 <- eigs_SPC2$values[1] * tcrossprod(eigs_SPC2$vectors[, 1])
      list(
        Heatmap=Heatmap,
        PC1=matrix(PC1, nrow=nrow(PC1), dimnames=dimnames(Heatmap)),
        PC2=matrix(PC2, nrow=nrow(PC2), dimnames=dimnames(Heatmap)),
        SPCA=matrix((SPC1 + SPC2) %>% replace(. == 0, NA), nrow=nrow(SPC1), dimnames=dimnames(Heatmap))
      ) %>%
        enframe
    }
  ),
  tar_file(
    fig.indrop.heatmap,
    save_figure(
      str_glue("figure/Midgut/Covar-Heatmap-{indrop.heatmaps$name}.pdf"),
      melt(indrop.heatmaps$value[[1]]) %>%
        ggplot(aes(Var1, Var2, fill=value, color=value))
        + geom_tile()
        + scale_x_discrete(name=NULL, breaks=NULL)
        + scale_y_discrete(name=NULL, breaks=NULL, limits=rev)
        + scale_fill_gradientn(
          limits=range(indrop.heatmap.colorscale$x)+c(0,1e-10),
          colors=indrop.heatmap.colorscale$y,
          na.value="transparent"
        )
        + scale_color_gradientn(
          limits=range(indrop.heatmap.colorscale$x)+c(0,1e-10),
          colors=indrop.heatmap.colorscale$y,
          na.value="transparent"
        )
        + annotate(
          "rect",
          xmin=0.5,
          xmax=12.5,
          ymin=0.5,
          ymax=12.5,
          color="black",
          fill="transparent"
        )
        + (
          if (indrop.heatmaps$name == "SPCA")
            annotate(
              "rect",
              # First arg: Outline the SPC1. Second arg: Outline the SPC2.
              xmin=c(0.5, 5.5),
              xmax=c(4.5, 9.5),
              ymin=c(12.5, 7.5),
              ymax=c(8.5, 3.5),
              color="black",
              fill="transparent"
            )
        )
        + coord_cartesian(expand=F)
        + theme(
          panel.background = element_rect(fill="transparent"), aspect.ratio = 1,
          legend.position = "none"
        ),
      width=2,
      height=2
    ),
    pattern = map(indrop.heatmaps)
  ),
  tar_file(
    fig.indrop.heatmap.legend,
    save_figure(
      "figure/Midgut/Covar-Legend.pdf",
      get_legend(
        ggplot(tibble(x=0, y=0, value=0.5), aes(x, y, color=value))
        + geom_point()
        + scale_color_gradientn(
          name = "Loading",
          limits=range(indrop.heatmap.colorscale$x)+c(0,1e-10),
          colors=indrop.heatmap.colorscale$y,
          na.value="transparent"
        )
      ),
      width=0.75,
      height=1.75
    )
  )
)

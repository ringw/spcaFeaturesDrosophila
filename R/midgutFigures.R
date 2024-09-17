midgut.colors = c(
  'ISC'="#FAD11B", # hcl(68, 95, 85)
  'EB'="#FF7349", # hcl(21, 120, 65)
  'dEC'="#FF82BA", # hcl(345, 80, 70)
  'EC'="#34CDC7", # hcl(187, 55, 75)
  'EC-like'="#A0E6EC", # hcl(10, 48, 74)
  'EE'="#BCA6FF", # hcl(274, 71, 73)
  'copper/iron'="#c58b81",
  'LFC'="#48F183", # hcl(136, 95, 85)
  'cardia'="#891416", # hcl(12, 84, 29)
  'bg'=hcl(c = 0, l = 87)
)
midgut.col = as.list(midgut.colors)

midgut.model.colors.bg = c(PCA=hcl(30, 8, 99), SPCA=hcl(129,6,99))
midgut.model.colors.legend = c(PCA=hcl(30, 12, 95), SPCA=hcl(129,11,95))

# Outlier cells in the -UMAP_2 direction in PCA-UMAP:
# > sum(indrop[['umap']]@cell.embeddings[, "UMAP_2"] < -9.5)
# 24
# > sum(indrop[['umap.spca']]@cell.embeddings[, "umapspca_2"] < -9.5)
# 0

plot_indrop_pca_annotations <- tribble(
  ~UMAP_1, ~UMAP_2, ~label,
  -10.95, 2.45, "aEC",
  -0.5, -8.25, "aEC",
  9.5, 7.25, "mEC",
  #Gs2+
  -4.35, 0, "pEC",
  # Above/left: LManVI+. Below: Mur29B+
  7.5, 0.6, "pEC"
)

plot_indrop_pca <- function(indrop, dpi=300, text_size=6.5) (
  indrop@meta.data %>% cbind(as.data.frame(indrop[['umap']]@cell.embeddings)) %>%
  arrange(pca_classif != 'unknown') %>%
  ggplot(aes(UMAP_1,UMAP_2, color=pca_classif))
  + rasterise(geom_point(shape=20, size=1e-3, show.legend = F), dpi=dpi)
  + theme_bw()
  + scale_color_manual(values=setNames(midgut.colors,NULL), guide=guide_legend(override.aes=list(size=3)))
  + scale_x_continuous(limits=c(-13,NA), expand=rep(0.02,2), breaks=pretty_breaks(4))
  + scale_y_continuous(limits=c(-9.22658,NA), expand=rep(0.02,2), breaks=pretty_breaks(4))
  + geom_text(
    aes(label=label), data=plot_indrop_pca_annotations, color="black", size=text_size
  )
  + labs(x=bquote("UMAP"[1]), y=bquote("UMAP"[2]))
  + theme(
    axis.ticks = element_line(color='transparent'),
    panel.background = element_rect(fill=midgut.model.colors.bg[1]),
    aspect.ratio = 3/4,

    # Double the relative size of text as we created an over-sized pdf.
    axis.title = element_text(size = rel(1.6)),
    axis.title.x = element_text(margin = margin(1, 0, 0, 0)),
    axis.title.y = element_text(margin = margin(0, -7, 0, 0)),
    axis.text = element_text(size = rel(1.6)),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1)
  )
)

plot_indrop_spca_annotations <- tribble(
  ~umapspca_1, ~umapspca_2, ~label,
  -2, -6.75, "aEC",
  # Amy-p+ cells
  -1.7, 7.8, "aEC",
  6.5, 6.5, "mEC",
  1.5, 6.35, "pEC",
  6.75, 1.75, "pEC"
)

plot_indrop_spca <- function(indrop, dpi=300, text_size=6.5) (
  indrop@meta.data %>%
    cbind(as.data.frame(indrop[['umap.spca']]@cell.embeddings)) %>%
    arrange(pca_classif != 'unknown') %>%
    ggplot(aes(umapspca_1,umapspca_2, color=spca_classif))
  + rasterise(geom_point(shape=20, size=1e-3, show.legend = F), dpi=dpi)
  + theme_bw()
  + scale_color_manual(values=setNames(midgut.colors,NULL), guide=guide_legend(override.aes=list(size=3)))
  + scale_x_continuous(expand=rep(0.02,2), breaks=pretty_breaks(4))
  + scale_y_continuous(limits=c(-7.547306,NA), expand=rep(0.02,2), breaks=pretty_breaks(4))
  + geom_text(
    aes(label=label), data=plot_indrop_spca_annotations, color="black", size=text_size
  )
  + labs(x=bquote("UMAP"[1]), y=bquote("UMAP"[2]))
  + theme(
    axis.ticks = element_line(color='transparent'),
    panel.background = element_rect(fill=midgut.model.colors.bg[2]),
    aspect.ratio = 3/4,

    # Double the relative size of text as we created an over-sized pdf.
    axis.title = element_text(size = rel(1.6)),
    axis.title.x = element_text(margin = margin(1, 0, 0, 0)),
    axis.title.y = element_text(margin = margin(0, -7, 0, 0)),
    axis.text = element_text(size = rel(1.6))
  )
)

plot_midgut_legend <- function(indrop) get_legend(
  indrop@meta.data %>% cbind(as.data.frame(indrop[['umap']]@cell.embeddings))
  %>% arrange(pca_classif != 'unknown')
  %>% ggplot(aes(UMAP_1,UMAP_2, color=pca_classif))
  + geom_point(size=2)
  + scale_color_manual(values=midgut.colors %>% setNames(names(.) %>% str_replace('bg','unknown')))
  + theme_bw() + theme(legend.direction='horizontal')
  + labs(color=NULL)
)

plot_midgut_feature <- function(indrop, bg_color, embedding, feature_name, limits=NULL, dpi=300, pt.size=0.25) (
  indrop[[embedding]]@cell.embeddings %*%
    matrix(diag(2), nrow=2, dimnames=list(NULL, c("UMAP_1", "UMAP_2"))) %>%
    as.data.frame %>%
    cbind(indrop@meta.data) %>%
    cbind(
      FetchData(indrop, feature_name) %>%
        pull(1) %>%
        matrix(ncol = 1, dimnames=list(names(.), "feature")) %>%
        `*`(
          if (grepl("SPARSE_", feature_name))
            1/sum(unlist(indrop[['spca']][, feature_name]))
          else 1
        )
    ) %>%
  ggplot(aes(UMAP_1, UMAP_2, color=feature))
  + rasterise(geom_point(stroke=NA, size=pt.size, show.legend = F), dpi=dpi)
  + theme_bw()
  + scale_x_continuous(
    limits = \(ll) c(pmax(ll[1], -13), ll[2]),
    expand=rep(0.02,2), breaks=pretty_breaks(4)
  )
  + scale_y_continuous(
    limits = \(ll) c(quantile(indrop[[embedding]]@cell.embeddings[,2], 0.005), ll[2]),
    expand=rep(0.02,2), breaks=pretty_breaks(4)
  )
  + scale_color_viridis_c(
    option='magma', end=0.9, limits=limits,
    oob=squish
  )
  + labs(x=bquote("UMAP"[1]), y=bquote("UMAP"[2]))
  + theme(
    axis.ticks = element_line(color='transparent'),
    panel.background = element_rect(fill=midgut.model.colors.bg[bg_color]),
    aspect.ratio = 3/4,

    # Double the relative size of text as we created an over-sized pdf.
    axis.title = element_text(size = rel(1.6)),
    axis.title.x = element_text(margin = margin(1, 0, 0, 0)),
    axis.title.y = element_text(margin = margin(0, -7, 0, 0)),
    axis.text = element_text(size = rel(1.6))
  )
)

plot_midgut_feature_legend <- function(indrop, feature_name="betaTry", limits=NULL, legend.direction="horizontal", legend.name="LogNormalize", guide=NULL) get_legend(
  plot_midgut_feature(indrop, "PCA", "umap", feature_name, limits) +
    geom_point() +
    guides(color = guide) +
    labs(color=legend.name) +
    theme_bw() +
    theme(legend.direction = legend.direction)
)

plot_midgut_model_background_legend <- function() get_legend(
  data.frame(
    name=c('PCA','SPCA')
  )
  %>%
  ggplot(aes(0,0,fill=name,color=name))
  # R shape 15 = square
  # + geom_point(size=4, shape=15)
  + geom_tile()
  + scale_fill_manual(values=midgut.model.colors.legend)
  + scale_color_manual(values=rep('black', 2), guide=guide_legend(override.aes = list()))
  + labs(color='model', fill='model')
  + theme_bw()
  + theme(legend.direction='horizontal')
)

midgut_dot_plot <- function(cpm_data, pct_data, bg_color=waiver()) {
  names(dimnames(cpm_data)) = c("gene", "cluster")
  names(dimnames(pct_data)) = c("gene", "cluster")
  scale_cpm_data <- cpm_data %>% t %>% scale %>% t
  plot_data <- melt(scale_cpm_data, value.name="scaleCPM") %>%
    inner_join(
      melt(pct_data, value.name="pct"),
      c("gene", "cluster")
    )
  plot_data$gene <- plot_data$gene %>% factor(rownames(cpm_data)) %>%
    display_gene_names
  plot_data <- plot_data %>% mutate(
    is_faint = ifelse(between(scaleCPM, 1, 1.75), "present", "absent"),
    is_faint_or_high = ifelse(scaleCPM > 1, "present", "absent")
  )

  ggplot(plot_data, aes(gene, cluster, color=is_faint_or_high, fill=scaleCPM, size=pct)) + geom_point(
    pch=21, stroke=0.1
  ) + scale_fill_distiller(
    type = "div", palette = "RdYlBu",
    limits = c(-0.5, 2.99), oob=squish,
    guide = guide_colorbar(title = "scale(CPM)", barheight = 3, barwidth = 1)
  ) + scale_color_manual(
    # Is faint (very light yellow) color or a higher color: Show a very dim stroke around
    # the point.
    values = c(absent="transparent", present="#33333333"),
    guide = guide_none()
  ) + scale_y_discrete(
    limits = rev
  ) + scale_size(
    labels = percent,
    range = c(0.1, 4),
    guide = guide_legend(override.aes = list(fill = "black"))
  ) + theme_bw() + theme(
    panel.background = element_rect(fill = bg_color),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(margin = margin(r = 2, t = 5.5, b = 5.5, l = -5)),
    legend.margin = margin(t = 20)
  ) + labs(
    x = NULL, y = NULL, size = "% Expressed"
  )
}

# Fig2 gtable plot. Needs to show the PCA-UMAP and SPCA-UMAP clusters,
# annotations, then plots of the E(spl)mb and SPC26 features annotated with
# inset box, then the inset panels.
plot_indrop_fig2 <- function(
  indrop
) {
  pt <- \(x) x * 25.4 / 72
  main_width <- unit(1.6, "in")
  main_height <- unit(1.2, "in")
  inset_width <- unit(1, "in")
  inset_height <- unit(0.75, "in")
  axis_theme <- theme(
    axis.text = element_text(size = unit(8, "pt")),
    axis.text.x = element_text(margin = margin(t = 0, b = 1)),
    axis.title = element_text(size = unit(8, "pt")),
    axis.title.x = element_text(margin = margin(1, 0, 0, 0)),
    axis.title.y = element_text(margin = margin(0, -4, 0, 0)),
    panel.border = element_rect(fill = NA, size = unit(0.25, "pt")),
    plot.margin = margin(t = 1, r = 5.5, b = 5, l = 5.5)
  )
  inset_theme <- axis_theme +
    theme(
      axis.title.y = element_text(margin = margin()),
      plot.margin = margin(t = -5, r = 1, b = 5, l = 1)
    )
  blank_ggplot <- (
    ggplot() +
      axis_theme +
      theme(panel.background = element_blank(), panel.border = element_blank())
  ) %>%
    set_panel_size(w = main_width, h = main_height)
  panel_limits <- list(
    pca = list(x = c(-13.35423, 12.19315), y = c(14.306092, -9.707613)),
    spca = list(x = c(-10.29552, 8.25750), y = c(13.074905, -7.551983))
  )
  annotate_label <- sapply(
    c("E(spl)mbeta-HLH", "SPC26"),
    \(n) sapply(
      panel_limits,
      \(limits) annotate(
        "text",
        limits$x[1] + 0.02 * diff(limits$x),
        limits$y[1] + 0.07 * diff(limits$y),
        label = n,
        hjust = 0,
        size = pt(8)
      )
    ),
    simplify = FALSE
  )
  make_inset <- \(gr) {
    container <- blank_ggplot
    container$grobs[[
      match("panel", container$layout$name)
    ]] <- gr
    container <- container
  }
  main_grid <- rbind(
    cbind(
      plot_indrop_pca(indrop, 1200, pt(8)) %>%
        `+`(axis_theme) %>%
        set_panel_size(w = main_width, h = main_height),
      plot_indrop_spca(indrop, 1200, pt(8)) %>%
        `+`(axis_theme) %>%
        set_panel_size(w = main_width, h = main_height),
      blank_ggplot
    ),
    cbind(
      plot_midgut_feature(
        indrop,
        bg_color="PCA",
        "umap",
        "E(spl)mbeta-HLH",
        limits=c(0, 4.5)
      ) %>%
        `+`(
          annotate(
            "rect", xmin=-3.8, ymin=7.7, xmax=-0.2, ymax=2.5, fill="transparent", color="black", linewidth=0.5
          )
        ) %>%
        `+`(annotate_label$`E(spl)mbeta-HLH`$pca) %>%
        `+`(axis_theme) %>%
        `+`(
          theme(
            plot.margin = axis_theme$plot.margin + margin(t = 5)
          )
        ) %>%
        set_panel_size(w = main_width, h = main_height),
      plot_midgut_feature(
        indrop,
        bg_color="SPCA",
        "umap.spca",
        "E(spl)mbeta-HLH",
        limits=c(0, 4.5)
      ) %>%
        `+`(
          annotate(
            "rect",xmin=-6.5, ymin=4.9, xmax=-2.5, ymax=2.3, fill="transparent", color="black", linewidth=0.5
          )
        ) %>%
        `+`(annotate_label$`E(spl)mbeta-HLH`$spca) %>%
        `+`(axis_theme) %>%
        set_panel_size(w = main_width, h = main_height),
      plot_midgut_feature(
        indrop,
        bg_color="SPCA",
        "umap.spca",
        "SPARSE_26",
        limits=c(0, 4.5)
      ) %>%
        `+`(
          annotate(
            "rect",xmin=-6.5, ymin=4.9, xmax=-2.5, ymax=2.3, fill="transparent", color="black", linewidth=0.5
          )
        ) %>%
        `+`(annotate_label$SPC26$spca) %>%
        `+`(axis_theme) %>%
        set_panel_size(w = main_width, h = main_height)
    ),
    cbind(
      plot_midgut_feature(
        indrop,
        bg_color="PCA",
        "umap",
        "E(spl)mbeta-HLH",
        limits=c(0, 4.5),
        pt.size = 0.5
      ) %>%
        `+`(
          coord_cartesian(c(-3.8, -0.2), c(2.5, 7.7), expand=F)
        ) %>%
        `+`(inset_theme) %>%
        set_panel_size(w = inset_width, h = inset_height) %>%
        make_inset(),
      plot_midgut_feature(
        indrop,
        bg_color="SPCA",
        "umap.spca",
        "E(spl)mbeta-HLH",
        limits=c(0, 4.5),
        pt.size = 0.5
      ) %>%
        `+`(
          coord_cartesian(c(-6.5, -2.5), c(2.3, 4.9), expand=F)
        ) %>%
        `+`(
          scale_y_continuous(breaks = c(3, 4))
        ) %>%
        `+`(inset_theme) %>%
        set_panel_size(w = inset_width, h = inset_height) %>%
        make_inset(),
      plot_midgut_feature(
        indrop,
        bg_color="SPCA",
        "umap.spca",
        "SPARSE_26",
        limits=c(0, 4.5),
        pt.size = 0.5
      ) %>%
        `+`(
          coord_cartesian(c(-6.5, -2.5), c(2.3, 4.9), expand=F)
        ) %>%
        `+`(
          scale_y_continuous(breaks = c(3, 4))
        ) %>%
        `+`(inset_theme) %>%
        set_panel_size(w = inset_width, h = inset_height) %>%
        make_inset()
    )
  )
}

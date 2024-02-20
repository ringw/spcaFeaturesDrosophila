midgut.colors = c(
  'ISC'=hsv(0.88, 0.75, 0.97),
  'EB'=hsv(0.166, 1, 0.95),
  'dEC'=hcl(137, 57, 71),
  'EC'=hcl(253, 74, 48),
  'EC-like'=hcl(230, 34, 72),
  'EE'=hcl(300, 72, 46),
  'copper/iron'=hcl(48, 63, 75),
  'LFC'=hcl(180, 38, 88),
  'cardia'=hcl(9, 66, 32),
  'bg'=hcl(c = 0, l = 87)
)
midgut.col = as.list(midgut.colors)

midgut.model.colors.bg = c(PCA=hcl(30, 8, 99), SPCA=hcl(129,7,98))
midgut.model.colors.legend = c(PCA=hcl(30, 12, 95), SPCA=hcl(129,11,95))

# Outlier cells in the -UMAP_2 direction in PCA-UMAP:
# > sum(indrop[['umap']]@cell.embeddings[, "UMAP_2"] < -9.5)
# 24
# > sum(indrop[['umap.spca']]@cell.embeddings[, "umapspca_2"] < -9.5)
# 0

plot_indrop_pca_annotations <- tribble(
  ~UMAP_1, ~UMAP_2, ~label,
  -8.75, 1.1, "aEC",
  0.8, -8.25, "aEC",
  7.75, 0.5, "aEC",
  11, 2, "mEC",
  4, -0.6, "pEC",
  0.75, 5.5, "pEC"
)

plot_indrop_pca <- function(indrop) (
  indrop@meta.data %>% cbind(as.data.frame(indrop[['umap']]@cell.embeddings)) %>%
  arrange(pca_classif != 'unknown') %>%
  ggplot(aes(UMAP_1,UMAP_2, color=pca_classif))
  + rasterise(geom_point(shape=20, size=1e-3, show.legend = F), dpi=300)
  + theme_bw()
  + scale_color_manual(values=setNames(midgut.colors,NULL), guide=guide_legend(override.aes=list(size=3)))
  + scale_y_continuous(limits=c(-9.5,NA), expand=rep(0.02,2), breaks=pretty_breaks(3))
  + scale_x_continuous(expand=rep(0.02,2), breaks=pretty_breaks(4))
  + geom_text(
    aes(label=label), data=plot_indrop_pca_annotations, color="black", size=2.5
  )
  + labs(x=bquote("UMAP"[1]), y=bquote("UMAP"[2]))
  + theme(
    axis.ticks = element_line(color='transparent'),
    panel.background = element_rect(fill=midgut.model.colors.bg[1])
  )
)

plot_indrop_spca_annotations <- tribble(
  ~umapspca_1, ~umapspca_2, ~label,
  -5.5, 1, "aEC",
  0, -7.5, "aEC",
  4.25, 9, "mEC",
  2.55, 2.05, "pEC"
)

plot_indrop_spca_connectors <- tribble(
  ~x, ~y, ~xend, ~yend, 2.4, 2.25, 1, 2.8, 2.8, 1.95, 6, 1
)

plot_indrop_spca <- function(indrop) (
  indrop@meta.data %>%
    cbind(as.data.frame(indrop[['umap.spca']]@cell.embeddings)) %>%
    arrange(pca_classif != 'unknown') %>%
    ggplot(aes(umapspca_1,umapspca_2, color=spca_classif))
  + rasterise(geom_point(shape=20, size=1e-3, show.legend = F), dpi=300)
  + theme_bw()
  + scale_color_manual(values=setNames(midgut.colors,NULL), guide=guide_legend(override.aes=list(size=3)))
  + scale_y_continuous(limits=c(-9.5,NA), expand=rep(0.02,2), breaks=pretty_breaks(3))
  + scale_x_continuous(expand=rep(0.02,2), breaks=pretty_breaks(4))
  + geom_text(
    aes(label=label), data=plot_indrop_spca_annotations, color="black", size=2.5
  )
  + geom_segment(
    aes(x=x, y=y, xend=xend, yend=yend),
    data=plot_indrop_spca_connectors, color="black", alpha=0.3
  )
  + labs(x=bquote("UMAP"[1]), y=bquote("UMAP"[2]))
  + theme(
    axis.ticks = element_line(color='transparent'),
    panel.background = element_rect(fill=midgut.model.colors.bg[2])
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

plot_midgut_feature <- function(indrop, bg_color, embedding, feature_name) (
  indrop[[embedding]]@cell.embeddings %*%
    matrix(diag(2), nrow=2, dimnames=list(NULL, c("UMAP_1", "UMAP_2"))) %>%
    as.data.frame %>%
    cbind(indrop@meta.data) %>%
    cbind(FetchData(indrop, feature_name) %>% pull(1) %>% matrix(ncol = 1, dimnames=list(names(.), "feature"))) %>%
  ggplot(aes(UMAP_1, UMAP_2, color=feature))
  + rasterise(geom_point(shape=20, size=1e-3, show.legend = F), dpi=300)
  + theme_bw()
  + scale_y_continuous(limits=c(-9.5,NA), expand=rep(0.02,2), breaks=pretty_breaks(3))
  + scale_x_continuous(expand=rep(0.02,2), breaks=pretty_breaks(4))
  + scale_color_viridis_c(option='magma')
  + labs(x=bquote("UMAP"[1]), y=bquote("UMAP"[2]))
  + theme(
    axis.ticks = element_line(color='transparent'),
    panel.background = element_rect(fill=midgut.model.colors.bg[bg_color])
  )
)

plot_midgut_feature_legend <- function(indrop, feature_name) get_legend(
  plot_midgut_feature(indrop, "PCA", "umap", feature_name)
    + geom_point()
    + labs(color="LogNormalize")
    + theme_bw()
    + theme(legend.direction = "horizontal")
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
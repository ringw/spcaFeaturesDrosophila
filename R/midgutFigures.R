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

plot_indrop_pca <- function(indrop) (
  indrop@meta.data %>% cbind(as.data.frame(indrop[['umap']]@cell.embeddings)) %>%
  arrange(pca_classif != 'unknown') %>%
  ggplot(aes(UMAP_1,UMAP_2, color=pca_classif))
  + rasterise(geom_point(shape=20, size=1e-3, show.legend = F), dpi=300)
  + theme_bw()
  + scale_color_manual(values=setNames(midgut.colors,NULL), guide=guide_legend(override.aes=list(size=3)))
  + scale_y_continuous(limits=c(-9.5,NA), expand=rep(0.02,2), breaks=pretty_breaks(3))
  + scale_x_continuous(expand=rep(0.02,2), breaks=pretty_breaks(4))
  + labs(x=bquote("UMAP"[1]), y=bquote("UMAP"[2]))
  + theme(
    axis.ticks = element_line(color='transparent'),
    panel.background = element_rect(fill=midgut.model.colors.bg[1])
  )
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
  + labs(x=bquote("UMAP"[1]), y=bquote("UMAP"[2]))
  + theme(
    axis.ticks = element_line(color='transparent'),
    panel.background = element_rect(fill=midgut.model.colors.bg[2])
  )
)
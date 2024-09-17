save_figure <- function(file, file_plot, width, height, device = cairo_pdf) {
  dir.create(dirname(dirname(file)), showW = FALSE)
  dir.create(dirname(file), showW = FALSE)
  ggplot2::ggsave(file, file_plot, width=width, height=height, dpi=120, device=device)
  return(file)
}

# SPCA tile plot of feature loadings. Custom palette for negative feature
# loadings, which blends into the viridis color palette for nonnegative feature
# loadings.
spc_tile_plot <- function(spca, column, fontface="bold", option="viridis", begin=0, end=1, fontsize=4) {
  if (begin < 0) {
    y <- c(
      seq_gradient_pal("#e0524d", viridis(10)[1])(
        seq(0, 1, length.out=abs(round(begin*100)))[-abs(round(begin*100))]
      ),
      viridis(101)[1:(1 + round(end * 100))]
    )
    scl <- scale_fill_gradientn(
      limits=c(begin, end+1e-10),
      colors=y,
      na.value='#000000'
    )
  } else {
    scl <- scale_fill_viridis_c(
      limits=c(begin,end), breaks=seq(begin,end,length.out=3), option=option, begin=begin, end=end,
      na.value='#000000'
    )
  }
  vals <- spca[, column, drop=T] %>% sort(dec=T) %>% subset(. != 0)
  if ('...' %in% rownames(spca@feature.loadings))
    vals <- c(
      vals[vals > 0],
      `...`=NA,
      vals[vals < 0]
    )
  list(spc.genes = names(vals)) %>% with(
    data.frame(
      gene = factor(spc.genes, levels=spc.genes) %>% display_gene_names,
      loading = as.numeric(vals)
    )
  ) %>%
    ggplot() +
    geom_tile(aes(x='', y=gene, fill=loading), color='black') +
    geom_text(aes(x='', y=gene, label=gene), size=fontsize, fontface=fontface, color='white') +
    scl +
    scale_x_discrete(breaks=NULL, expand=c(0,0)) +
    scale_y_discrete(breaks=NULL, limits=rev, expand=c(0,0)) +
    labs(x=NULL, y=NULL)
}

# Largest PCA (dense, thousands of values) feature loadings. The purpose is to
# plot the 4 most positive and 4 most negative genes, while sharing code with
# the SPCA feature loadings plot.
pc_tile_plot <- function(pca, column, ...) {
  loadings <- pca[, column, drop=T]
  loadings.keep <- union(
    head(order(loadings), 4),
    tail(order(loadings), 4)
  )
  loadings[setdiff(seq_along(loadings), loadings.keep)] <- 0
  pca@feature.loadings[, column] <- loadings
  pca@feature.loadings <- rbind(
    pca@feature.loadings,
    matrix(
      NA,
      ncol = ncol(pca),
      dimnames=list("...", colnames(pca))
    )
  )
  spc_tile_plot(pca, column, ...)
}

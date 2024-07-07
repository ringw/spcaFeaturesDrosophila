save_figure <- function(file, file_plot, width, height, device = cairo_pdf) {
  dir.create(dirname(dirname(file)), showW = FALSE)
  dir.create(dirname(file), showW = FALSE)
  ggplot2::ggsave(file, file_plot, width=width, height=height, dpi=120, device=device)
  return(file)
}

spc_tile_plot <- function(spca, column, fontface="bold", option="viridis", begin=0, end=1, fontsize=4) (list(spc.genes = spca[, column, drop=T] %>% sort(dec=T) %>% subset(. != 0) %>% names) %>% with(
  data.frame(
    gene = factor(spc.genes, levels=spc.genes) %>% display_gene_names,
    loading = spca[spc.genes, column, drop=T]
  )
) %>%
  ggplot
  + geom_tile(aes(x='', y=gene, fill=loading), color='black')
  + geom_text(aes(x='', y=gene, label=gene), size=fontsize, fontface=fontface, color='white')
  + scale_fill_viridis_c(limits=c(begin,end), breaks=seq(begin,end,length.out=3), option=option, begin=begin, end=end)
  + scale_x_discrete(breaks=NULL, expand=c(0,0))
  + scale_y_discrete(breaks=NULL, limits=rev, expand=c(0,0))
  + labs(x=NULL, y=NULL)
  # + spc.theme
  # + theme(legend.position = 'bottom')
)

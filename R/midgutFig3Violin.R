# Violin plot ----

# Stroke color (reduce lightness of EB color so that the line shows up).
midgut.colors.stroke = c(
  'ISC'=hsv(0.88, 0.75, 0.97),
  'EB'=hsv(0.166, 1, 0.85),
  'dEC'=hcl(137, 57, 71),
  'EC'=hcl(253, 74, 48),
  'EC-like'=hcl(230, 34, 72),
  'EE'=hcl(300, 72, 46),
  'copper/iron'=hcl(48, 63, 75),
  'LFC'=hcl(180, 38, 88),
  'cardia'=hcl(9, 66, 32),
  'bg'=hcl(c = 0, l = 87)
)

tiny_violin_plot <- function(seurat, column_to_use, levels_to_use, feature_to_use) {
  data <- seurat %>%
    FetchData(c(column_to_use, feature_to_use)) %>%
    subset(as.character(.[, 1]) %in% levels_to_use)
  colnames(data) <- c("cluster", "embedding")

  # Display 500 of the embedding values.
  jitter_data <- data %>%
    group_by(cluster) %>%
    summarise(embedding = embedding %>% sample(pmin(length(.), 500)))

  rasterise(
    ggplot(data, aes(cluster, embedding, fill=cluster))
    + geom_violin(linewidth=0.5, scale="width")
    + geom_jitter(data=jitter_data, shape=16, stroke=NA, size=0.6),
    dpi = 300
  ) + theme_cowplot() + scale_fill_manual(
    values = c(midgut.colors, setNames(midgut.colors["EC-like"], "")), guide=guide_none()
  ) + scale_y_continuous(
    # Expand the y-axis less since it is a tiny plot
    expand = rep(0.02, 2)
  ) + labs(
    x = NULL, y = NULL
  ) + theme(
    plot.margin = margin(1, 1, 1, 1), aspect.ratio = 0.618
  )
}

# ECDF plot ----

plot_spca_cdf <- function(seurat, column_to_use, levels_to_use, feature_to_use) {
  data <- seurat %>% FetchData(c(column_to_use, feature_to_use)) %>%
    subset(as.character(.[, 1]) %in% levels_to_use)
  spca_samples <- (
    data[, 2] %>%
      split(data[, 1])
  )[levels_to_use] %>%
    sapply(ecdf, simplify=FALSE) %>%
    enframe(name="cluster")
  plot_data <- spca_samples %>%
    full_join(
      data.frame(embedding = seq(min(data[, 2]) + 1e-6, max(data[, 2]), length.out=100)),
      by=character()
    ) %>%
    rowwise %>%
    mutate(
      value = value(embedding)
    )
  # Insert zeroes into the plot data
  plot_data <- rbind(
    data.frame(cluster = levels_to_use, value = 0, embedding = -Inf),
    plot_data
  )

  ggplot(
    plot_data, aes(embedding, value, color=cluster, group=cluster)
  ) + geom_line() + scale_color_manual(
    values=midgut.colors.stroke, guide=guide_none()
  ) + coord_cartesian(
    NULL, c(0,1),
    expand = FALSE
  ) + scale_y_continuous(
    labels=percent
  ) + labs(
    x=NULL, y=NULL
  ) + theme_cowplot() + theme(
    plot.margin = margin(1, 1, 1, 1), aspect.ratio = 0.618,
    # y-labels still didn't have enough space at the top with 1pt margin.
    axis.text.y = element_text(vjust = 0.8)
  )
}
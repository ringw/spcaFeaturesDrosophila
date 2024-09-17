elbow_plot <- function(values) {
  enframe(
    values, "Component", "Explained Variance"
  ) %>%
    ggplot(aes(Component, `Explained Variance`)) +
    geom_point(size = 0.75, color = "#141c68") +
    scale_x_continuous(expand = expansion(add = 0.9)) +
    coord_cartesian(NULL, c(1, NA)) +
    theme_cowplot() +
    theme(
      aspect.ratio = 1 / 3,
      plot.margin = margin(12, 4, -2, 4),
      axis.text = element_text(size = 14)
    )
}

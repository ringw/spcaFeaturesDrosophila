make_dimreduc_img <- function(dimreduc, d=48, k=8) {
  features <- apply(
    dimreduc@feature.loadings[, seq(d)],
    2,
    \(v) v %>%
      sort(decreasing=T) %>%
      head(k) %>%
      c(rep(0, k - length(.))) %>%
      setNames(NULL)
  ) %>%
    t %>%
    pmax(0)
  features <- features %>% replace(. == 0, NA)
  dimnames(features) <- list(component=NULL, feature=NULL)
  melt(features, value.name = "loading") %>%
    ggplot(
      aes(feature, component, fill=loading)
    ) + geom_tile() + scale_y_reverse() + coord_cartesian(
    expand=F
  ) + scale_fill_viridis_c(
    limits=c(0, 1), na.value=hsv(0, 0.05, 0.85), guide=guide_none()
  ) + labs(
    x = 'Support of PC (sorted genes)', y = 'PCs (descending variance)'
  )
}
clamp <- function(numeric, lower, upper) pmin(pmax(numeric, lower), upper)

covarianceToGraph <- function(mat, targetSparsity=1e-3) {
  off_diagonals <- mat[-seq(1, nrow(mat)^2, by=1+nrow(mat))]
  cov_scale_factor <- 1 / quantile(abs(off_diagonals), 1-targetSparsity)
  clamp(mat * cov_scale_factor, 0, 1)
}

graphLaplacian <- function(mat) {
  diagAdjust <- rowSums(mat)
  -mat + diag(x = diagAdjust)
}

plotGraphLaplacianMatrices <- function(matlist, ...) {
  plot_x <- seq(0, 50, by=0.1)
  ggarrange(
    plotlist = mapply(
      \(m, tag) m %>%
        covarianceToGraph %>%
        graphLaplacian %>%
        eigen(sym=T) %>%
        with(values) %>%
        density(bw=1) %>%
        with(approx(x, y, plot_x)) %>%
        as.data.frame %>%
        ggplot(aes(x, y))
        + geom_line()
        + coord_cartesian(expand=F)
        + labs(x = "N = 1000 Bandwidth = 1", y = "Density", tag = tag)
      ,
      matlist,
      LETTERS %>% head(length(matlist)),
      SIMPLIFY=FALSE
    ) # ,
    # ...
  )
}
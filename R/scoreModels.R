score_univariate_mixture_model <- function(assay, features_to_use, feature_loadings, categorical) {
  scale.data <- assay %>%
    ScaleData(features = features_to_use, verb = F) %>%
    GetAssayData(slot = "scale.data") %>%
    t %>%
    subset(select = features_to_use)
  univariate_mixture <- scale.data %>%
    as.data.frame %>%
    split(categorical) %>%
    mapply(
      \(name, data) list(
        ident = name,
        mean = colMeans(data),
        s2 = eigen(var(data)) %>% with(
            vectors %*% Diagonal(x = pmax(values, 1e-10)) %*% t(vectors)
          ) %>%
          as.matrix %>%
          matrix(
            nrow = nrow(.),
            dimnames = rep(list(features_to_use), 2)
          )
      ),
      levels(categorical),
      .,
      SIMPLIFY=FALSE
    )

  pc1 <- enframe(feature_loadings) %>%
    inner_join(data.frame(name = features_to_use, i = seq_along(features_to_use)), "name") %>%
    with(sparseVector(x = value, i = i, length = length(features_to_use)))
  pc1_residuals <- scale.data - ((scale.data %*% pc1) %*% t(pc1)) %>% as.matrix
  univariate_mixture %>%
    sapply(
      \(cluster_stats) cluster_stats %>%
        with(
          log(mean(as.character(categorical) == ident))
          + (
            pc1_residuals %>%
              # subset(as.character(categorical) == ident) %>%
              dmvnorm(
                mean,
                s2,
                log = TRUE
              )
          )
        )
    ) %>%
    # sapply places each result (mixture component) into a column: For each
    # cell, take the mixture likelihoods (pre-multiplied by component
    # probability); this is not performed in log space.
    rowLogSumExps %>%
    # We now have log-likelihood for each cell in the residuals mixture model.
    # Add log-likelihoods in log space.
    sum
}
# with(as.data.frame(acc[['spca']]@cell.embeddings), cor(SPARSE_1, SPARSE_14))
# with(as.data.frame(acc[['spca']]@cell.embeddings), cor(SPARSE_3, SPARSE_26))
# with(as.data.frame(acc[['spca']]@cell.embeddings), cor(SPARSE_12, SPARSE_11))

acc_blend_features <- function(acc, feature1, feature2) {
  n_colors <- 256
  scale1 = viridis(n_colors, end=0.9, option="viridis")
  scale2 = viridis(n_colors, begin=0.1, option="magma")
  blend <- mixcolor(
    0.6,
    hex2RGB(scale1) %>% as("LAB"), hex2RGB(scale2) %>% as("LAB")
  ) %>% hex

  features <- as.data.frame(acc[['umap.spca']]@cell.embeddings) %>%
    cbind(
      acc %>% FetchData(c(feature1, feature2))
    )

  feature_maxs <- colQuantiles(features %>% as.matrix, probs=0.95)
  feature_mapping <- as.matrix(features) %>%
    `%*%`(matrix(diag(1 / feature_maxs), ncol=ncol(.), dimnames=rep(list(colnames(.)), 2)))
  feature_mapping <- round(
    feature_mapping * (n_colors - 1) + 1
  ) %>% pmin(n_colors)
  colors <- mixcolor(
    0.6,
    hex2RGB(scale1[feature_mapping[, feature1]]) %>% as("LAB"),
    hex2RGB(scale2[feature_mapping[, feature2]]) %>% as("LAB")
  ) %>% hex
  features <- features %>% cbind(color = colors)

  ggplot(features, aes(UMAP_1, UMAP_2)) + geom_point(color=features$color, size=2) + geom_tile(aes(color=C1), data.frame(UMAP_1=0, UMAP_2=0, C1=0), alpha=0, width=0, height=0) + scale_color_viridis_c(limits = c(0, feature_maxs[feature1]), end=0.9) + new_scale_color() + geom_tile(aes(color=C1), data.frame(UMAP_1=0, UMAP_2=0, C1=0), alpha=0, width=0, height=0) + scale_color_gradientn(limits = c(0, min(feature_maxs[c(feature1,feature2)])), colors = blend) + new_scale_color() + geom_tile(aes(color=C1), data.frame(UMAP_1=0, UMAP_2=0, C1=0), alpha=0, width=0, height=0) + scale_color_viridis_c(limits = c(0, feature_maxs[feature2]), option = "magma", begin=0.1)
}
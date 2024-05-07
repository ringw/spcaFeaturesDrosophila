seurat_elasticnet_spca <- function(seurat, varnum, n_components, assay = "RNA") {
  features <- matrix(
    nrow = nrow(GetAssayData(seurat, "scale.data", assay=assay)),
    ncol = 0,
    dimnames = list(
      rownames(GetAssayData(seurat, "scale.data", assay=assay)),
      NULL
    )
  )
  scale.data <- t(GetAssayData(seurat, assay=assay, slot="scale.data"))
  for (i in seq(n_components)) {
    elasticnet_obj <- elasticnet::spca(
      scale.data, 1, varnum, sparse="varnum"
    )
    v <- elasticnet_obj$loadings[, 1] %>%
      `*`(elasticnet_obj$loadings[, 1] %>% subset(. != 0) %>% sign %>% median)
    features <- features %>% cbind(v)
    scale.data <- scale.data - (scale.data %*% v) %*% t(v)
    colnames(features)[i] <- str_glue("ELASTICNET_{i}")
  }
  features
}
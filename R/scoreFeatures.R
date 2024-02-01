# Select SPCA features given a continuous or discrete cell feature of interest.
querySpcaFeatures <- function(spca, query, cells.select = NULL) {
  # Expect seurat with DimReduc named "spca", feature as the arithmetic mean of query features.
  if (type_sum(query) == "chr") {
    scaleNonneg <- function(mat) {
      mat <- as.matrix(mat)
      (
        DelayedArray(mat) * RleArray(
          Rle(
            1 / colSds(mat),
            rep(nrow(mat), ncol(mat))
          ),
          dim(mat)
        )
      ) %>% as.matrix
    }
    query <- FetchData(spca, query) %>% scaleNonneg %>% rowMeans
  }
  if (attr(spca, "class") == "Seurat")
    spca <- spca[["spca"]]
  if (!is.null(cells.select)) {
    if (min(query) < 0) {
      warning(
        "Centered feature passed to querySpcaFeatures. Making feature nonnegative."
      )
      query <- (query - min(query))
    }
    query <- query * (rownames(spca) %in% cells.select)
  }
  scores <- data.frame(
    cor = cor(spca@cell.embeddings, query)
  )
  scores$score <- scores$cor * spca@stdev
  scores %>% arrange(desc(abs(score)))
}
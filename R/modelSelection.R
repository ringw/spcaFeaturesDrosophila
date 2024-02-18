spcaStats <- function(seurat, reduction.name = "spca") {
  spca <- seurat[[reduction.name]]
  apply(
    spca@feature.loadings,
    2,
    \(fe) {
      fe <- fe %>% subset(. != 0)
      # TODO: Re-run spca, stop using tar_cue never, and ensure that it has the
      # correct feature names!
      fe <- fe[names(fe) %in% rownames(seurat[['RNA']])]
      abundance <- colMeans(
        seurat[['RNA']]@counts[names(fe), ] %>% as.matrix %>% t
        /
        seurat$nCount_RNA
      )
      sim <- rpois(
        length(fe) * ncol(seurat),
        abundance %*% t(seurat$nCount_RNA)
      ) %>%
        matrix(
          length(fe), ncol(seurat),
          dimnames = list(names(fe), Cells(seurat))
        ) %>%
        rbind(
          others = seurat$nCount_RNA - colSums(.)
        ) %>%
        NormalizeData(verbose = F)
      vars <- rowVars(sim[-nrow(sim), ]) %>%
        `/`(seurat[['RNA']]@meta.features[names(fe), "vst.variance"])
      sqrt(sum(vars * fe^2))
    }
  ) %>% zoo::rollapply(10, median, fill=c("extend","extend","extend"))
}
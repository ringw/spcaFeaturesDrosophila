run_fgsea_multiple_times <- function(
  indrop.deg,
  indrop.gsea.genes,
  gene_association,
  go_basic_biological_process
) {
  pathways <- with(gene_association, split(V3, V5)) %>%
    sapply(unique, simplify=F)
  pathways <- pathways[sort(unique(go_basic_biological_process$id))]
  sapply(
    indrop.deg,
    \(model) fgsea(
      pathways,
      setNames(
        model$map[indrop.gsea.genes, 2],
        names(indrop.gsea.genes)
      ),
      10
    ),
    simplify=FALSE
  )
}

plot_fgsea_multiple_times <- function(
  indrop.deg,
  indrop.gsea.genes,
  gene_association
) {
  pathways <- with(gene_association, split(V3, V5)) %>%
    sapply(unique, simplify=F)
  pathways <- pathways[order(names(pathways))]
  plot_data <- sapply(
    pathways,
    \(pathway) sapply(
      indrop.deg,
      \(model) {
        df <- plotEnrichment(
          pathway,
          setNames(model$map[indrop.gsea.genes, 2], names(indrop.gsea.genes))
        )$data
        # df is a series of enrichment score before-after points on the
        # empirical cdf. We don't want to draw a stroke on the odd-located
        # df rows (which are x = 0, or rank - 1).
        df$tick <- c(FALSE, TRUE)
        df
      },
      simplify=FALSE
    ),
    simplify=FALSE
  )
}

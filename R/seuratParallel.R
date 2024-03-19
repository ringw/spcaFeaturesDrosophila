replicate_spca_model_parallel <- function(
  indrop, memory_cgroups, seed,
  varnum, npcs, eigen_gap, search_cap, do.correct.elbow
) {
  plan(multisession, workers=4)

  if (is.null(memory_cgroups)) memory_cgroups <- rep(list(NULL), 4)
  result = list()
  covar <- indrop@misc$covar
  for (rep in 1:4) {
    result[[rep]] <- future(
      with_seed(
        seed[rep],
        seurat_spca_compute_feature_loadings(
          covar,
          varnum=varnum,
          npcs=npcs,
          search_cap=search_cap,
          eigen_gap=eigen_gap,
          cgroup=memory_cgroups[rep]
        )
      )
    )
  }
  for (rep in 1:4) {
    result[[rep]] <- value(result[[rep]]) %>%
      seurat_spca_from_feature_loadings(indrop, "RNA", do.correct.elbow = TRUE)
  }
  result
}
seurat_joint_spca <- function(
  seurat, alpha, ridge_alpha, n_components, assay = "RNA"
) {
  obs_data <- t(chol(seurat@misc$covar + diag(x=rep(1e-8, 1000))))
  basiliskProc <- basiliskStart(scikit_learn_env)
  calculate_spca_sklearn <- function() {
    np <- reticulate::import("numpy")
    np$random$seed(4L)
    skdecomp <- reticulate::import("sklearn.decomposition")
    fit_obj <- skdecomp$SparsePCA(
      alpha=alpha, ridge_alpha=ridge_alpha, n_components=as.integer(n_components)
    )
    fit_obj$fit(obs_data)
    fit_obj$components_
  }
  componentsMat <- tryCatch(
    basiliskRun(basiliskProc, calculate_spca_sklearn),
    finally = basiliskStop(basiliskProc)
  )
  componentsMat %>%
    matrix(
      nrow = nrow(.),
      ncol = ncol(.),
      dimnames = list(
        str_glue("SPARSE_{1:n_components}"),
        colnames(seurat@misc$covar)
      )
    )
}
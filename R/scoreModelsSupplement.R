# Supplemental figures for scoring model selection.
score_models_supplement <- list(
  tar_target(
    fig.laplacian,
    save_figure(
      "figure/Model-Selection/Graph-Laplacian-Spectrum.pdf",
      plotGraphLaplacianMatrices(
        list(A=indrop@misc$covar, B=tenx@misc$covar, C=acc@misc$covar)
      ),
      width = 6, height = 4.5
    )
  )
)
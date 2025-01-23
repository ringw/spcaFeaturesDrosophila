human_cell_landscape_pca <- list(
  tar_target(
    human.cell.landscape.pattern,
    tribble(
      ~label, ~filename, ~param.k,
      "AdultAscendingColon", "new/AdultAscendingColon1.rmbatchdge.txt.gz", 3,
      "AdultCerebellum", "new/AdultCerebellum1.rmbatchdge.txt.gz", 7,
      "AdultMuscle", "new/AdultEsophagus1.rmbatchdge.txt.gz", 6,
      "AdultPancreas", "new/AdultPancreas1.rmbatchdge.txt.gz", 10,
      "AdultPleura", "new/AdultPleura1.rmbatchdge.txt.gz", 3,
      "AdultTemporalLobe", "new/AdultTemporalLobe1.rmbatchdge.txt.gz", 3,
      "AdultUterus", "new/AdultUterus1.rmbatchdge.txt.gz", 3,
      "FetalMuscle", "new/FetalMuscle1.rmbatchdge.txt.gz", 3,
      "Placenta", "new/Placenta1.rmbatchdge.txt.gz", 3,
    )
  ),
  tar_file(
    human.cell.landscape.matrix,
    human.cell.landscape.pattern$filename,
    pattern = map(human.cell.landscape.pattern)
  ),
  tar_target(
    human.cell.landscape.pca,
    {
      seurat <- human.cell.landscape.matrix %>%
        read.table() %>%
        CreateSeuratObject() %>%
        `Idents<-`(value = human.cell.landscape.pattern$label) %>%
        NormalizeData(verb = F) %>%
        FindVariableFeatures(nfeatures = 1500) %>%
        ScaleData(verb = F) %>%
        RunPCA(verb = F) %>%
        RunUMAP(dims = 1:10)
      seurat@misc$covar <- tcrossprod(
        seurat[["RNA"]]@scale.data[VariableFeatures(seurat), ]
      ) /
        ncol(seurat)
      seurat
    },
    pattern = map(human.cell.landscape.pattern, human.cell.landscape.matrix)
  ),
  tar_target(
    human.cell.landscape.param.k,
    tibble(k = 1:20)
  ),
  tar_target(
    human.cell.landscape.select.k,
    tibble(
      human.cell.landscape.pattern,
      human.cell.landscape.param.k,
      model = human.cell.landscape.pca %>%
        seurat_spca(
          "covar",
          varnum=k,
          npcs=3,
          eigen_gap=0.01, search_cap=500000
        ) %>%
        list()
    ),
    pattern = cross(map(human.cell.landscape.pattern, human.cell.landscape.pca), human.cell.landscape.param.k)
  ),
  tar_target(
    human.cell.landscape.spca,
    tibble(
      human.cell.landscape.pattern,
      model = human.cell.landscape.pca %>%
        seurat_spca(
          "covar",
          varnum=human.cell.landscape.pattern$param.k,
          npcs=60,
          eigen_gap=0.01, search_cap=150000
        ) %>%
        list()
    ),
    pattern = map(human.cell.landscape.pattern, human.cell.landscape.pca)
  ),
  tar_target(
    human.cell.landscape.k,
    human.cell.landscape.select.k %>%
      arrange(label, k) %>%
      rowwise() %>%
      summarise(label, k, explained_variance = model@stdev[1]^2)
  ),
  tar_target(
    human.cell.landscape,
    tibble(
      human.cell.landscape.pattern["label"],
      seurat = {
        seurat <- human.cell.landscape.pca
        seurat[["spca"]] <- human.cell.landscape.spca$model[[1]]
        seurat %>%
          RunUMAP(reduction = "spca", dims = 1:50, reduction.name = "umap.spca") %>%
          FindNeighbors(reduction = "spca", dims = 1:50) %>%
          FindClusters()
      } %>%
        list(),
    ),
    pattern = map(human.cell.landscape.pattern, human.cell.landscape.pca, human.cell.landscape.spca)
  )
)

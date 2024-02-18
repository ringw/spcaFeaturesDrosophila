readMouseExpressionMatrix <- function(hdf5_path) {
  cellNames <- DelayedArray(HDF5ArraySeed(hdf5_path, "/data/samples")) %>% as.character
  geneNames <- DelayedArray(HDF5ArraySeed(hdf5_path, "/data/gene")) %>% as.character
  counts <- DelayedArray(HDF5ArraySeed(hdf5_path, "/data/counts"))
  dimnames(counts) <- list(cell = cellNames, gene = geneNames)
  counts
}

filterMouseExpressionMatrix <- function(counts) {
  ncells <- colSums(counts != 0)
  numi <- rowSums(counts)
  ncells_min <- 1000
  numi_min <- 5000
  counts <- counts[numi >= numi_min, ncells >= ncells_min]
}

computeMouseHVGs <- function(counts) {
  # cells <- sample(rownames(counts)) %>%
  #   split(cut(seq_along(.), seq(0, length(.) + shard_size - 1, by=10000))) %>%
  #   sapply(\(v) rownames(counts) %>% subset(. %in% v))
  # names(cells) <- paste0("shard", seq_along(cells))

  # Switch to Seurat conventions: genes in rows.
  counts <- t(counts)
  cells_to_analyze <- 100000
  if (ncol(counts) > cells_to_analyze)
    counts <- counts[, sample(ncol(counts), cells_to_analyze)]
  counts <- counts %>% as("sparseMatrix")
  # change this to 10K
  counts %>%
    NormalizeData(scale.factor = 1000 * 1000) %>%
    FindVariableFeatures
}

scale_genes_from_list <- function(counts, numi, cell_list, gene_list) {
  # write(setdiff(gene_list, colnames(counts)), stdout())
  # stopifnot(0 == length(setdiff(gene_list, colnames(counts))))
  data <- counts[cell_list, gene_list] %>% as.matrix
  numi <- numi[cell_list]
  # Cells in rows; broadcast the numi along each row.
  data <- log1p(data * 1000 * 1000 / numi)
  # Scale by sd; do not center as var() will apply centering.
  data <- as.matrix(
    as(data, "RleArray") * RleArray(
      Rle(1 / colSds(data), rep(nrow(data), ncol(data))),
      dim(data),
      dimnames(data)
    )
  )
  data %>% as("sparseMatrix")
}

assign_matrix_blocks <- function(namevec, blocks) {
  m <- matrix(NA, nrow = length(namevec), ncol = length(namevec), dim = rep(list(namevec), 2))
  for (b in blocks)
    m[rownames(b), colnames(b)] <- b
  m
}

seurat_spca_only <- function(feature.loadings, cell.embeddings, rna.counts = NULL) {
  assay.features <- paste0("SPC", seq(ncol(feature.loadings)))
  colnames(feature.loadings) <- assay.features
  colnames(cell.embeddings) <- assay.features

  # TODO: Fix the sign of the features
  seurat <- CreateAssayObject(data = t(cell.embeddings)) %>%
    CreateSeuratObject(assay = "SPCA")

  dimreduc.features <- paste0("SPARSE_", seq(ncol(feature.loadings)))
  colnames(feature.loadings) <- dimreduc.features
  colnames(cell.embeddings) <- dimreduc.features
  seurat[['spca']] <- CreateDimReducObject(
    loadings = feature.loadings, embeddings = cell.embeddings,
    assay = 'SPCA'
  )

  if (!is.null(rna.counts)) {
    seurat[['RNA']] <- CreateAssayObject(counts = rna.counts)
  }

  seurat
}

seurat_normalize_by_numi <- function(seurat, numi) {
  seurat[['RNA']]@data <- as.matrix(
    as(as.matrix(seurat[['RNA']]@counts), "RleArray")
    * RleArray(
      Rle(
        10000 / numi,
        rep(nrow(seurat[['RNA']]@counts), ncol(seurat[['RNA']]@counts))
      ),
      dim = dim(seurat[['RNA']]@counts)
    )
  ) %>% log1p
  seurat
}
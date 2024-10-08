make_rownames.acc <- function(counts_file) {
  acc <- (
    read.table(
      counts_file, row.names = 1, header = T
    ) %>%
      subset(select = -c(Chr, Start, End, Strand, Length))
  )
  mapIds(
    org.Hs.eg.db, rownames(acc) %>% str_replace("\\.[0-9]+", ""), "SYMBOL", "ENSEMBL"
  ) %>%
    replace(
      is.na(.), names(.)[is.na(.)]
    ) %>%
    replace(
      duplicated(.), names(.)[duplicated(.)]
    ) %>%
    setNames(NULL)
}

# Write chr, start, end for our genes.
write_gene_order.acc <- function(counts_file, rownames.acc, output_file) {
  gene_order <- (
    read.table(
      counts_file, row.names = 1, header = T
    ) %>%
      subset(select = c(Chr, Start, End))
  )
  dimnames(gene_order) = list(
    rownames.acc,
    c('chr', 'start', 'stop')
  )
  gene_order$chr <- gene_order$chr %>% str_extract("[^;]+")
  gene_order[, 2:3] <- gene_order[, 2:3] %>%
    apply(
      1,
      \(chars) strsplit(chars, ";") %>%
        unlist %>%
        as.numeric %>%
        range
    ) %>% t
  write.table(gene_order, output_file, col.names=F, quote = FALSE, sep = "\t")
  output_file
}

preprocess_acc <- function(counts_file, rownames.acc) {
  acc <- (
    read.table(
      counts_file, row.names = 1, header = T
    ) %>%
      subset(select = -c(Chr, Start, End, Strand, Length))
  )
  colnames(acc) <- colnames(acc) %>% str_replace(
    str_escape(".sort.bam.gene_counts.tsv"),
    ""
  )
  mt.genes.logical <- mapIds(
    org.Hs.eg.db, rownames(acc) %>% str_replace("\\.[0-9]+", ""), "CHR", "ENSEMBL"
  ) %>%
    `==`("MT") %>%
    replace(is.na(.), F)
  rownames(acc) <- rownames.acc
  acc <- acc %>% CreateSeuratObject()
  acc$pct.mito <- acc %>%
    PercentageFeatureSet(features = rownames(acc)[mt.genes.logical])
  acc <- acc %>% subset(
    cells = Cells(acc) %>% subset(acc$nCount_RNA >= 2000 & acc$pct.mito <= 30)
  )

  acc$individual <- str_extract(Cells(acc), "ACC[0-9]+") %>%
    factor(paste0("ACC", c(2, 5, 7, 15, 19, 21, 22)))
  acc$recurrence <- str_extract(Cells(acc), "(?<=ACC5.)P[0-9]+") %>%
    recode(P1="primary", P2="primary", P11="recurrence")

  # log2(TPM+1) / log2(E)
  # -> ln(TPM+1)
  acc = acc %>% NormalizeData(scale.factor = 1000 * 1000)
  acc = acc %>% SetAssayData('data', GetAssayData(acc[['RNA']], 'data') / log(2))
  interesting.features <- 'TP63'
  acc = acc %>% FindVariableFeatures %>%
    ScaleData(features = union(VariableFeatures(.), interesting.features))
  spca.features <- VariableFeatures(acc) %>% head(1000)
  req.features <- setdiff(interesting.features, spca.features)
  spca.features <- spca.features %>% head(1000 - length(req.features)) %>%
    union(req.features)
  acc@misc$covar = acc[['RNA']]@scale.data[spca.features, ] %>%
    t %>%
    as.matrix %>%
    var

  acc
}

process_acc_spca <- function(acc, spca) {
  acc <- acc %>%
    spca_with_centered_umap(spca, seed.use=12, umap_transform=diag(c(-1, -1)))
  # Major cell types of interest on the left-hand side instead:
  # acc[['umap.spca']]@cell.embeddings[, 'UMAP_1'] <- (
  #   -acc[['umap.spca']]@cell.embeddings[, 'UMAP_1']
  # )
  # The paper refers to NOTCH1 in luminal cells.
  # We will use feature loadings with a sum of squares of 1, rather than a sum
  # of 1 (which could produce a weighted mean).
  acc[['N1N2']] <- colMeans(acc[['RNA']][c('NOTCH1','NOTCH2'), ]) * sqrt(2)
  # Parikh preferred DLL1, JAG2, by correlation with the overall myoepithelial
  # transcriptional program (with ACTA2, TP63 having large feature loadings) and
  # rejecting JAG1. We will start here from a place of less knowledge and apply
  # all canonical Notch ligands.
  acc[['ligands']] <- colMeans(acc[['RNA']][c('DLL1','DLL3','DLL4','JAG1','JAG2'), ]) * sqrt(5)
  # NOTCH3 response, and related genes.
  acc@meta.data$`N3+` <- colMeans(
    acc[['RNA']][c('NOTCH3','HES4','HEY1','HEY2','NRARP'),]
  ) * sqrt(5)
  acc[['targets']] <- colMeans(
    acc[['RNA']][c('NOTCH3','HES4','HEY1','HEY2','NRARP'),]
  ) * sqrt(5)
  acc$spca_clusters = (
    acc %>% FindNeighbors(red = "spca", dims = 1:33) %>% FindClusters(res = 0.5)
  )$seurat_clusters
  acc$spca_coarse = (
    acc %>% FindNeighbors(red = "spca", dims = 1:33) %>% FindClusters(res = 0.1)
  )$seurat_clusters

  RNA.scaled <- acc[["RNA"]]
  Key(RNA.scaled) <- "rnascaled_"
  acc[["RNA.scaled"]] <- RNA.scaled
  DefaultAssay(acc) <- "RNA.scaled"
  acc <- acc %>% ScaleData(vars.to.regress = c("nFeature_RNA", "pct.mito"))
  acc <- acc %>% RunPCA(verb = F)
  acc <- acc %>% RunUMAP(dims = 1:15)
  acc$pca_clusters = (
    acc %>% FindNeighbors(dims = 1:10) %>% FindClusters(res = 0.5)
  )$seurat_clusters
  acc$pca_coarse = (
    acc %>% FindNeighbors(dims = 1:10) %>% FindClusters(res = 0.1)
  )$seurat_clusters

  acc
}

write_seurat_column <- function(seurat, column_name, output_path) {
  df <- data.frame(
    cell = Cells(seurat),
    ident = seurat@meta.data[, column_name]
  )
  write.table(
    df,
    output_path,
    row.names = F,
    col.names = F,
    quote = FALSE,
    sep = "\t"
  )
  output_path
}

spcaFeatureTable <- function(acc, group.by="individual") {
  features <- acc[['spca']]@cell.embeddings >= 5
  features <- features[, 1:30]
  apply(
    features %>% cbind(total_cells = TRUE),
    2,
    \(v) table(
      is_active = v,
      feature = FetchData(acc, group.by) %>% as.data.frame %>% pull(group.by)
    )["TRUE", ]
  )
}

acc_to_caf_spca <- function(acc) {
  # Rectangular area identified from accFigures.
  caf <- acc %>% subset(
    cells = Cells(acc) %>% subset(
      FetchData(acc, 'nCount_RNA') %>%
        cbind(acc[['umap.spca']]@cell.embeddings) %>%
        with(
          between(UMAP_1, -2.75, 6.5) & between(UMAP_2, -9.75, -4)
          & nCount_RNA >= 500000
        )
    )
  )
  caf <- caf %>%
    RunUMAP(reduction='spca', dims=1:50, reduction.name='umap.spca') %>%
    FindNeighbors(reduction='spca', dims=1:30) %>%
    FindClusters(graph.name = 'RNA_snn', res = 0.5)
  caf[['spca']]@misc$pct.expl <- manova(
    spca ~ ident,
    list(spca = caf[['spca']]@cell.embeddings, ident = Idents(caf))
  ) %>% with(
    colSds(fitted.values, useNames = TRUE) / (colSds(fitted.values) + colSds(residuals))
  )
  subset_dims <- (caf[['spca']]@misc$pct.expl >= 0.55) %>%
    which %>% names
  spca.subset <- CreateDimReducObject(
    caf[['spca']]@cell.embeddings[, subset_dims],
    caf[['spca']]@feature.loadings[, subset_dims]
  )
  colnames(spca.subset@cell.embeddings) <- colnames(
    spca.subset@cell.embeddings
  ) %>%
    str_replace('SPARSE_', 'SPC')
  colnames(spca.subset@feature.loadings) <- colnames(
    spca.subset@feature.loadings
  ) %>%
    str_replace('SPARSE_', 'SPC')
  Key(spca.subset) <- 'spcasubset_'
  caf[['spca.subset']] <- spca.subset
  caf <- caf %>% AddMetaData(as.data.frame(spca.subset@cell.embeddings))
  Idents(caf) <- Idents(caf) %>%
    fct_recode(pCAF="4", iCAF2="0", iCAF="1", dCAF="2", mCAF="3") %>%
    fct_relevel(c("dCAF", "iCAF", "iCAF2", "mCAF", "pCAF"))
  caf
}

acc_call_idents <- function(acc, caf) {
  acc <- acc %>% acc_add_components_from_umap
  Idents(acc) <- Idents(acc) %>% fct_recode(otherCAF="CAF")
  Idents(acc, cells = Cells(caf)) <- Idents(caf)
  only.acc <- acc %>% subset(idents = "ACC")
  only.acc <- only.acc %>% FindNeighbors(reduction="spca", dims=1:30) %>%
    FindClusters(graph.name = "RNA_snn", res = 0.1)
  Idents(only.acc) <- Idents(only.acc) %>% fct_recode(
    luminal="1", myoepithelial="0", otherACC="2"
  )
  Idents(acc, cells = Cells(only.acc)) <- Idents(only.acc)
  acc
}

acc_glm <- function(counts, colData) {
  g <- glm_gp(
    counts,
    ~ 0 + ident,
    colData,
    on_disk = FALSE,
    verbose = TRUE,
    size_factors = colData$size_factor
  )
  # Dense to sparse for storage
  assay(g$data) <- counts
  g
}

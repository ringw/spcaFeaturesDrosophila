preprocess_acc <- function(counts_file) {
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
  rownames(acc) <- mapIds(
    org.Hs.eg.db, rownames(acc) %>% str_replace("\\.[0-9]+", ""), "SYMBOL", "ENSEMBL"
  ) %>%
    replace(
      is.na(.), names(.)[is.na(.)]
    ) %>%
    replace(
      duplicated(.), names(.)[duplicated(.)]
    ) %>%
    setNames(NULL)
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
  # Load SPCA from file
  acc.spca.feature <- read.csv(spca_file, check.names=F) %>% as.matrix
  spca.emb = matrix(nrow = ncol(acc), ncol = nrow(acc.spca.feature), dimnames = list(Cells(acc), paste0("SPCA_", 1:nrow(acc.spca.feature))))
  spca.data = GetAssayData(acc, slot = 'scale.data')
  spca.data = spca.data[colnames(acc.spca.feature), ]
  for (i in 1:nrow(acc.spca.feature)) {
    # PCA loadings given the subset of k variables only.
    u = t(spca.data) %*% acc.spca.feature[i, ]
    # Deflate only the k variables
    spca.emb[, i] = u
    # (I - x %*% t(x)) %*% spca.data
    spca.data = spca.data - acc.spca.feature[i, ] %*% t(u)
  }
  acc[["spca"]] = CreateDimReducObject(embeddings = spca.emb, key = "SPCA_", assay = "RNA")
  acc[["spca"]]@stdev = apply(spca.emb, 2, \(x) norm(x, type='2')) / sqrt(ncol(acc))
  k = sum(acc.spca.feature[1,] != 0)

  acc = acc %>% RunUMAP(reduction='spca', dims=1:50, reduction.name='umap.spca')
  acc$Notch_Score <- colMeans(
    acc[['RNA']][c('NRARP','NOTCH3','HES4','HEY1','HEY2'),]
  ) %>% scale(scale=F)
  acc$Notch_1_2 <- colMeans(acc[['RNA']][c('NOTCH1','NOTCH2'), ]) %>% scale(scale=F)
  acc$spca_clusters = (
    acc %>% FindNeighbors(red = "spca", dims = 1:33) %>% FindClusters(res = 0.5)
  )$seurat_clusters

  RNA.scaled <- acc[["RNA"]]
  Key(RNA.scaled) <- "rnascaled_"
  acc[["RNA.scaled"]] <- RNA.scaled
  DefaultAssay(acc) <- "RNA.scaled"
  acc <- acc %>% ScaleData(vars.to.regress = c("nFeature_RNA", "pct.mito"))
  acc <- acc %>% RunPCA(verb = F)
  acc <- acc %>% RunUMAP(dims = 1:15)

  # Y = beta M^T
  # beta -> beta * betamat.transform
  # M -> M %*% solve(t(betamat.transform))
  # We want to subtract the average of the other cell types from cluster '2'.
  # A gene is a row of "beta", and we want column operations.
  mm = model.matrix(~ spca_clusters, acc@meta.data)
  betamat.transform = diag(nrow=12)
  dimnames(betamat.transform) = rep(list(colnames(mm)), 2)
  subtract_coefs = c('spca_clusters1','spca_clusters7','spca_clusters8')
  betamat.transform[subtract_coefs, 'spca_clusters2'] = -1 / (length(subtract_coefs)+1)
  mm = mm %*% solve(t(betamat.transform))
  colnames(mm)[3] = 'cls_2_vs_average'
  lm(
    spca ~ 0 + mm,
    list(spca = acc[['spca']]@cell.embeddings, mm = mm)
  )
  lm(
    spca ~ recurrence,
    list(spca = acc[['spca']]@cell.embeddings, recurrence = acc$recurrence),
    subset = !is.na(acc$recurrence)
  )

  acc
}
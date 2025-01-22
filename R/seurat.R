build_meta_features <- function(counts_file, flybase_file, gtf_file) {
  flybase = read.table(
    flybase_file,
    col.names=c('symbol', 'organism', 'flybase', 'secondary_flybases', 'ann', 'secondary_ann'),
    sep='\t',
    quote=''
  )
  counts = read.csv(counts_file, row.names=1, check.names=F)

  # Now guess the transcript lengths for the genes (longest isoform).
  reference = read.table(
    gtf_file,
    sep = '\t',
    col.names = c('chr', 'source', 'type', 'start', 'end', 'sc', 'strand', 'fr', 'annotation'),
    header = F,
    quote = ''
  ) %>% subset(grepl('RNA', type))
  reference$gene_id = str_extract(reference$annotation, 'gene_id "([^"]+)"', group=1)
  reference$gene_symbol = str_extract(reference$annotation, 'gene_symbol "([^"]+)"', group=1)

  meta.features = data.frame(gene_id = rownames(counts)) %>% left_join(
    data.frame(
      gene_symbol = flybase$symbol,
      gene_id = flybase$flybase
    ),
    'gene_id'
  )
  meta.features = meta.features %>% left_join(
    reference %>% group_by(gene_id) %>% summarise(
      transcript_length = max(abs(end - start) + 1)
    ),
    'gene_id'
  )
  meta.features$gene_symbol = meta.features$gene_symbol %>% replace(
    grepl('\\\\', .),
    meta.features$gene_id[grepl('\\\\', .)]
  )

  meta.features$gene_display = factor(meta.features$gene_symbol) %>% recode(
    Dl='Delta',
    `E(spl)malpha-BFM`='E(spl)mα',
    `E(spl)mbeta-HLH`='E(spl)mβ',
    alphaTry='α-Try',
    betaTry='β-Try'
  ) %>% as.character
  meta.features %>% column_to_rownames('gene_symbol')
}

midgut_seurat_for_technology <- function(counts_file, metadata_file, metafeatures, technology, all_size_factors, spca_genes=NULL) {
  counts = read.csv(counts_file, row.names=1, check.names=F)
  metadata = read.csv(metadata_file, row.names=1, check.names=F)
  counts = counts[, rownames(metadata) %>% subset(metadata$technology == technology)]
  stopifnot(all.equal(rownames(counts), metafeatures$gene_id))
  rownames(counts) = rownames(metafeatures)
  technology. = technology
  metadata = metadata %>% subset(technology == technology., select = -c(UMAP_1, UMAP_2))

  seurat = CreateSeuratObject(counts, meta.data=metadata)
  seurat[['RNA']]@meta.features[, colnames(metafeatures)] = metafeatures
  seurat$pctRibo = seurat %>% PercentageFeatureSet('^Rp[SL]')
  seurat = seurat %>% NormalizeData
  spca.features = seurat %>% FindVariableFeatures(nfeatures=1000) %>% VariableFeatures
  spca.features = union(spca_genes, spca.features) %>% head(1000)
  seurat = (
    seurat
    %>% FindVariableFeatures(nfeatures=2000)
    # We decided not to regress out pctMito. ISCs are high in mito content in
    # this experiment, so removing all variance explained by pctMito definitely
    # would make it harder to form separate clusters of ISC and EB.
    %>% ScaleData
    %>% RunPCA(verb=F)
    %>% RunUMAP(dims=1:30)
  )
  # Fix orientation of the UMAP (stem cells in top-left quadrant).
  seurat[["umap"]]@cell.embeddings <- seurat[["umap"]]@cell.embeddings %*% (
    matrix(diag(c(-1, -1)), nrow=2, dimnames=list(NULL, c("UMAP_1", "UMAP_2")))
  )
  seurat@misc$covar = var(
    t(
      seurat[['RNA']]@scale.data[spca.features,]
    )
  )
  seurat$size_factor = all_size_factors[colnames(seurat)]

  seurat
}

runif_uint64 <- function() {
  # Julia will accept an UInt64 random seed.
  julia_seed = do.call(
    paste0,
    append(
      list("0x"),
      as.character(
        as.hexmode(
          runif(2, max = bitwShiftL(1, 30) * 4)
          + bitwShiftL(-1, 30) * 2)
      )
    )
  )
}

seurat_spca <- function(seurat, matrix_name, varnum, npcs, search_cap, eigen_gap, assay='RNA', do.correct.elbow = F, cgroup = NULL) {
  covar = seurat@misc[[matrix_name]]
  seurat_spca_compute_feature_loadings(covar, varnum, npcs, search_cap, eigen_gap, cgroup) %>%
    seurat_spca_from_feature_loadings(seurat, assay, do.correct.elbow)
}

seurat_spca_compute_feature_loadings <- function(
  covar, varnum, npcs, search_cap, eigen_gap, cgroup
) {
  julia_seed = runif_uint64()
  feature_loadings = run_optimal_spca(
    covar, K=varnum, D=npcs, search_cap=search_cap, eigen_gap=eigen_gap,
    uint64_seed=julia_seed, cgroup=cgroup
  )
  # Fix the signs of feature loadings using the median sign.
  feature_loadings_heatmap = feature_loadings %>% apply(1, \(v) v %>% subset(. != 0)) %>% t
  feature_loadings = feature_loadings * rowMedians(sign(feature_loadings_heatmap))
  t(feature_loadings)
}

seurat_spca_from_feature_loadings <- function(
  feature_loadings, seurat, assay, do.correct.elbow,
  do.rename.features = TRUE
) {
  if (do.rename.features)
    colnames(feature_loadings) <- paste0('SPARSE_', seq(ncol(feature_loadings)))
  dataMeans = rowMeans(
    seurat[[assay]]@data[rownames(feature_loadings), ]
  )
  dataSds = rowSds(
    seurat[[assay]]@data[rownames(feature_loadings), ]
  )
  # data should be nonnegative, so we will update sd to a pseudo-sd of "1" if
  # the mean is 0.
  dataSds = dataSds %>% replace(dataMeans == 0, 1)
  scale_data_command_name <- str_glue("ScaleData.{assay}")
  if (scale_data_command_name %in% names(seurat@commands))
    scale_data_command <- seurat@commands[[scale_data_command_name]]
  else
    scale_data_command <- list(do.center=TRUE, do.scale=TRUE)
  scale_data_from_zero = (
    seurat[[assay]]@scale.data[rownames(feature_loadings), ]
    + (
      (
        if (scale_data_command$do.center)
        dataMeans
        else 0
      ) / (
        if (scale_data_command$do.scale)
        dataSds
        else 1
      )
    )
  )
  cell_embeddings = t(scale_data_from_zero) %*% feature_loadings
  stdev = colSds(cell_embeddings)
  obj = CreateDimReducObject(
    cell_embeddings,
    feature_loadings,
    stdev=stdev,
    assay=assay,
    key='SPARSE_'
  )
  if (do.correct.elbow) {
    obj.perm = order(obj@stdev, decreasing=T)
    obj.names = colnames(obj@cell.embeddings)
    # un-permuted feature loadings and their permutation
    obj@misc$correction.perm = obj.perm
    obj@misc$search.feature.loadings = obj@feature.loadings
    obj@cell.embeddings = obj@cell.embeddings[, obj.perm, drop=FALSE]
    colnames(obj@cell.embeddings) = obj.names
    obj@feature.loadings = obj@feature.loadings[, obj.perm, drop=FALSE]
    colnames(obj@feature.loadings) = obj.names
    obj@stdev = obj@stdev[obj.perm]
    names(obj@stdev) = obj.names
  }
  obj
}

seurat_spca_from_feature_loadings_nocenter <- function(
  feature_loadings, seurat, assay, do.correct.elbow
) {
  colnames(feature_loadings) = paste0('SPARSE_', seq(ncol(feature_loadings)))
  dataMeans = rowMeans(
    seurat[[assay]]@data[rownames(feature_loadings), ]
  )
  dataSds = rowSds(
    seurat[[assay]]@data[rownames(feature_loadings), ]
  )
  # data should be nonnegative, so we will update sd to a pseudo-sd of "1" if
  # the mean is 0.
  dataSds = dataSds %>% replace(dataMeans == 0, 1)
  cell_embeddings = t(seurat[[assay]]@scale.data[rownames(feature_loadings), ]) %*% feature_loadings
  stdev = colSds(cell_embeddings)
  obj = CreateDimReducObject(
    cell_embeddings,
    feature_loadings,
    stdev=stdev,
    assay=assay,
    key='SPARSE_'
  )
  if (do.correct.elbow) {
    obj.perm = order(obj@stdev, decreasing=T)
    obj.names = colnames(obj@cell.embeddings)
    # un-permuted feature loadings
    obj@misc$search.feature.loadings = obj@feature.loadings
    obj@cell.embeddings = obj@cell.embeddings[, obj.perm]
    colnames(obj@cell.embeddings) = obj.names
    obj@feature.loadings = obj@feature.loadings[, obj.perm]
    colnames(obj@feature.loadings) = obj.names
    obj@stdev = obj@stdev[obj.perm]
    names(obj@stdev) = obj.names
  }
  obj
}

# Function on Seurat object. To be released later in an spcaFeatures library.
RunSparsePCA <- function(seurat, ...) {
  more_args <- list(...)
  if (!("matrix_name" %in% names(more_args))) {
    # User did not create a "covar" matrix so we will create it.
    scale.data <- seurat[[
      if (is.character(more_args$assay)) more_args$assay else "RNA"
    ]]@scale.data
    if ("nfeatures" %in% names(more_args)) {
      scale.data <- scale.data[head(VariableFeatures(seurat), more_args$nfeatures), ]
      more_args <- more_args[names(more_args) != "nfeatures"]
    }
    seurat@misc$covar <- tcrossprod(scale.data) / (nrow(scale.data) - 1)
    more_args <- more_args %>% append(list(matrix_name="covar"))
  }
  seurat[['spca']] = do.call(seurat_spca, append(list(seurat), more_args))
  seurat
}

spca_with_centered_umap <- function(seurat, spca, dims=1:50, seed.use=1, umap_transform=diag(2)) {
  seurat[['spca']] = spca
  seurat[['spca.centered']] <- CreateDimReducObject(
    embeddings = spca@cell.embeddings %>% scale(scale = FALSE, center = TRUE),
    key = "SPCACEN_",
    assay = "RNA"
  )
  seurat <- seurat %>%
    RunUMAP(
      reduction = "spca.centered", reduction.name = "umap.spca", dims = dims,
      seed.use = seed.use
    )
  dimnames(umap_transform) <- rep(list(colnames(seurat[['umap.spca']])), 2)
  seurat[['umap.spca']]@cell.embeddings <- (
    seurat[['umap.spca']]@cell.embeddings %*% umap_transform
  )
  seurat
}

midgut_classify_cell_types <- function(seurat, colname) {
  midgut_levels = c(
    'ISC',
    'EB',
    'dEC',
    'EC',
    'EC-like',
    'EE',
    'copper/iron',
    'LFC',
    'cardia'
  )
  midgut_subclassif_levels = c(
    'ISC',
    'EB',
    'dEC',
    'aEC',
    'EC-like',
    'mEC',
    'pEC',
    'EE',
    'copper/iron',
    'LFC',
    'cardia'
  )
  seurat@meta.data[, str_replace(colname, 'clusters', 'classif')] = (
    seurat@meta.data[, colname]
    %>% fct_relabel(\(names) str_replace(names, '\\..*|[0-9]+', ''))
    %>% fct_relevel(midgut_levels)
  )
  seurat@meta.data[, str_replace(colname, 'clusters', 'subclassif')] = (
    seurat@meta.data[, colname]
    %>% fct_relabel(\(names) str_replace(names, '(?<=EE|others)\\..*|[0-9]+', ''))
    %>% fct_recode(aEC="EC.anterior", mEC="EC.meso", pEC="EC.posterior") %>%
      fct_relevel(midgut_subclassif_levels)
  )
  seurat
}

compute_sct_clusters <- function(indrop.sct.pca) {
  indrop.sct.snn <- FindNeighbors(indrop.sct.pca@cell.embeddings[, 1:9])$snn
  sct_clusters <- FindClusters(
    indrop.sct.snn, res=1.42, random.seed=5
  ) %>%
    rownames_to_column %>%
    pull(res.1.42, rowname)
  sct_clusters %>% recode(
    `0`='EC-like1',
    `1`='EC-like2',
    `2`='dEC1',
    `3`='copper/iron',
    `4`='ISC',
    `5`='dEC2',
    `6`='EE1',
    `7`='aEC',
    `8`='pEC',
    `9`='aEC2',
    `10`='others.1',
    `11`='aEC3',
    `12`='EE2',
    `13`='aEC4',
    `14`='cardia',
    `15`='aEC5',
    `16`='others.2',
    `17`='pEC2',
    `18`='LFC',
    `19`='aEC6',
    `20`='aEC7',
    `21`='aEC8',
    `22`='EB',
    `23`='mEC',
    `25`='EC-like3',
    `26`='EC-like4',
    `27`='pEC3',
    `28`='pEC4',
    `29`='pEC5',
    `30`='others.3',
    `31`='cardia.2',
    `32`='mEC2'
  )
}

build_midgut_misc_stats <- function(indrop, indrop.sct.pca, npcs=25) {
  stem.explained = mapply(
    \(ident, reduction) {
      ident = indrop@meta.data[, ident]
      stem = reduction@cell.embeddings %>% subset(ident %in% c('ISC','EB'))
      ISC = reduction@cell.embeddings %>% subset(ident == 'ISC')
      EB = reduction@cell.embeddings %>% subset(ident == 'EB')
      result = c(
        ISC=sum(colVars(ISC[, 1:npcs])) * nrow(ISC) / nrow(stem),
        EB=sum(colVars(EB[, 1:npcs])) * nrow(EB) / nrow(stem)
      )
      c(
        ISC.sample=sum(colVars(ISC[, 1:npcs])),
        EB.sample=sum(colVars(EB[, 1:npcs])),
        result,
        SSB=sum(colVars(stem[, 1:npcs])) - sum(result)
      )
    },
    c('pca_clusters', 'spca_clusters', 'sct_clusters'),
    list(indrop[['pca']], indrop[['spca']], indrop.sct.pca)
  )
  pca.to.spca.shrink.var = (
    stem.explained[,'pca_clusters'] / stem.explained[,'spca_clusters']
  )
  list(
    stem.explained=stem.explained,
    pca.to.spca.shrink.var=pca.to.spca.shrink.var
  )
}

# Replace a singular vector with the negation, etc, for illustrative purposes.
# Will be used for cosmetic tweaks to the non-sparse PCA model.
replace_pca_embedding_feature <- function(seurat, column_name, fn) {
  seurat[['pca']]@cell.embeddings[, column_name] <- seurat[[
    'pca'
  ]]@cell.embeddings[, column_name] %>% fn
  seurat
}

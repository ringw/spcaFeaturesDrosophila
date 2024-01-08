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
  metadata = metadata %>% subset(technology == technology.)

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
  seurat@misc$covar = var(
    t(
      seurat[['RNA']]@scale.data[spca.features,]
    )
  )
  seurat$size_factor = all_size_factors[colnames(seurat)]

  seurat
}

seurat_spca <- function(seurat, matrix_name, varnum, npcs, search_cap, eigen_gap, assay='RNA', do.correct.elbow = F) {
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
  covar = seurat@misc[[matrix_name]]
  feature_loadings = run_optimal_spca(covar, K=varnum, D=npcs, search_cap=search_cap, eigen_gap=eigen_gap, uint64_seed=julia_seed)
  # Fix the signs of feature loadings using the median sign.
  feature_loadings_heatmap = feature_loadings %>% apply(1, \(v) v %>% subset(. != 0)) %>% t
  feature_loadings = feature_loadings * rowMedians(sign(feature_loadings_heatmap))
  feature_loadings = t(feature_loadings)
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
  scale_data_from_zero = (
    seurat[['RNA']]@scale.data[rownames(feature_loadings), ]
    + (
      (
        if (seurat@commands[[
          paste('ScaleData', assay, sep='.')
        ]]$do.center)
        dataMeans
        else 0
      ) / (
        if (seurat@commands[[
          paste('ScaleData', assay, sep='.')
        ]]$do.scale)
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
    obj.names = paste0(
      'SPARSE_',
      seq(ncol(feature_loadings)),
      ifelse(obj.perm == seq_along(obj.perm), '\'', '')
    )
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
  seurat[['spca']] = seurat_spca(seurat, ...)
  seurat
}

spca_with_centered_umap <- function(seurat, spca, dims=1:50, seed.use=1) {
  seurat[['spca']] = spca
  seurat[['umap.spca']] = RunUMAP(
    spca@cell.embeddings[, dims] %>% scale(scale=F, center=T), seed.use=seed.use
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
  seurat@meta.data[, str_replace(colname, 'clusters', 'classif')] = (
    seurat@meta.data[, colname]
    %>% fct_relabel(\(names) str_replace(names, '\\..*|[0-9]+', ''))
    %>% fct_relevel(midgut_levels)
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
    `2`='EC1',
    `3`='copper/iron',
    `4`='ISC',
    `5`='dEC1',
    `6`='EE1',
    `7`='EC2',
    `8`='EC3',
    `9`='EC4',
    `10`='others.1',
    `11`='EC5',
    `12`='EE2',
    `13`='EC6',
    `14`='cardia',
    `15`='EC7',
    `17`='EC8',
    `18`='LFC',
    `19`='EC9',
    `20`='EC10',
    `21`='EC11',
    `22`='EB',
    `23`='EC',
    `25`='EC-like3',
    `26`='EC-like4',
    `27`='EC12',
    `28`='EC13',
    `29`='EC14',
    `32`='EC15'
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

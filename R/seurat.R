midgut_seurat_for_technology <- function(counts_file, metadata_file, gtf_file, technology) {
  counts = read.csv(counts_file, row.names=1, check.names=F)
  metadata = read.csv(metadata_file, row.names=1, check.names=F)
  counts = counts[, rownames(metadata) %>% subset(metadata$technology == technology)]
  technology. = technology
  metadata = metadata %>% subset(technology == technology.)

  # Now map counts rownames using the gtf file.
  reference = read.table(
    gtf_file,
    sep = '\t',
    col.names = c('chr', 'source', 'type', 'start', 'end', 'sc', 'strand', 'fr', 'annotation'),
    header = F,
    quote = ''
  ) %>% subset(grepl('RNA', type))
  reference$gene_id = str_extract(reference$annotation, 'gene_id "([^"]+)"', group=1)
  reference$gene_symbol = str_extract(reference$annotation, 'gene_symbol "([^"]+)"', group=1)
  meta.features = reference %>% group_by(gene_id) %>% summarise(
    gene_symbol = max(gene_symbol),
    transcript_length = max(abs(end - start) + 1)
  )
  meta.features = data.frame(gene_id=rownames(counts)) %>% left_join(
    meta.features, 'gene_id'
  )
  rownames(meta.features) = meta.features$gene_symbol %>% replace(
    duplicated(.) | is.na(.),
    meta.features$gene_id[duplicated(.) | is.na(.)]
  )
  rownames(counts) = rownames(meta.features)

  seurat = CreateSeuratObject(counts, meta.data=metadata)
  seurat[['RNA']]@meta.features[,c('gene_id','gene_symbol','transcript_length')] = meta.features
  seurat$pctRibo = seurat %>% PercentageFeatureSet('^Rp[SL]')
  seurat = seurat %>% NormalizeData
  spca.features = seurat %>% FindVariableFeatures(nfeatures=1000) %>% VariableFeatures
  seurat = (
    seurat
    %>% FindVariableFeatures(nfeatures=2000)
    %>% ScaleData(vars.to.regress = 'pctMito')
    %>% RunPCA(verb=F)
    %>% RunUMAP(dims=1:30)
  )
  seurat@misc$covar = var(
    t(
      seurat[['RNA']]@scale.data[spca.features,]
    )
  )
  
  seurat
}

seurat_spca <- function(seurat, matrix_name, varnum, npcs, search_cap, assay='RNA') {
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
  feature_loadings = run_optimal_spca(covar, K=varnum, D=npcs, search_cap=search_cap, uint64_seed=julia_seed)
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
  obj
}

# Function on Seurat object. To be released later in an spcaFeatures library.
RunSparsePCA <- function(seurat, matrix_name, varnum, npcs, search_cap, assay='RNA') {
  seurat[['spca']] = seurat_spca(seurat, matrix_name, varnum, npcs, search_cap, assay=assay)
  seurat
}

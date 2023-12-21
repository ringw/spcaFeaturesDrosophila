add_pooled_size_factors <- function(
    counts_file, metadata_file, output_path='GSE120537_augmented.csv') {
  counts = read.csv(counts_file, row.names=1, check.names=F)
  metadata = read.csv(metadata_file, row.names=1, check.names=F)
  clusters = quickCluster(counts, block=metadata$technology)
  metadata$size = pooledSizeFactors(counts, clusters=clusters)
  write.csv(metadata, output_path, row.names=T)
  output_path
}

midgut_pooled_size_factors <- function(counts_file, metadata_file) {
  counts = read.csv(counts_file, row.names=1, check.names=F)
  metadata = read.csv(metadata_file, row.names=1, check.names=F)
  clusters = quickCluster(counts, block=metadata$technology, min.size=300)
  setNames(
    pooledSizeFactors(counts, clusters=clusters),
    colnames(counts))
}
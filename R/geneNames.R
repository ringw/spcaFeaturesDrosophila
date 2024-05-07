gene_names <- list(
  # Symbol change to Delta from latest Flybase release
  Dl='Delta',
  `E(spl)malpha-BFM`="E(spl)m\u03B1",
  `E(spl)mbeta-HLH`="E(spl)m\u03B2",
  `E(spl)m3-HLH`="E(spl)m3",
  alphaTry="\u03B1Try",
  betaTry="\u03B2Try",
  lrRNA="mt:lrRNA",
  # From latest Flybase release
  CG2556="mdu"
)

feature_names <- gene_names %>%
  append(
    paste0('SPARSE_', 1:50) %>%
      sapply(\(n) n %>% str_replace('SPARSE_', 'SPC'), simplify=F)
  )

display_gene_names <- function(fct) {
  do.call(
    recode,
    append(
      list(fct),
      feature_names
    )
  )
}
gene_names <- list(
  Dl='Delta',
  `E(spl)malpha-BFM`="E(spl)m\u03B1",
  `E(spl)mbeta-HLH`="E(spl)m\u03B2",
  betaTry="\u03B2Try"
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
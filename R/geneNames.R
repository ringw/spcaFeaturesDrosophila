gene_names <- list(
  Dl='Delta',
  `E(spl)malpha-BFM`="E(spl)m\u03B1",
  `E(spl)mbeta-HLH`="E(spl)m\u03B2"
)

display_gene_names <- function(fct) {
  do.call(
    recode,
    append(
      list(fct),
      gene_names
    )
  )
}
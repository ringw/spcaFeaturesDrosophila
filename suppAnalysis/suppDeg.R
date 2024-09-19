display_pca_spca_deg_data <- function(deg_data, midgut.metafeatures, ...) {
  pca <- display_deg_data(deg_data$pca_clusters)
  spca <- display_deg_data(deg_data$spca_clusters)
  tbl <- full_join(
    summarise(
      pca, rowname, nExpressed_ISCEB_PCA = total_num_expressed, Mu_ISC_PCA = baseMean,
      `L2FC95-PCA`=log2FoldChange - qnorm(0.975) * lfcSE,
      L2FC_PCA=log2FoldChange,
      `L2FC95+PCA`=log2FoldChange + qnorm(0.975) * lfcSE,
      FDR_PCA=q_val,
    ),
    summarise(
      spca, rowname, nExpressed_ISCEB_SPCA = total_num_expressed,
      Mu_ISC_SPCA = baseMean,
      `L2FC95-SPCA`=log2FoldChange - qnorm(0.975) * lfcSE,
      L2FC_SPCA=log2FoldChange,
      `L2FC95+SPCA`=log2FoldChange + qnorm(0.975) * lfcSE,
      FDR_SPCA=q_val,
    ),
    "rowname"
  )
  tbl <- tibble(
    symbol_r6.27 = pull(tbl, 1),
    flybase_r6.27 = midgut.metafeatures[pull(tbl, 1), "gene_id"],
    tbl[-1]
  )
  tbl[(tbl$FDR_PCA >= 0.05) %>% replace(is.na(.), FALSE), 3:7] <- NA
  tbl[(tbl$FDR_SPCA >= 0.05) %>% replace(is.na(.), FALSE), 9:13] <- NA
  tbl[
    (tbl[3:7] %>% as.matrix %>% is.na %>% `!`() %>% rowAlls) |
      (tbl[9:13] %>% as.matrix %>% is.na %>% `!`() %>% rowAlls),
  ] %>%
    arrange(
      ifelse(
        (sign(`L2FC95-SPCA`) == sign(`L2FC95+SPCA`)) %>%
          replace(is.na(.), 0),
        pmin(abs(`L2FC95-SPCA`), abs(`L2FC95+SPCA`)) * sign(L2FC_SPCA),
        0
      ),
      ifelse(
        sign(`L2FC95-PCA`) == sign(`L2FC95+PCA`),
        pmin(abs(`L2FC95-PCA`), abs(`L2FC95+PCA`)) * sign(L2FC_PCA),
        0
      )
    )
}

tar_load(indrop.deg | tenx.deg | midgut.metafeatures)
print(display_pca_spca_deg_data(indrop.deg, midgut.metafeatures))
print(display_pca_spca_deg_data(tenx.deg, midgut.metafeatures))

build_glms <- function(seurat, all_sizes, columns_to_fit) {
  sapply(
    columns_to_fit,
    \(column_name) {
      colData = data.frame(
        cluster = seurat@meta.data[, column_name]
      )
      g = glm_gp(
        seurat[['RNA']]@counts,
        ~ cluster,
        colData,
        size_factors = seurat$size_factor,
        on_disk = F,
        verbose = T
      )
      # Use sparseMatrix for counts, not matrix.
      assay(g$data) = seurat[['RNA']]@counts
      g
    },
    simplify=F
  )
}

build_glm_cpm <- function(seurat, column_to_fit, dispersion_trend) {
  size_factors <- FetchData(seurat, c("size_factor", "nCount_RNA", column_to_fit)) %>%
    rownames_to_column %>%
    list %>%
    append(rlang::syms(column_to_fit)) %>%
    do.call(group_by, .) %>%
    mutate(
      # Normalization factor: Take all cells in the cluster (the group) - the
      # estimate for these cells coming from the summed size factor. Add 1MM
      # multiplier. Normalize out the sum of UMI counts for the cells.
      cpm_scalar = sum(size_factor) * 1000 * 1000 / sum(nCount_RNA),
      # Take away the normalization factor (divide out) when building the model.
      # Then the coefficients will be ln-CPM (log link function).
      cpm_size_factor = size_factor / cpm_scalar
    )
  g <- glm_gp(
    seurat[['RNA']]@counts,
    ~ 0 + ident,
    seurat@meta.data %>% mutate(ident = FetchData(seurat, column_to_fit)[, 1]),
    size_factors = size_factors$cpm_size_factor,
    overdispersion = dispersion_trend,
    overdispersion_shrinkage = F,
    on_disk = F,
    verbose = T
  )
  # Use sparseMatrix for counts, not matrix.
  assay(g$data) = seurat[['RNA']]@counts
  g
}

predict_glm_cpm <- function(glm.cpm) {
  cpm <- predict(
    glm.cpm,
    diag(ncol(glm.cpm$Beta)) %>%
      matrix(
        nrow=nrow(.),
        ncol=ncol(.),
        dimnames=rep(list(colnames(glm.cpm$Beta)), 2)
      ),
    offset=0,
    type="response"
  )
  colnames(cpm) <- colnames(cpm) %>% str_replace("^ident", "")
  cpm
}

build_present_gene_list <- function(seurat, columns) {
  grid = expand.grid(column=columns, level=c('ISC','EB'))
  grid$level = as.character(grid$level)
  # Simplify this. For our two groups, if each group is sparse with at most one
  # count in the regression model, then the log-likelihood seems to be giving an
  # error. We can keep genes with no counts in one class, but at a minimum we
  # should quantify the gene in 2 cells.
  setNames(
    rowMins(sapply(
      columns,
      \(column_name) as.data.frame(
        t(seurat[['RNA']]@counts)
      ) %>% split(seurat@meta.data[,column_name]) %>% subset(names(.) %in% c('ISC','EB'))
      %>% sapply(\(counts) colSums(counts != 0)) %>% rowMaxs
    )) >= 2,
    rownames(seurat))
}

build_de_data <- function(glms, present_genes) {
  mapply(
    \(column_name, model) {
      keep_cells = colData(model$data)[, 'cluster'] %in% c('ISC','EB')
      glm.binary = glm_gp(
        model$data[present_genes, keep_cells],
        design = ~ cluster,
        size_factors = model$size_factors[keep_cells],
        overdispersion = model$overdispersion_shrinkage_list$dispersion_trend[present_genes],
        overdispersion_shrinkage = F,
        on_disk = F
      )
      message('Predict s.e. of log fold change')
      pred = with(
        predict(glm.binary, matrix(c(0,1), nrow=1), offset=0, se.fit=T),
        cbind(fit, se.fit)
      )
      pred[!is.finite(pred[,2]), 1] = NA
      message('Shrink log fold change')
      a = apeglm(
        glm.binary$data,
        glm.binary$model_matrix,
        log.lik=NULL,
        coef=2,
        param=glm.binary$overdispersions,
        mle=pred,
        method='nbinomCR',
        offset=glm.binary$Offset
      )
      message('Perform F-test of original model')
      a$mle.test = test_de(model, 'clusterEB')
      a
    },
    names(glms),
    glms,
    SIMPLIFY = F
  )
}

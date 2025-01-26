simple_violin_plot <- function(data, size=0.6, ...) {
  colnames(data) <- c("cluster", "embedding", "weight")
  data$cols <- "black"
  data$cols[c(1:9, 37:45)] <- "#990000"
  data <- data %>% arrange(desc(cols))
  (
    ggplot(data, aes(cluster, embedding, fill=cluster))
    + geom_violin(aes(weight=weight), linewidth=0.5, ...)
    + geom_jitter(aes(color=cols), shape=16, stroke=NA, size=size) +
    scale_color_identity()
  ) + theme_cowplot() + scale_fill_hue() +
    scale_y_continuous(
      expand = c(0, 0)
    ) + labs(
      x = NULL, y = NULL
    ) + theme(
      plot.margin = margin(1, 1, 1, 1), aspect.ratio = 1
    )
}

plot_logumi_celltypes <- function(indrop) {
  red_names <- c(
    pca_classif = "LN PCA",
    sct_classif = "SCT PCA",
    spca_classif = "LN SPCA"
  )
  indrop$logUMI <- log(indrop$nUMI) / log(10)
  gg <- sapply(
    names(red_names),
    \(n) tiny_violin_plot(
      indrop,
      n,
      c("EB", "ISC", "EC", "EC-like", "EE"),
      "logUMI",
      pt.size = 0.1
    ) +
      geom_boxplot(outlier.shape = NA, fill = "transparent") +
      scale_x_discrete(
        red_names[n]
      ) +
      scale_y_continuous(
        bquote(log[10]*"(nUMI)"), breaks = seq(2.5, 4.5, by = 0.5),
        labels = c("", "3", "", "4", "")
      ),
    simplify=FALSE
  )
  gg[[2]] <- gg[[2]] + scale_y_continuous(NULL, labels=NULL)
  gg[[3]] <- gg[[3]] + scale_y_continuous(NULL, labels=NULL)
  do.call(
    cbind,
    sapply(
      gg,
      \(g) g %>% set_panel_size(w = unit(1.5, "in"), h = unit(0.75, "in"))
    )
  )
}

plot_pctMito_celltypes <- function(indrop) {
  red_names <- c(
    pca_classif = "LN PCA",
    sct_classif = "SCT PCA",
    spca_classif = "LN SPCA"
  )
  gg <- sapply(
    names(red_names),
    \(n) tiny_violin_plot(
      indrop,
      n,
      c("EB", "ISC", "EC", "EC-like", "EE"),
      "pctMito",
      pt.size = 0.1
    ) +
      geom_boxplot(outlier.shape = NA, fill = "transparent") +
      scale_x_discrete(
        red_names[n]
      ) +
      scale_y_continuous(
        "Mito Fraction", labels=\(v) percent(v/100)
      ),
    simplify=FALSE
  )
  gg[[2]] <- gg[[2]] + scale_y_continuous(NULL, labels=NULL)
  gg[[3]] <- gg[[3]] + scale_y_continuous(NULL, labels=NULL)
  do.call(
    cbind,
    sapply(
      gg,
      \(g) g %>% set_panel_size(w = unit(1.5, "in"), h = unit(0.75, "in"))
    )
  )
}

plot_logumi_corr <- function(indrop, indrop.sct.pca, var.test) {
  logUMI <- log(indrop$nUMI) / log(10)
  pctMito <- indrop$pctMito
  pca.var <- indrop[["pca"]]@stdev^2
  sct.var <- indrop.sct.pca@stdev^2
  spca.var <- indrop[["spca"]]@stdev^2
  npca <- 36
  nspca <- 36
  ylim <- c(0, 0.2)
  rbind(
    tibble(x="LN PCA", v=as.numeric(cor(FetchData(indrop, str_glue("PC_{1:npca}")), get(var.test)))^2, w=head(pca.var, npca)),
    tibble(x="SCT PCA", v=as.numeric(cor(indrop.sct.pca@cell.embeddings[, 1:npca], get(var.test)))^2, w=head(sct.var, npca)),
    tibble(x="LN SPCA", v=as.numeric(cor(FetchData(indrop, str_glue("SPARSE_{1:nspca}")), get(var.test)))^2, w=head(spca.var, nspca))
  ) %>%
    mutate(
      x = x %>% factor(unique(.))
    ) %>%
    simple_violin_plot(size = 2, bw = 0.02) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    # scale_y_continuous(trans="sqrt") +
    scale_fill_hue(h = c(270, 450) + 15) +
    labs(y = bquote(R^2)) +
    theme(
      axis.text = element_text(size = unit(10, "pt")),
      axis.text.x = element_text(angle = 45, hjust = 0.8),
      axis.title = element_text(size = unit(12, "pt")),
      legend.position = "none"
    )
}

plot_logumi_bar <- function(indrop.expl.var.stats.logumi) {
  red_names <- c(
    pca="LN PCA",
    pca.sct="SCT PCA",
    spca2pca="LN SPCA"
  )
  indrop.expl.var.stats.logumi <- indrop.expl.var.stats.logumi %>%
    subset(
      ident %in% c("stem-like", "EC", "EE") &
        reduction %in% names(red_names)
    )
  # ylim <- c(0, 1.02 * max(indrop.expl.var.stats.logumi$pct.expl))
  ylim <- c(0, 0.2)
  gg <- sapply(
    names(red_names),
    \(n) indrop.expl.var.stats.logumi %>%
      subset(ident %in% c("stem-like", "EC", "EE") & reduction == n) %>%
      tibble(
        group = ident %>% factor(., .) %>% fct_recode(ISCEB="stem-like")
      ) %>%
      ggplot(aes(group, pct.expl, fill = group)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values = midgut.colors) +
      scale_x_discrete(
        # labels = c("ISCEB  ", "EC", "EE")
      ) +
      scale_y_continuous(limits = ylim, expand = c(0, 0)) +
      coord_cartesian(
        xlim = c(1, 3) + c(-1, 1) * 0.45
      ) +
      labs(
        x = red_names[n],
        y = bquote("Multiple"*" "*R^2),
      ) +
      theme_cowplot() +
      theme(
        axis.text = element_text(size = unit(10, "pt")),
        axis.text.x = element_text(angle = 45, hjust = 0.8),
        axis.title = element_text(size = unit(12, "pt")),
        legend.position = "none",
        plot.margin = margin(1, 0.2 * 72 / 25.4, 0, 0.2 * 72 / 25.4)
      ),
    simplify=FALSE
  )
  for (i in 2:3) {
    gg[[i]] <- gg[[i]] +
      scale_y_continuous(
        name = NULL,
        breaks = NULL,
        limits = ylim,
        expand = c(0, 0)
      ) +
      theme(
        axis.line.y = element_blank()
      )
  }
  do.call(
    cbind,
    sapply(
      gg,
      \(g) g %>% set_panel_size(w = unit(0.75, "in"), h = unit(1.5, "in")),
      simplify = FALSE
    )
  )
}

midgut_explained_variance_dims = c(pca=25, pca.sct=25, pca.subset=9, spca=34, spca2pca=25)
midgut_expl_var_stem_like <- function(
  indrop, indrop.sct.pca,
  variable = log(indrop$nCount_RNA)
) {
  # Create the "coarse" factor levels (combine ISC+EB) for ANOVA.
  for (n in c("pca", "spca", "sct"))
    indrop[[str_glue("{n}_coarse")]] <- indrop[[str_glue("{n}_classif")]] %>%
      deframe %>%
      recode(ISC='stem-like', EB='stem-like', `EC-like`='EC')

  spca_dims <- indrop@commands$FindNeighbors.RNA.spca$dims
  indrop[["SCT"]] <- CreateAssayObject(counts = GetAssayData(indrop, "counts"))
  indrop[["pca.sct"]] <- indrop.sct.pca

  indrop <- indrop %>% RunPCA(
    reduction.key = "PCSUBSET_",
    reduction.name = "pca.subset",
    features = indrop[["spca"]]@feature.loadings[
      ,
      spca_dims
    ] %>%
      rowSums %>%
      as.logical %>%
      which %>%
      names,
    verb = FALSE,
    approx = FALSE
  )
  indrop[["spca2pca"]] <- RunPCA(
    indrop[['spca']]@cell.embeddings[, spca_dims] %>% t,
    assay = "RNA", reduction.key = "SPARSE2DENSE_",
    verb = FALSE,
    approx = FALSE
  )
  expl.var.clusters = c(
    pca='pca_coarse',
    pca.sct='sct_coarse',
    pca.subset = 'pca_coarse',
    spca='spca_coarse',
    spca2pca='spca_coarse'
  )
  indrop.pct.var.logUMI = expand.grid(
    ident = levels(indrop$pca_coarse),
    reduction = c('pca', 'pca.sct', 'pca.subset', 'spca', 'spca2pca')
  ) %>% as.matrix %>% as.data.frame %>%
    rowwise %>%
    mutate(
      mn = manova(
        spca ~ nuisance_variable,
        list(
          spca = indrop@reductions[[reduction]]@cell.embeddings[, 1:midgut_explained_variance_dims[reduction]],
          nuisance_variable = variable
        ),
        subset = FetchData(indrop, expl.var.clusters[as.character(reduction)]) %>%
          pull(1) %>%
          `==`(as.character(ident))
      ) %>%
        list
    ) %>%
    mutate(
      residuals = sum(colVars(mn$residuals)),
      fitted.values = sum(colVars(mn$fitted.values)),
      pct.expl = fitted.values / (residuals + fitted.values),
      .keep = "unused"
    )
}

plot_midgut_expl_var <- function(
  expl_var_stem_like,
  reduction = "pca",
  ident_levels = c("stem-like", "EC", "EE")
) {
  reduction. <- reduction
  (
    expl_var_stem_like %>% subset(ident %in% ident_levels & reduction == reduction.) %>%
    within(ident <- ident %>% fct_relevel(ident_levels))
    %>% ggplot(aes(x=ident, fill=ident, y=pct.expl))
    + geom_bar(stat='identity')
    + coord_flip()
    + scale_x_discrete(limits=rev)
    + scale_y_continuous(expand=c(0,0), labels=percent, breaks=scales::breaks_pretty(3), limits=c(0,0.25))
    + scale_fill_manual(
      values = c(midgut.colors, `stem-like`="goldenrod"),
      guide = guide_none()
    )
    + labs(x=NULL, y=NULL)
  )
}

midgut_expl_var_stem_between <- function(
  indrop,
  variable = log(indrop$nCount_RNA),
  legend.position = "right"
) {
  mn.pca = manova(
    pca ~ pca_clusters,
    append(list(pca=indrop[['pca']]@cell.embeddings[,1:midgut_explained_variance_dims["pca"]]), indrop@meta.data),
    subset=indrop$pca_clusters %in% c('ISC','EB')
  )
  mn.pca.between = manova(
    fitted.values ~ variable.,
    data = mn.pca["fitted.values"] %>%
      append(list(variable. = variable[rownames(.$fitted.values)]))
  )
  mn.pca.isc = manova(
    residuals ~ variable.,
    data = mn.pca["residuals"] %>%
      append(list(variable. = variable[rownames(.$residuals)])),
    subset = (indrop$pca_clusters %>% subset(. %in% c('ISC','EB'))) == 'ISC'
  )
  mn.pca.eb = manova(
    residuals ~ variable.,
    data = mn.pca["residuals"] %>%
      append(list(variable. = variable[rownames(.$residuals)])),
    subset = (indrop$pca_clusters %>% subset(. %in% c('ISC','EB'))) == 'EB'
  )

  mn.spca = manova(
    spca ~ spca_clusters,
    append(list(spca=indrop[['spca']]@cell.embeddings[,1:midgut_explained_variance_dims["spca"]]), indrop@meta.data),
    subset=indrop$spca_clusters %in% c('ISC','EB')
  )
  mn.spca.between = manova(
    fitted.values ~ variable.,
    data = mn.spca["fitted.values"] %>%
      append(list(variable. = variable[rownames(.$fitted.values)]))
  )
  mn.spca.isc = manova(
    residuals ~ variable.,
    data = mn.spca["residuals"] %>%
      append(list(variable. = variable[rownames(.$residuals)])),
    subset = (indrop$spca_clusters %>% subset(. %in% c('ISC','EB'))) == 'ISC'
  )
  mn.spca.eb = manova(
    residuals ~ variable.,
    data = mn.spca["residuals"] %>%
      append(list(variable. = variable[rownames(.$residuals)])),
    subset = (indrop$spca_clusters %>% subset(. %in% c('ISC','EB'))) == 'EB'
  )
  sum_explained_variance <- \(reducedDim) sum(colVars(reducedDim))
  stem.explained = cbind(
    data.frame(
      reduction = rep(c('LN-PCA','LN-SPCA'), each=4),
      data = c('total','between','ISC','EB') %>% factor(., ordered=T),
      weight = c(
        1,
        1,
        mean('ISC' == (indrop$pca_clusters %>% subset(. %in% c('ISC','EB')))),
        mean('EB' == (indrop$pca_clusters %>% subset(. %in% c('ISC','EB')))),
        1,
        1,
        mean('ISC' == (indrop$spca_clusters %>% subset(. %in% c('ISC','EB')))),
        mean('EB' == (indrop$spca_clusters %>% subset(. %in% c('ISC','EB'))))
      )
    ),
    sapply(
      list(mn.pca,mn.pca.between,mn.pca.isc,mn.pca.eb,mn.spca,mn.spca.between,mn.spca.isc,mn.spca.eb),
      \(mn) mn %>% subset(grepl('residuals|fitted.values', names(.))) %>% sapply(sum_explained_variance)
    ) %>% t %>% as.data.frame
  )
  stem.graph = stem.explained[stem.explained$data != 'total', ] %>% subset(select=-weight)
  stem.graph[,3:4] = stem.graph[,3:4] * stem.explained[stem.explained$data != 'total', ] %>% pull(weight)
  stem.graph = stem.graph %>% melt(id.vars=c('reduction','data'))
  stem.graph$xmin = stem.graph %>% with(
    ifelse(reduction == 'LN-PCA', 1, 2)
    + ifelse(variable == 'fitted.values', -0.25, -0.45)
  )
  stem.graph$xmax = stem.graph %>% with(
    ifelse(reduction == 'LN-PCA', 1, 2) + 0.45
  )
  stem.graph$height = stem.graph$value
  stem.graph[stem.graph$variable == 'residuals', 'height'] = stem.graph[stem.graph$variable == 'residuals', 'height'] + stem.graph[stem.graph$variable == 'fitted.values', 'height']
  stem.graph$ymin = NA
  for (ind in seq(nrow(stem.graph))) {
    stem.graph[ind, 'ymin'] = sum(
      stem.graph %>% subset(
        variable == 'residuals'
        & reduction == stem.graph[ind, 'reduction']
        & data < stem.graph[ind, 'data']
      ) %>% pull(height)
    )
  }
  stem.graph$ymax = stem.graph$ymin + stem.graph$height
  stem.graph$data = factor(
    stem.graph$data, c(levels(stem.graph$data), 'logUMI')
  )
  stem.graph$data[stem.graph$variable == 'fitted.values'] = 'logUMI'

  ymax_detail <- max(stem.graph %>% subset(reduction == "LN-SPCA") %>% pull(ymax)) * 1.05
  ymax_detail <- 15
  ggarrange(
    stem.graph
    %>%
      ggplot
    +
      scale_x_continuous(
        breaks=c(1,2),
        labels=c('PCA', 'SPCA')
      )
    + scale_y_continuous(
      # mimic "expand" the top but not the y-intercept of the plot
      limits=c(0, max(stem.graph$ymax) * 1.05),
      expand=c(0,0)
    )
    + geom_rect_pattern(
      aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill=data,
          pattern_fill=data),
      pattern_color=NA,
      pattern_angle=45,
      pattern_spacing=0.25,
      pattern_density=0.5,
      pattern_key_scale_factor=0.1,
      show.legend = F
    )
    # Re-draw the rect because it lies over the patterned rect.
    + geom_rect(
      aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill=data),
      \(df) df %>% subset(data == 'logUMI'),
      show.legend = F
    )
    + scale_fill_manual(
      values=c(between=midgut.col$ISC, EB=midgut.col$EB, ISC=midgut.col$ISC, logUMI=hsv(0.42, 0.55, 0.72))
    )
    + scale_pattern_fill_manual(
      values=c(between=midgut.col$EB, EB=midgut.col$EB, ISC=midgut.col$ISC, logUMI=hsv(0.42, 0.55, 0.72))
    )
    # + geom_rect(
    #   aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill=data),
    #   \(df) df %>% subset(data != 'between')
    # )
    + geom_rect(
      aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
      data.frame(xmin=1.5,xmax=2.5,ymin=0,ymax=ymax_detail),
      color='black',
      fill='transparent'
    )
    + theme(
      legend.position = legend.position
    )
    ,
    stem.graph
    %>% subset(xmin >= 1.5)
    %>%
      ggplot
    +
      scale_x_continuous(
        breaks=c(1,2),
        labels=c('PCA', 'SPCA')
      )
    + scale_y_continuous(expand=c(0,0), limits=c(0,ymax_detail))
    + geom_rect_pattern(
      aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill=data,
          pattern_fill=data),
      pattern_color=NA,
      pattern_angle=45,
      pattern_spacing=0.25,
      pattern_density=0.5,
      pattern_key_scale_factor=0.1,
      pattern_yoffset=0.15
    )
    # Re-draw the rect because it lies over the patterned rect.
    + geom_rect(
      aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill=data),
      \(df) df %>% subset(data == 'logUMI'),
      show.legend = F
    )
    + scale_fill_manual(
      values=c(between=midgut.col$ISC, EB=midgut.col$EB, ISC=midgut.col$ISC, logUMI=hsv(0.42, 0.55, 0.72)),
      guide=guide_legend(title='Var. Source')
    )
    + scale_pattern_fill_manual(
      values=c(between=midgut.col$EB, EB=midgut.col$EB, ISC=midgut.col$ISC, logUMI=hsv(0.42, 0.55, 0.72)),
      guide=guide_legend(title='Var. Source')
    )
    + theme(
      legend.position = legend.position,
      plot.margin = list(right = NULL, bottom = margin(5.5, 100, 5.5, 5.5))[[legend.position]]
      ,
      legend.title = list(
        right=NULL,
        bottom=element_text(size=rel(0.8))
      )[[legend.position]]
    )
    ,
    widths=list(right = c(1,1.55), bottom = c(1.35,1))[[legend.position]]
  )
}
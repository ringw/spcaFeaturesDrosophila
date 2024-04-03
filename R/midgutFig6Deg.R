display_deg_data <- function(deg_data, min_num_expressed = 25) {
  deg_data %>%
    with(
      tibble(
        rowname = rownames(map),
        displayName = display_gene_names(rowname),
        baseMean = exp(map[, 1]),
        log2FoldChange = map[, 2] / log(2),
        lfcSE = sd[, 2] / log(2),
        p_val = mle.test$pval[match(rownames(map), mle.test$name)],
        q_val = mle.test$adj_pval[match(rownames(map), mle.test$name)],
        total_num_expressed = total_num_expressed
      )
    ) %>%
    subset(total_num_expressed >= min_num_expressed) %>%
    mutate(q_val = p.adjust(p_val, "BH"))
}

plot_arrange_deg_segments <- function(deg_models) {
  results <- deg_models %>% sapply(display_deg_data, simplify=FALSE)
  panels <- full_join(
    enframe(results, "model", "values"),
    tibble(sgn = c(-1, 1)),
    by=character()
  ) %>%
    rowwise %>%
    mutate(
      pl = (
        plot_midgut_de_panel(values, sgn)
        + labs(
          # unicode subscript 2
          x = "log\u2082(EB/ISC) fold change posterior 95% CI",
          y = if (sgn == -1) "Ranked Genes" else NULL
        )
      ) %>%
        list
    ) %>%
    group_by(sgn) %>%
    mutate(
      limits = c(
        0,
        layer_scales(pl[[1]])$x$range$range,
        layer_scales(pl[[2]])$x$range$range
      ) %>%
        range %>%
        list
    ) %>%
    rowwise %>%
    mutate(
      pl = (
        pl + expand_limits(x=limits[[1]])
      ) %>%
        list
    )
  panels$pl[[1]] <- panels$pl[[1]] + labs(tag = "A")
  panels$pl[[2]] <- panels$pl[[2]] + labs(tag = " ")
  panels$pl[[3]] <- panels$pl[[3]] + labs(tag = "B")
  panels$pl[[4]] <- panels$pl[[4]] + labs(tag = " ")
  ggarrange(
    plotlist = panels$pl,
    common.legend = T,
    legend = "bottom"
  )
}

de.labels = factor(
  c(
    'Amy-p'='Tissue (enzyme)',
    'MtnC'='Tissue (enzyme)',
    'Hr4'='Tissue (endocrine)',
    'Dl'='Stem cell',
    'His3.3A'='Stem cell',
    'Larp4B'='Stem cell',
    'trol'='Cells during development (various)',
    'Debcl'='Cells during development (various)',
    'Skp2'='Cells during development (various)',
    'pirk'='Immune response or development',
    'Rel'='Immune response or development',
    'E(spl)mbeta-HLH'='Cells during development (various)',
    'E(spl)malpha-BFM'='Cells during development (various)',
    'E(spl)m3-HLH'='Cells during development (various)',
    'Myc'='Cells during development (various)',
    'mtgo'='Cells during development (various)',
    'Dtg'='Cells during development (various)',
    'klu'='Cells during development (various)',
    'path'='Cells during development (various)',
    'LanA'='Cells during development (various)',
    'ImpL2'='Cells during development (various)',
    'Jon65Aiii'='Tissue (enzyme)',
    'Sox21a'='Cells during development (various)',

    # Myoblast fusion and somatic muscle development
    'Ldh'='Cells during development (various)',
    'luna'='Cells during development (various)',

    'Jon25Biii'='Tissue (enzyme)'
  ),
  c(
    'Stem cell',
    'Cells during development (various)',
    'Immune response or development',
    'Tissue (endocrine)',
    'Tissue (enzyme)'
  )
)

plot_midgut_de_panel <- function(de_data, sgn=1, extent=3, overlap=0.25, limit=12, sigma=qnorm(0.975), sigma.plot=sigma, title=NULL) {
  de_data <- de_data %>% mutate(
    label = displayName %>% paste0(
      ' (',
      cut(de_data$q_val, c(0, 0.0001, 0.001, 0.01, 0.05))
      %>% fct_recode(
        '*' = '(0.01,0.05]',
        '**' = '(0.001,0.01]',
        '***' = '(0.0001,0.001]',
        '****' = '(0,0.0001]'
      ),
      ')'
    ),
    annotation = factor(
      de.labels[rowname],
      levels=c("Not characterized", levels(de.labels))
    ) %>%
      replace(is.na(.), "Not characterized")
  )
  # Retain lrRNA, srRNA (MT genes)
  de_data = subset(de_data, q_val < 0.05 & sign(log2FoldChange) == sgn & !grepl('[^ls]rRNA', label))
  de_data = de_data[head(order(with(de_data, log2FoldChange-sgn*sigma*lfcSE), decreasing=sgn>0), limit), ]
  if (nrow(de_data) > 0) de_data$number = seq(nrow(de_data))
  else de_data = cbind(de_data, data.frame(number = numeric(0)))
  (ggplot(de_data)
    + geom_segment(
      aes(
        y=number,
        yend=number,
        x=pmax(log2FoldChange-sigma.plot*lfcSE, ifelse(log2FoldChange>0, 0, -Inf)),
        xend=pmin(log2FoldChange+sigma.plot*lfcSE, ifelse(log2FoldChange<0, 0, Inf)),
        color=annotation
      ),
      linewidth=0.35
    )
    + geom_text(
      aes(
        x=sign(log2FoldChange) * pmax(0.25, abs(log2FoldChange)),
        y=number-0.5,
        label=label
      ),
      family='Noto Serif',
      size=2.5
    )
    + geom_point(aes(x=log2FoldChange, y=number, color=annotation), size=0.5)
    # Include all factor levels in every plot's legend.
    + geom_blank(aes(x=0, y=1, color=annotation), data.frame(annotation=de.labels))
    + scale_x_continuous(
      expand=rep(0.01,2)
    )
    + scale_y_reverse(breaks=NULL)
    + scale_color_manual(
      breaks = levels(de_data$annotation),
      values = c(
        '#000000',
        hsv(0.135, 0.95, 0.82),
        hsv(0.59, 0.95, 0.65),
        hsv(0.52, 0.55, 0.75),
        hsv(0.9, 0.6, 0.8),
        hsv(0.95, 1, 0.6)
      ),
      guide = guide_legend(title = 'Expression in:')
    )
    + labs(
      x = "log2(EB/ISC) fold change ± σ",
      y = "Ranked Genes"
    )
    + ggtitle(title)
    )
}
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

de.labels = factor(
  c(
    'Amy-p'='Absorptive EC enzyme',
    'MtnC'='Absorptive EC enzyme',
    'Dl'='Stem cell',
    'His3.3A'='Stem cell',
    # Development 2020
    'baf'='Stem cell',
    'trol'='Development',
    'Debcl'='Development',
    'Skp2'='Development',
    'Rel'='Development',
    'E(spl)mbeta-HLH'='Development',
    'E(spl)malpha-BFM'='Development',
    'E(spl)m3-HLH'='Development',
    'Myc'='Development',
    'mtgo'='Development',
    'Dtg'='Development',
    'klu'='Development',
    'path'='Development',
    'LanA'='Development',
    'ImpL2'='Development',
    'Jon65Aiii'='Absorptive EC enzyme',
    'Sox21a'='Development',

    # Myoblast fusion and somatic muscle development
    'Ldh'='Development',
    'luna'='Development',

    'Jon25Biii'='Absorptive EC enzyme',
    'Jon65Aiv'='Absorptive EC enzyme'
  ),
  c(
    'Stem cell',
    'Development',
    'Absorptive EC enzyme'
  )
)

plot_arrange_deg_model_color_panels <- function(
  deg_model,
  sigma=qnorm(0.975),
  sigma.plot=sigma,
  limits=c(-4.5, 4.5)
) {
  de_data <- deg_model %>%
    display_deg_data %>%
    mutate(
      label = displayName %>% paste0(
        ' (',
        cut(q_val, c(0, 0.0001, 0.001, 0.01, 0.05))
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
    ) %>%
    # Retain lrRNA, srRNA (MT genes)
    filter(q_val < 0.05, !grepl('[^ls]rRNA', label))
  de_data <- de_data %>%
    arrange(desc(pmax(0, abs(log2FoldChange) - sigma*lfcSE)), desc(abs(log2FoldChange))) %>%
    group_by(sign(log2FoldChange)) %>%
    dplyr::slice(1:12) %>%
    mutate(y = factor(seq_along(label)))

  (
    ggplot(de_data)
    + geom_rect(
      aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, fill=celltype),
      tribble(
        ~xmin, ~ymin, ~xmax, ~ymax, ~celltype,
        -Inf, -Inf, 0, Inf, "ISC",
        0, -Inf, Inf, Inf, "EB"
      )
    )
    + geom_segment(
      aes(
        y,
        yend=y,
        x=pmax(log2FoldChange-sigma.plot*lfcSE, ifelse(log2FoldChange>0, 0, -Inf)),
        xend=pmin(log2FoldChange+sigma.plot*lfcSE, ifelse(log2FoldChange<0, 0, Inf)),
        color=annotation
      ),
      linewidth=0.25
    )
    + geom_point(
      aes(log2FoldChange, y, color=annotation),
      size=0.25
    )
    + geom_text(
      aes(
        x=sign(log2FoldChange) * pmax(0.25, abs(log2FoldChange)),
        y=max(as.numeric(y)) + 1 - as.numeric(y) + 0.5,
        label=label
      ),
      family='Noto Serif',
      size=2
    )
    + scale_x_continuous(
      name=bquote(log[2]*"(EB/ISC) posterior 95% CI"),
      breaks=pretty_breaks(5)
    )
    + scale_y_discrete(
      name="DE Genes",
      limits=rev,
      breaks=NULL,
      expand=c(0,0)
    )
    + coord_cartesian(
      limits,
      c(0.6, 0.9 + length(levels(de_data$y))),
      # Clip off, for annotate custom tag on the plot later.
      clip = "off"
    )
    + scale_color_manual(
      breaks = levels(de_data$annotation),
      values = c(
        '#000000',
        # Stem cell
        "#DEB417", # hcl(65, 85, 75)
        # Development
        "#F17D50", # hcl(25, 100, 65),
        # Absorptive EC enzyme
        "#3B9F9B" # hcl(187, 40, 60)
      ),
      na.value = "#000000",
      limits = factor(c("Stem cell", "Development", "Absorptive EC enzyme")),
      guide = guide_legend(title = 'Characterized as a marker:')
    )
    + scale_fill_manual(
      guide=NULL,
      values=muted(midgut.colors[1:2], 95, 20)
    )
    + theme(
      axis.ticks=element_line(linewidth=0.1),
      legend.position = "bottom",
      panel.ontop = TRUE,
      panel.background = element_rect(fill="transparent")
    )
  )
}
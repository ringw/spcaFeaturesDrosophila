# NOTCH2 / NOTCH3
# MYB / MYC
# with(as.data.frame(acc[['spca']]@cell.embeddings), cor(SPARSE_1, SPARSE_14))
# with(as.data.frame(acc[['spca']]@cell.embeddings), cor(SPARSE_3, SPARSE_26))
# with(as.data.frame(acc[['spca']]@cell.embeddings), cor(SPARSE_12, SPARSE_11))

acc_annotated_figure <- function(acc, variable = "proportion_cnv", guide = "P(CNV)") {
  # Add some useful columns here
  acc$proportion_cnv = 1 - rowProds(1 - (acc %>% FetchData(paste0('proportion_cnv_chr', 1:22)) %>% as.matrix))
  for (n in levels(acc$individual))
    acc@meta.data <- acc@meta.data %>%
      cbind(matrix(as.numeric(acc$individual == n), ncol=1, dimnames=list(NULL, n)))
  acc[['umap.spca']]@cell.embeddings %>%
    cbind(acc@meta.data) %>%
    ggplot(aes(UMAP_1, UMAP_2, color=get(variable))) + rasterize(
      geom_point(size = 0.08), dpi = 600
    ) + scale_color_viridis_c(
      option = 'rocket', begin = 0.05, end = 0.8,
      limits = c(0, 1), labels = scales::percent,
      guide = guide_colorbar(title = guide)
    ) + annotate(
      "rect", xmin = -2.75, xmax = 5, ymin = -0.75, ymax = 8.5, color = 'black', fill = 'transparent'
    ) + annotate(
      "text", 1.15, 9.5, label = "ACC"
    ) + annotate(
      "rect", xmin = -2.75, xmax = 6.5, ymin = -9.75, ymax = -4, color = 'black', fill = 'transparent'
    ) + annotate(
      "text", 1.75, -3, label = "CAF"
    ) + coord_cartesian(
      c(NA, NA), c(-11, 10)
    ) + theme_bw()
}

acc_fetch_feature <- function(acc, feature) {
  data <- cbind(as.data.frame(acc[['umap.spca']]@cell.embeddings), LogNormalize =  acc %>% FetchData(feature) %>% pull(feature))
  bind_rows(
    list(
      # Rotate 90 degrees clockwise so that DLL1+ cells are on the left.
      ACC = data %>%
        subset(between(UMAP_1, -2.75, 5) & between(UMAP_2, -0.75, 8.5)) %>%
        as.matrix %>%
        `%*%`(
          matrix(c(0,1,0, -1,0,0, 0,0,1), ncol=3, dimnames=list(NULL, colnames(.)))
        ) %>%
        as.data.frame,
      CAF = data %>%
        subset(between(UMAP_1, -2.75, 6.5) & between(UMAP_2, -9.75, -4))
    ),
    .id = "inset"
  )
}

acc_fetch_f2 <- function(acc, feature1, feature2) {
  sapply(
    c(feature1, feature2),
    \(n) acc_fetch_feature(acc, n),
    simplify = F
  ) %>%
    bind_rows(.id = "feature") %>%
    ggplot(aes(UMAP_1, UMAP_2)) + facet_grid(vars(inset), vars(1), scales='free')
}

acc_blend_features <- function(acc, feature1, feature2) {
  n_colors <- 256
  scale1 = viridis(n_colors, end=0.9, option="viridis")
  scale2 = viridis(n_colors, begin=0.1, option="magma")
  blend <- mixcolor(
    0.6,
    hex2RGB(scale1) %>% as("LAB"), hex2RGB(scale2) %>% as("LAB")
  ) %>% hex

  features <- as.data.frame(acc[['umap.spca']]@cell.embeddings) %>%
    cbind(
      acc %>% FetchData(c(feature1, feature2))
    )

  feature_maxs <- colQuantiles(features %>% as.matrix, probs=0.95)
  feature_mapping <- as.matrix(features) %>%
    `%*%`(matrix(diag(1 / feature_maxs), ncol=ncol(.), dimnames=rep(list(colnames(.)), 2)))
  feature_mapping <- round(
    feature_mapping * (n_colors - 1) + 1
  ) %>% pmin(n_colors)
  colors <- mixcolor(
    0.6,
    hex2RGB(scale1[feature_mapping[, feature1]]) %>% as("LAB"),
    hex2RGB(scale2[feature_mapping[, feature2]]) %>% as("LAB")
  ) %>% hex
  features <- features %>% cbind(color = colors)

  ggplot(features, aes(UMAP_1, UMAP_2)) + geom_point(color=features$color, size=2) + geom_tile(aes(color=C1), data.frame(UMAP_1=0, UMAP_2=0, C1=0), alpha=0, width=0, height=0) + scale_color_viridis_c(limits = c(0, feature_maxs[feature1]), end=0.9) + new_scale_color() + geom_tile(aes(color=C1), data.frame(UMAP_1=0, UMAP_2=0, C1=0), alpha=0, width=0, height=0) + scale_color_gradientn(limits = c(0, min(feature_maxs[c(feature1,feature2)])), colors = blend) + new_scale_color() + geom_tile(aes(color=C1), data.frame(UMAP_1=0, UMAP_2=0, C1=0), alpha=0, width=0, height=0) + scale_color_viridis_c(limits = c(0, feature_maxs[feature2]), option = "magma", begin=0.1)
}

acc_feature1_colors <- function(n=256) viridis(n, begin=0.1, end=0.9, option="magma")
acc_feature2_colors <- function(n=256) viridis(n, end=0.75, option="viridis")

acc_blend_features_inset <- function(acc, feature1, feature2, cells = NULL, show.legend = TRUE) {
  data <- acc_fetch_feature(acc, feature1)
  data$feature2 <- acc_fetch_feature(acc, feature2)[, 'LogNormalize']

  # Insetting already causes us to subset the cells in acc_fetch_feature
  cells = cells %>% intersect(rownames(data))

  feature_maxs <- colQuantiles(data %>% subset(select=-c(inset,UMAP_1,UMAP_2)) %>% as.matrix, probs=0.95)

  bounding_box <- data %>%
    subset(select=c(inset,UMAP_1,UMAP_2)) %>%
    group_by(inset) %>%
    summarise_all(list(center = ~ mean(range(.)), width = ~ diff(range(.))))
  colnames(bounding_box) <- colnames(bounding_box) %>% str_replace("_center", "")

  if (!is.null(cells)) data <- data[cells, ]

  n_colors <- 256
  blend <- mixcolor(
    0.6,
    hex2RGB(acc_feature1_colors()) %>% as("LAB"),
    hex2RGB(acc_feature2_colors()) %>% as("LAB")
  ) %>% hex

  # Features should have a min of 0; we are taking a large quantile here and
  # dividing.
  data[, c('LogNormalize','feature2')] <- data[, c('LogNormalize','feature2')] %>%
    `-`(
      rep(colMins(as.matrix(.)), each=nrow(.))
    )
  feature_mapping <- data %>% subset(select=-c(inset,UMAP_1,UMAP_2)) %>% as.matrix %>%
    `%*%`(matrix(diag(1 / feature_maxs), ncol=ncol(.)))
  feature_mapping <- round(
    feature_mapping * (n_colors - 1) + 1
  ) %>% pmin(n_colors)
  colors <- mixcolor(
    feature_mapping[, 2] / rowSums(feature_mapping),
    hex2RGB(acc_feature1_colors()[feature_mapping[, 1]]) %>% as("LAB"),
    hex2RGB(acc_feature2_colors()[feature_mapping[, 2]]) %>% as("LAB")
  ) %>% hex
  data <- data %>% cbind(color = colors)

  g <- ggplot(data, aes(UMAP_1, UMAP_2)) + facet_wrap(
    vars(inset), ncol = 1, scales='free'
  ) + rasterize(
    geom_point(
      color=data$color, size=0.16
    ),
    dpi = 600
  ) + theme_bw() + scale_x_continuous(
    labels = NULL
  ) + scale_y_continuous(
    labels = NULL
  ) + labs(
    x = "UMAP (rotate)", y = "UMAP (rotate)"
  )
  if (show.legend)
    g <- g + geom_tile(
      aes(color=c(0, feature_maxs[1]), width=UMAP_1_width, height=UMAP_2_width),
      bounding_box %>% cbind(C1 = 0),
      fill = "transparent", linewidth = 0
    ) + scale_color_viridis_c(
      breaks = scales::pretty_breaks(2),
      limits = c(0, feature_maxs[feature2]), option = "magma", begin=0.1, end=0.9,
      guide = guide_colorbar(title = display_gene_names(feature1), barheight = 3)
    ) + new_scale_color() + geom_tile(
      aes(color=C1), data.frame(UMAP_1=-Inf, UMAP_2=-Inf, C1=c(0, min(feature_maxs))), alpha=0,
      width=0, height=0
    ) + scale_color_gradientn(
      breaks = scales::pretty_breaks(2),
      limits = c(0, min(feature_maxs[c(feature1,feature2)])),
      colors = blend, guide = guide_colorbar(title = "both", barheight = 3)
    ) + new_scale_color() + geom_tile(
      aes(color=C1), data.frame(UMAP_1=-Inf, UMAP_2=-Inf, C1=c(0, feature_maxs[2])), alpha=0,
      width=0, height=0
    ) + scale_color_viridis_c(
      breaks = scales::pretty_breaks(2),
      limits = c(0, feature_maxs[feature1]), end=0.75,
      guide = guide_colorbar(title = display_gene_names(feature2), barheight = 3)
    )
  g
}

acc_arrange_figure <- function(acc, individual = NULL) {
  if (is.null(individual))
    cells <- NULL
  else
    cells <- Cells(acc) %>% subset(as.character(acc$individual) == individual)
  ggarrange(
    if (is.null(individual)) acc_annotated_figure(acc)
    else acc_annotated_figure(acc, individual, guide = individual),
    acc_blend_features_inset(acc, "NOTCH3", "DLL1", cells = cells),
    acc_blend_features_inset(acc, "N3+", "N1N2", cells = cells),
    acc_blend_features_inset(acc, "SPARSE_14", "SPARSE_1", cells = cells),
    acc_blend_features_inset(acc, "SPARSE_26", "SPARSE_3", cells = cells),
    acc_blend_features_inset(acc, "SPARSE_11", "SPARSE_12", cells = cells),
    nrow = 1,
    widths = c(2,1,1,1,1,1)
  )
}
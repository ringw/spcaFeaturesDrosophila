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

# Add the idents from the annotated figure.
acc_add_components_from_umap <- function(acc) {
  Idents(acc) <- with(
    acc[['umap.spca']]@cell.embeddings %>% as.data.frame,
    between(UMAP_1, -2.75, 5) & between(UMAP_2, -0.75, 8.5)
  ) %>%
    ifelse(
      "ACC",
      with(
        acc[['umap.spca']]@cell.embeddings %>% as.data.frame,
        between(UMAP_1, -2.75, 6.5) & between(UMAP_2, -9.75, -4)
      ) %>%
        ifelse(
          "CAF",
          "others"
        )
    )
  acc
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

scaleNonneg <- function(mat) {
  mat <- as.matrix(mat)
  (
    DelayedArray(mat) * RleArray(
      Rle(
        1 / colSds(mat),
        rep(nrow(mat), ncol(mat))
      ),
      dim(mat)
    )
  ) %>% as.matrix %>% round(digits = 5)
}

acc_make_feature_profiles <- function(acc, features) {
  features <- acc %>% FetchData(features) %>% as.matrix %>% scaleNonneg %>% as.data.frame
  cdf_rowData <- expand.grid(
    x = seq(-1, 8, by = 0.05),
    ident = c("ACC", "CAF")
  )[, c(2, 1)]
  cdf_matrix <- apply(
    features,
    2,
    \(y) cdf_rowData %>%
      left_join(
        sapply(
          c("ACC", "CAF"),
          \(n) data.frame(
            ident = n,
            x = seq(-1, 8, by = 0.05),
            y = y %>%
              subset(as.character(Idents(acc)) == n) %>%
              cut(c(-Inf, seq(-1, 8, by = 0.05))) %>%
              table %>%
              cumsum %>% `/`(
                sum(as.character(Idents(acc)) == n)
              )
          ),
          simplify = F
        ) %>%
          bind_rows(.id = "ident"),
        by = join_by(ident, x)
      ) %>%
      pull(y)
  )
  rownames(cdf_matrix) <- interaction(cdf_rowData$ident, cdf_rowData$x)
  cl = kmeans(t(cdf_matrix), 3, ns=100L)
  cdf_colData <- data.frame(profile = factor(cl$cluster, sort(unique(cl$cluster))))
  ident_profile <- sapply(
    c("ACC", "CAF"),
    \(n) sapply(
      levels(cdf_colData$profile),
      \(l) mean(
        cdf_matrix[
          which.max(
            cdf_rowData$ident == n
            & cdf_rowData$x == 4
          ),
          as.character(cdf_colData$profile) == l
        ]
      )
    ) %>% sort %>% head(1) %>% names
  )
  names(ident_profile) <- names(ident_profile) %>% paste0("_high")
  ident_profile <- c(
    ident_profile,
    others = levels(cdf_colData$profile)[
      which.max(!(levels(cdf_colData$profile) %in% ident_profile))
    ]
  )
  cdf_colData$profile <- cdf_colData$profile %>% list %>% append(
    as.list(ident_profile)
  ) %>% do.call(fct_recode, .) %>% fct_relevel(
    "ACC_high", "CAF_high"
  )
  SingleCellExperiment(
    cdf_matrix,
    colData = cdf_colData,
    rowData = cdf_rowData
  )
}

acc_make_feature_profiles_density <- function(acc, features) {
  features <- acc %>% FetchData(features) %>% as.matrix %>% scaleNonneg %>% as.data.frame
  cdf_rowData <- expand.grid(
    x = seq(-1, 8, by = 0.05),
    ident = c("ACC", "CAF")
  )[, c(2, 1)]
  cdf_matrix <- apply(
    features,
    2,
    \(y) cdf_rowData %>%
      left_join(
        sapply(
          c("ACC", "CAF"),
          \(n) data.frame(
            ident = n,
            bw = 0.2,
            x = seq(-1, 8, by = 0.05),
            y = (
              y %>%
                subset(as.character(Idents(acc)) == n) %>%
                density(from = -1, to = 8, n = 9 * 20 + 1)
            )$y
          ),
          simplify = F
        ) %>%
          bind_rows(.id = "ident"),
        by = join_by(ident, x)
      ) %>%
      pull(y)
  )
  rownames(cdf_matrix) <- interaction(cdf_rowData$ident, cdf_rowData$x)
  cl = kmeans(t(cdf_matrix), 3, ns=100L)
  cdf_colData <- data.frame(profile = factor(cl$cluster, sort(unique(cl$cluster))))
  ident_profile <- sapply(
    c("ACC", "CAF"),
    \(n) sapply(
      levels(cdf_colData$profile),
      \(l) mean(
        cdf_matrix[
          which.max(
            cdf_rowData$ident == n
            & cdf_rowData$x == 4
          ),
          as.character(cdf_colData$profile) == l
        ]
      )
    ) %>% sort %>% head(1) %>% names
  )
  names(ident_profile) <- names(ident_profile) %>% paste0("_high")
  ident_profile <- c(
    ident_profile,
    others = levels(cdf_colData$profile)[
      which.max(!(levels(cdf_colData$profile) %in% ident_profile))
    ]
  )
  cdf_colData$profile <- cdf_colData$profile %>% list %>% append(
    as.list(ident_profile)
  ) %>% do.call(fct_recode, .) %>% fct_relevel(
    "ACC_high", "CAF_high"
  )
  SingleCellExperiment(
    cdf_matrix,
    colData = cdf_colData,
    rowData = cdf_rowData
  )
}

acc_plot_feature_profiles <- function(acc.spca.profiles) {
  acc.spca.profiles %>% assay %>%
    melt(c("label", "feature"), value.name = "y") %>%
    left_join(
      acc.spca.profiles %>% colData %>%
        as.data.frame %>%
        rownames_to_column("feature"),
      by = "feature"
    ) %>%
    mutate(
      ident = str_extract(label, "(ACC|CAF)"),
      x = str_extract(label, "\\.(.+)", group = 1) %>% as.numeric
    ) %>%
    ggplot(aes(x, y, color=ident)) + facet_wrap(vars(profile)) + geom_line(aes(group=interaction(feature, ident)), alpha=0.25) + geom_smooth(linewidth = 2, method="loess")
}

acc_pca_profiles <- function(acc) {
  apply(
    acc[['pca']]@cell.embeddings[, 1:10],
    2,
    \(v) tribble(
      ~ident, ~sign,
      "ACC", 1,
      "ACC", -1,
      "CAF", 1,
      "CAF", -1
    )[
      which.max(
        as.numeric(
          cor(
            cbind(
              positive = v %>% pmax(0),
              negative = v %>% pmin(0)
            ),
            cbind(
              ACC = Idents(acc) == "ACC",
              CAF = Idents(acc) == "CAF"
            )
          )
        )
      ),
    ] %>%
      as.list
  )
  with(
    acc[['pca']]@cell.embeddings %>% as.data.frame,
    list(
      CAF=cbind(
        PC1=PC_1,
        PC2=-PC_2,
        PC8=PC_8
      ) %>% scale,
      ACC=cbind(
        PC3=PC_3,
        PC5=-PC_5,
        PC6=-PC_6,
        PC9=-PC_9
      ) %>% scale
    )
  )
}

acc_spca_profiles <- function(acc) {
  # tar_load(acc.spca.profiles)
  # acc.spca.profiles$profile
  caf_spcs <- c(1, 3, 7, 12, 13, 17, 23, 25, 29, 36, 38, 45)
  acc_spcs <- c(11, 14, 26, 37, 41, 49)
  list(
    CAF=acc[['spca']]@cell.embeddings[, caf_spcs] %>% scale,
    ACC=acc[['spca']]@cell.embeddings[, acc_spcs] %>% scale
  )
}

acc_plot_feature_profiles <- function(acc, m) {
  names(dimnames(m)) <- c("cell", "feature")
  x <- seq(-5, 5, by=0.05)
  profiles <- sapply(
    c("ACC", "CAF"),
    \(n) apply(
      m[as.character(Idents(acc)) == n, ],
      2,
      \(y) data.frame(
        x = x,
        y = ecdf(y)(x)
      ),
      simplify = FALSE
    ) %>% bind_rows(.id = "feature"),
    simplify = FALSE
  ) %>% bind_rows(.id = "ident")
  profiles
  ggplot(
    profiles,
    aes(x, y, group=interaction(feature, ident), color=ident)
  ) + geom_line()
}

acc_nonneg_feature_plot_gradient <- c(hcl(0, 0, 87), hcl(85, 33, 80), hcl(68, 37, 75), hcl(52, 40, 70), hcl(50, 57, 70), hcl(47, 64, 68), hcl(30, 84, 56), hcl(12, 103, 45))
nonneg_feature_plot_annotate <- function(acc, feature, max_scale, annotations=NULL, subset=NULL) {
  g <- cbind(
    data.frame(feature = FetchData(acc, feature) %>% pull(feature)),
    acc[['umap.spca']]@cell.embeddings
  ) %>% ggplot(aes(x0 = UMAP_1, y0 = UMAP_2, color=feature)) + rasterize(
    geom_circle(
      aes(r = 0.005),
      data = subset
    ),
    dpi = 300
  ) + scale_color_gradientn(
    colors = acc_nonneg_feature_plot_gradient,
    limits=c(0, if (is.null(max_scale)) NA else max_scale), oob=scales::squish,
    guide = guide_colorbar(title = feature %>% display_gene_names)
  ) + coord_cartesian(
    c(-3, 6.5), c(-10, 9), expand=F
  ) + scale_y_continuous(
    breaks = c(-5, 0, 5)
  ) + theme_bw() + theme(
    axis.title.y = element_text(margin = margin()),
    axis.text.y = element_text(margin = margin(t = 5.5, r = 2, b = 5.5))
  )
  for (an in annotations) {
    # g <- g + append(
      # list(color="yellow", fill="transparent", linewidth=2),
      # as.list(an)
    # ) %>% do.call(annotate, .)
    # g <- g + annotate(color="yellow", fill="transparent", linewidth=2, "tile", 0, 0, width=3, height=3)
    g <- g + annotate(
      color="yellow",
      fill="transparent",
      linewidth=2,
      geom=an$geom,
      x=an$x,
      y=an$y,
      width=an$width,
      height=an$height
    )
  }
  g
}

nonneg_feature_plot_query <- function(
  acc, feature, max_scale, ident
) {
  nonneg_feature_plot_annotate(
    acc, feature, max_scale,
    subset=\(df) df[Idents(acc) == ident, ]
  ) + new_scale_color() + rasterize(
    geom_circle(
      aes(r = 0.005, color = feature),
      cbind(
        data.frame(feature = rep(0, ncol(acc))),
        acc[['umap.spca']]@cell.embeddings
      ) %>%
        subset(Idents(acc) != ident)
    ),
    dpi = 300
  ) + scale_color_gradient2(
    mid = hcl(97, 12, 93)
  )
}

acc_dim_plot_individual <- function(acc, feature) {
  g <- cbind(
    data.frame(feature = FetchData(acc, feature) %>% pull(feature)),
    acc[['umap.spca']]@cell.embeddings
  ) %>% ggplot(aes(x = UMAP_1, y = UMAP_2, color=feature)) + rasterize(
    geom_point(
      size = 0.1,
      data = subset
    ),
    dpi = 300
  ) + scale_color_viridis_d(
    option = "turbo", begin = 0.1, end = 0.9,
    guide = guide_legend(
      override.aes = list(size = 3)
    )
  ) + coord_cartesian(
    c(-3, 6.5), c(-10, 9), expand=F
  ) + scale_y_continuous(
    breaks = c(-5, 0, 5)
  ) + theme_bw() + theme(
    axis.title.y = element_text(margin = margin()),
    axis.text.y = element_text(margin = margin(t = 5.5, r = 2, b = 5.5))
  )
  g
}
pickClusteringResolution <- function(indrop, snn_obj, res_range = 100:200, by = 0.01) {
  eval_fct_isc_eb <- function(fct) {
    indrop.copy <- indrop
    indrop.copy$feature <- fct
    aves <- AverageExpression(
      indrop.copy,
      features=c(
        # nuisance features (found in spca clustering) to exclude clusters:
        "IA-2", "MtnA",
        # nuisance features found in pca clustering
        "Amy-p",
        # features of interest: Dl, esg (stem), pros (EE)
        "Dl", "esg", "pros"
      ),
      gr="feature")$RNA
    aves <- aves[, order(aves["esg", ], decreasing=TRUE)]
    aves[, 1:5] %>% round(3)
    ifelse(
      sum(
        # Nuisance genes...
        # aves['IA-2',] < 2 & aves['MtnA',] < 5 & aves['Amy-p',] < 5 &
        aves['pros',] < 1
        # Marker genes
        & aves['esg',] > 5 & aves['Dl',] > 0.5
      ) >= 2,
      2 * var(aves["Dl", aves["pros",] < 1 & aves["esg",] > 2]),
      -1
    )
  }
  # This is not actually a log2 value, and the pseudo-log fold change can
  # actually be quite large. Set it to an optimistic value and we will try to
  # attain this log-fold change value for Delta.
  target_delta_fold_change <- 2.5
  gtools::binsearch(
    \(res_int) sapply(
      0:9,
      \(random.seed) {
        delta_fold_change <- FindClusters(
          snn_obj,
          res = res_int * by,
          random.seed = random.seed,
          verb = F
        ) %>%
          pull(1) %>%
          eval_fct_isc_eb
        if (delta_fold_change < 0)
          delta_fold_change
        else
          (delta_fold_change - target_delta_fold_change)^2
      }
    ) %>% mean,
    range(res_range)
  )
}

eval_isc_clusters_in_fct <- function(indrop, fct) {
  indrop.copy <- indrop
  indrop.copy$feature <- fct
  aves <- AverageExpression(
    indrop.copy,
    features=c(
      # nuisance features (found in spca clustering) to exclude clusters:
      "IA-2", "MtnA",
      # nuisance features found in pca clustering
      "Amy-p",
      # features of interest: Dl, esg (stem), pros (EE)
      "Dl", "esg", "pros"
    ),
    gr="feature")$RNA
  aves <- aves[
    ,
    # Nuisance EE gene
    aves['pros', ] < 1
    # esg (high in ISC and EB) and Dl (ISC high, EB low, others off)
    & aves['esg', ] > 5 & aves['Dl', ] > 0.5,
    drop=FALSE
  ]
  colnames(aves)[order(aves['Dl', ], decreasing=TRUE)]
}

find_isc_clustering_res <- function(
  indrop, snn_obj, random.seed = 0:9, res_range = 20:200, by = 0.01, intended_count = 2
) {
  gtools::binsearch(
    \(res_int) sapply(
      random.seed,
      \(random.seed) {
        FindClusters(
          snn_obj,
          res = res_int * by,
          random.seed = random.seed,
          verb = F
        ) %>%
          pull(1) %>%
          eval_isc_clusters_in_fct(indrop, .) %>%
          length
      }
    ) %>% mean - intended_count + 0.2,
    range(res_range)
  )
}

countNumIscClusters <- function(indrop, snn_obj, res_range = 100:200, by = 0.01) {
  eval_fct_isc_eb <- function(fct) {
    indrop.copy <- indrop
    indrop.copy$feature <- fct
    aves <- AverageExpression(
      indrop.copy,
      features=c(
        # nuisance features (found in spca clustering) to exclude clusters:
        "IA-2", "MtnA",
        # nuisance features found in pca clustering
        "Amy-p",
        # features of interest: Dl, esg (stem), pros (EE)
        "Dl", "esg", "pros"
      ),
      gr="feature")$RNA
    aves <- aves[, order(aves["esg", ], decreasing=TRUE)]
    aves[, 1:5] %>% round(3)
    ifelse(
      sum(
        # Nuisance genes...
        # aves['IA-2',] < 2 & aves['MtnA',] < 5 & aves['Amy-p',] < 5 &
        aves['pros',] < 1
        # Marker genes
        & aves['esg',] > 5 & aves['Dl',] > 0.5
      ) >= 2,
      2 * var(aves["Dl", aves["pros",] < 1 & aves["esg",] > 2]),
      -1
    )
  }
  # This is not actually a log2 value, and the pseudo-log fold change can
  # actually be quite large. Set it to an optimistic value and we will try to
  # attain this log-fold change value for Delta.
  target_delta_fold_change <- 2.5
  gtools::binsearch(
    \(res_int) sapply(
      0:9,
      \(random.seed) {
        delta_fold_change <- FindClusters(
          snn_obj,
          res = res_int * by,
          random.seed = random.seed,
          verb = F
        ) %>%
          pull(1) %>%
          eval_fct_isc_eb
        if (delta_fold_change < 0)
          delta_fold_change
        else
          (delta_fold_change - target_delta_fold_change)^2
      }
    ) %>% mean,
    range(res_range)
  )
}

binsearchChooseX <- function(binsearch) {
  if (all(binsearch$value < 0)) return(NA)
  where <- binsearch$where[
    order(
      # put FALSEs (positive value) before TRUEs (-1)
      binsearch$value < 0,
      # put small value (least squares) first
      abs(binsearch$value)
    )
  ]
  where[1]
}

clusterAndLabelIscEb <- function(indrop, snn_obj, res, random.seed) {
  fct <- snn_obj %>%
    FindClusters(res = res, random.seed = random.seed, verb = F) %>%
    pull(1)
  labelIscEbFct(indrop, fct)
}

labelIscEbFct <- function(indrop, fct) {
  indrop.copy <- indrop
  indrop.copy$feature <- fct
  aves <- AverageExpression(
    indrop.copy,
    features=c(
      "Dl", "esg", "pros"
    ),
    gr="feature")$RNA
  levels.of.interest <- (
    aves["pros",] < 1 & aves["esg",] > 2
  ) %>%
    which %>%
    names %>%
    head(2)
  if (length(levels.of.interest) != 2)
    return(rep(NA, ncol(indrop)))
  fct <- fct %>%
    list %>%
    append(
      setNames(
        if (diff(aves["Dl", levels.of.interest]) > 0)
          rev(levels.of.interest)
        else
          levels.of.interest,
        c("ISC", "EB")
      )
    ) %>%
    do.call(fct_recode, .) %>%
    fct_relevel("ISC", "EB")
  fct
}

guessRandIndex <- function(fct1, fct2, mc = 50) {
  levels1 <- setdiff(levels(fct1), c("ISC", "EB"))
  levels2 <- setdiff(levels(fct2), c("ISC", "EB"))
  if (length(levels2) > length(levels1)) {
    items <- list(fct1, fct2)
    lvls <- list(levels1, levels2)
    fct1 <- items[[2]]
    fct2 <- items[[1]]
    levels1 <- lvls[[2]]
    levels2 <- lvls[[1]]
  }

  # levels1 can be matched to the fct2 "ISC" or "EB" cluster if necessary, but
  # the fct2 "ISC" and "EB" cells will each be matched to the fct1 ISC and EB.
  tbl <- table(fct1, fct2)[levels1, ] %>%
    `/`(rowSums(.))

  sapply(
    seq(mc),
    \(seed.use) seed.use %>%
      with_seed(
        {
          assignment1 <- apply(
            tbl,
            1,
            \(v) levels(fct2)[sample(length(v), 1, prob = v)]
          )
          reassign1 <- append(list(fct1), assignment1) %>% do.call(dplyr::recode, .)
          adj.rand.index(reassign1, fct2)
        }
      )
  ) %>%
    max
}

corIscEb <- function(fct1, fct2) {
  makeIscEbIndicator <- \(v) ifelse(
    v == "EB",
    1,
    ifelse(
      v == "ISC",
      -1,
      0
    )
  )
  cor(makeIscEbIndicator(fct1), makeIscEbIndicator(fct2))
}
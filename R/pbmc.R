pbmc_to_annotated_signac <- function(counts_h5, metadata, fragments_tsv, hsapiens_annotations) {
  counts <- Read10X_h5(filename = counts_h5)
  metadata <- metadata %>%
    read.csv(
      header = TRUE,
      row.names = 1
    )

  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = fragments_tsv,
    min.cells = 10,
    min.features = 200
  )

  pbmc <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )

  Annotation(pbmc) <- hsapiens_annotations
  pbmc
}

annotated_signac_to_activity_seurat <- function(pbmc, annotations) {
  activity <- pbmc %>%
    GeneActivity %>%
    CreateSeuratObject(counts = ., assay = "ACTIVITY") %>%
    NormalizeData %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData
  activity@misc$covar <- activity[["ACTIVITY"]]@scale.data[
    VariableFeatures(activity),
  ] %>%
    t %>%
    var
  activity@misc$covar <- activity@misc$covar * (
    create_genes_non_overlapping_mask(
      rownames(activity@misc$covar), annotations
    )
  )
  activity
}

create_genes_non_overlapping_mask <- function(gene_list, annotations) {
  annotations <- annotations %>%
    subset(gene_name %in% gene_list)
  chrs <- annotations %>% split(seqnames(annotations))
  chrAntiMasks <- sapply(
    chrs,
    \(granges) {
      antiMaskAnnList <- findOverlaps(granges, granges) %>% as("List")
      chrGeneList <- unique(granges$gene_name)
      antiMaskIndList <- (
        antiMaskAnnList %>%
          split(granges$gene_name) %>%
          mapply(
            \(gRangeName, gRangeInds) granges$gene_name[
              gRangeInds %>%
                do.call(c, .) %>%
                unique
            ] %>%
              unique %>%
              setdiff(gRangeName) %>%
              match(chrGeneList) %>%
              unique,
            names(.),
            .,
            SIMPLIFY=F
          )
      )[chrGeneList]
      antiMask <- sparseMatrix(
        i = do.call(c, antiMaskIndList),
        j = rep(seq_along(antiMaskIndList), sapply(antiMaskIndList, length)),
        dims = rep(length(chrGeneList), 2),
        dimnames = rep(list(chrGeneList), 2)
      )
    },
    simplify=FALSE
  )
  antiMask <- .bdiag(chrAntiMasks)
  dimnames(antiMask) <- rep(
    do.call(c, sapply(chrAntiMasks, rownames, simplify=F)) %>% list,
    2
  )
  as.matrix(1 - antiMask[gene_list, gene_list])
}

pbmc_urls <- list(
  tar_download(
    pbmc.atac.filtered.peak.tar.gz,
    "https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.tar.gz",
    "10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.tar.gz",
    cue = tar_cue("never")
  ),
  tar_file(
    pbmc.atac.filtered.10x,
    tibble(
      output_path = "10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix",
      bak_path = paste0(output_path, "~"),
      apply_rename = if (file.exists(output_path)) file.rename(output_path, bak_path),
      dir_create = dir.create(output_path),
      extract_tar = run(
        "tar",
        c(
          "xf",
          pbmc.atac.filtered.peak.tar.gz,
          "-C",
          output_path,
          "--strip-components",
          "1"
        )
      ) %>%
        list
    )$output_path
  ),
  tar_file(
    grch38.p13.gtf,
    # accessed via:
    # https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/?HistoryId=MCID_6615bf3cd88ba65c89400b4c&QueryKey=5&ReleaseType=RefSeq&FileType=GENOME_GTF&Flat=true
    "GCF_000001405.39_GRCh38.p13_genomic.gtf.gz"
  ),
  tar_download(
    pbmc.atac.h5,
    "https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5",
    "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5",
    cue = tar_cue("never")
  ),
  tar_download(
    pbmc.atac.metadata,
    "https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv",
    "atac_v1_pbmc_10k_singlecell.csv",
    cue = tar_cue("never")
  ),
  tar_download(
    pbmc.atac.fragments.file,
    "https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz",
    "atac_v1_pbmc_10k_fragments.tsv.gz",
    cue = tar_cue("never")
  ),
  tar_download(
    pbmc.atac.fragments.index,
    "https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi",
    "atac_v1_pbmc_10k_fragments.tsv.gz.tbi",
    cue = tar_cue("never")
  ),
  tar_file(
    name = pbmc.atac.fragments,
    # the "index" is a dependency for the fragments file.
    c(pbmc.atac.fragments.file, pbmc.atac.fragments.index)[1]
  )
)

pbmc_targets <- list(
  tar_target(
    hsapiens_annotations,
    {
      annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
      seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
      genome(annotations) <- "hg19"
      annotations
    },
    packages = tar_option_get("packages") %>% c("EnsDb.Hsapiens.v75", "Signac")
  ),
  tar_target(
    pbmc.atac.activity,
    pbmc_to_annotated_signac(
      pbmc.atac.h5, pbmc.atac.metadata, pbmc.atac.fragments, hsapiens_annotations
    ) %>%
      annotated_signac_to_activity_seurat(
        hsapiens_annotations
      ),
    packages = tar_option_get("packages") %>% c("Signac")
  ),
  # Feature loadings - the linear transformation from 
  tar_target(
    pbmc.atac.spca.linear.transformation,
    seurat_spca_compute_feature_loadings(
      pbmc.atac.activity@misc$covar, varnum=5, npcs=80, eigen_gap=0.001,
      search_cap=500000, cgroup = memory_cgroups[1]
    )
  ),
  tar_target(
    pbmc.atac.spca.linear.transformation.k8,
    with_seed(
      -1588785896,
      seurat_spca_compute_feature_loadings(
        pbmc.atac.activity@misc$covar, varnum=8, npcs=50, eigen_gap=0.001,
        search_cap=500000, cgroup = memory_cgroups[2]
      )
    )
  ),
  tar_target(
    pbmc.atac.spca.linear.transformation.k4,
    with_seed(
      -1588785896,
      seurat_spca_compute_feature_loadings(
        pbmc.atac.activity@misc$covar, varnum=4, npcs=120, eigen_gap=0.001,
        search_cap=500000, cgroup = memory_cgroups[3]
      )
    )
  ),
  tar_target(
    pbmc.atac.spca,
    seurat_spca_from_feature_loadings(
      pbmc.atac.spca.linear.transformation,
      pbmc.atac.activity,
      "ACTIVITY",
      do.corr = T
    )
  )
)
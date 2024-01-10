library(ggpubr)
acc$is.pca.selected <- with(
  acc[['umap']]@cell.embeddings %>% as.data.frame,
  between(umap_1, -5,7) & between(umap_2, -2,8)
)
acc[['spca']][,'SPARSE_11'] %>% subset(. != 0) %>% abs %>% as.data.frame %>% arrange(desc(SPARSE_11)) %>% rownames
acc[['spca']][,'SPARSE_14'] %>% subset(. != 0) %>% abs %>% as.data.frame %>% arrange(desc(SPARSE_14)) %>% rownames
ggarrange(
  FeaturePlot(acc, 'Notch_Score') + annotate('rect', xmin=-5, ymin=-2, xmax=7, ymax=8, fill="#ffff0022")
  + labs(x = "PCA-UMAP, Parikh gene circuit: 'NRARP','NOTCH3','HES4','HEY1','HEY2'")
  + scale_color_viridis_c(option = "magma"),
  DimPlot(acc, red='umap.spca', cells.highlight = Cells(acc) %>% subset(acc$is.pca.selected)) + annotate("rect", xmin=-3, ymin=-1, xmax=5.5, ymax=9, fill="transparent", color="black"),
  FeaturePlot(acc, 'Notch_Score', red='umap.spca') + coord_cartesian(c(-3,5.5), c(-1,9))
  + labs(x = "PCA-UMAP, Parikh gene circuit: 'NRARP','NOTCH3','HES4','HEY1','HEY2'")
  + scale_color_viridis_c(option = "magma"),
  FeaturePlot(acc, 'SPARSE_11', red='umap.spca') + coord_cartesian(c(-3,5.5), c(-1,9))
  + labs(x = "SPC11: GABRP,AZGP1,CLDN4,CLDN3,FXYD3,CALML5,WFDC2,SLPI,KRT19,PRR15L")
  + scale_color_viridis_c(option = "magma"),
  FeaturePlot(acc, 'SPARSE_14', red='umap.spca') + coord_cartesian(c(-3,5.5), c(-1,9))
  + labs(x = "SPC14: KRT16,KRT14,KRT5,KRT16P6,TACSTD2,KRT17,KRT15,KRT6B,KRT7,MMP7")
  + scale_color_viridis_c(option = "magma"),
  ggplot(),
  ggplot(),
  FeaturePlot(acc, 'Notch_Score', red='umap.spca')
  + labs(x = "PCA-UMAP, Parikh gene circuit: 'NRARP','NOTCH3','HES4','HEY1','HEY2'")
  + scale_color_viridis_c(option = "magma"),
  FeaturePlot(acc, 'SPARSE_11', red='umap.spca')
  + labs(x = "SPC11: GABRP,AZGP1,CLDN4,CLDN3,FXYD3,CALML5,WFDC2,SLPI,KRT19,PRR15L")
  + scale_color_viridis_c(option = "magma"),
  FeaturePlot(acc, 'SPARSE_14', red='umap.spca')
  + labs(x = "SPC14: KRT16,KRT14,KRT5,KRT16P6,TACSTD2,KRT17,KRT15,KRT6B,KRT7,MMP7")
  + scale_color_viridis_c(option = "magma"),
  ncol = 5,
  nrow = 2
)
ggsave(
  'AccExample.pdf',
  width=40,
  height = 16
)

cbind(acc[['spca']]@cell.embeddings, Notch_Score=acc$Notch_Score) %>% scale %>% as.data.frame %>% lm(Notch_Score ~ SPARSE_11 * SPARSE_14, .) %>% summary

cbind(acc[['spca']]@cell.embeddings, Notch_Score=acc$Notch_Score) %>% scale %>% cbind(acc[['umap.spca']]@cell.embeddings) %>% as.data.frame %>% lm(Notch_Score ~ SPARSE_11 * SPARSE_14, ., subset = between(.$UMAP_1, -3,5.5) & between(.$UMAP_2, -1,9)) %>% summary

acc[['RNA']][acc[['spca']][,'SPARSE_11'] %>% subset(. != 0) %>% abs %>% as.data.frame %>% arrange(desc(SPARSE_11)) %>% rownames,] %>% t %>% as.matrix %>% subset(with(acc[['umap.spca']]@cell.embeddings %>% as.data.frame, between(UMAP_1, -3,5.5) & between(UMAP_2, -1,9))) %>% cor
acc[['RNA']][acc[['spca']][,'SPARSE_14'] %>% subset(. != 0) %>% abs %>% as.data.frame %>% arrange(desc(SPARSE_14)) %>% rownames,] %>% t %>% as.matrix %>% subset(with(acc[['umap.spca']]@cell.embeddings %>% as.data.frame, between(UMAP_1, -3,5.5) & between(UMAP_2, -1,9))) %>% cor

acc.label <- acc
Idents(acc.label) <- "unknown"
acc.label <- acc.label %>% CellSelector(DimPlot(., red='umap.spca'), ., "malignant")
acc.label <- acc.label %>% CellSelector(DimPlot(., red='umap.spca'), ., "tissue")
acc.label <- acc.label %>% CellSelector(DimPlot(., red='umap.spca'), ., "others1")
acc.label <- acc.label %>% CellSelector(DimPlot(., red='umap.spca'), ., "others2")
acc.label <- acc.label %>% CellSelector(DimPlot(., red='umap.spca'), ., "others3")
acc.label <- acc.label %>% CellSelector(DimPlot(., red='umap.spca'), ., "others4")
table(Idents(acc.label))
Idents(acc.label) <- Idents(acc.label) %>% fct_relevel(c('malignant','tissue','others1','others2','others3','others4'))

VlnPlot(acc.label, paste0('SPARSE_', 1:12))
ggsave('AccExample-SPCA-Violin-1.pdf', width=8, height=12)
VlnPlot(acc.label, paste0('SPARSE_', 13:24))
ggsave('AccExample-SPCA-Violin-2.pdf', width=8, height=12)

VlnPlot(acc.label, paste0('PC_', 1:12))
ggsave('AccExample-PCA-Violin.pdf', width=8, height=12)

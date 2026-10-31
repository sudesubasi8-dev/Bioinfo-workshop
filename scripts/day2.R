library(Seurat)
library(tidyverse)
library(patchwork)
counts <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19")
dim(counts)
class(counts)
counts[c("CD3D", "MS4A1", "LYZ"), 1:6]
format(object.size(counts), units = "MB")
sprintf("as an ordinary (dense) matrix: %.0f MB", prod(dim(counts)) * 8 / 1e6)
zero_share <- 1 - Matrix::nnzero(counts) / prod(dim(counts))
zero_share
pbmc <- CreateSeuratObject(counts = counts, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
head(pbmc[[]])
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
summary(pbmc[[]][, c("nCount_RNA", "nFeature_RNA", "percent.mt")])
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
meta <- pbmc[[]] |> as_tibble(rownames = "barcode")

meta |>
  pivot_longer(c(nFeature_RNA, nCount_RNA, percent.mt), names_to = "metric") |>
  ggplot(aes(value)) + geom_histogram(bins = 80) +
  facet_wrap(~ metric, scales = "free") +
  theme_bw()
quantile(meta$nFeature_RNA, c(0.01, 0.05, 0.5, 0.95, 0.99))
quantile(meta$percent.mt, c(0.5, 0.9, 0.95, 0.99))
thresholds <- expand_grid(max_features = c(2000, 2500, 3000, Inf), max_mt =
                            c(5, 10, 20)) |>
  mutate(cells_kept = map2_int(max_features, max_mt,
                               \(f, m) sum(meta$nFeature_RNA > 200 &
                                             meta$nFeature_RNA < f & meta$percent.mt < m)))
thresholds |> pivot_wider(names_from = max_mt, values_from = cells_kept,
                          names_prefix = "percent_mt_below_")
platelet_like <- LayerData(pbmc, layer = "counts")["PPBP", ] >= 10
tibble(min_features = c(200, 300, 400, 500)) |>
  mutate(cells_kept = map_int(min_features, \(m) sum(meta$nFeature_RNA > m &
                                                       meta$nFeature_RNA < 2500 & meta$percent.mt < 5)),
         platelet_like_kept = map_int(min_features, \(m) sum(platelet_like &
                                                               meta$nFeature_RNA > m)))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor
                      = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
top10
LabelPoints(plot = VariableFeaturePlot(pbmc), points = top10, repel = TRUE)
all_genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all_genes)
pbmc
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
print(pbmc[["pca"]], dims = 1:3, nfeatures = 5)
DimPlot(pbmc, reduction = "pca") + ElbowPlot(pbmc, ndims = 30)
DimHeatmap(pbmc, dims = 1:9, cells = 500, balanced = TRUE)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
table(Idents(pbmc))
res_table <- map_dfr(c(0.2, 0.5, 0.8, 1.2), \(r) {
  cl <- FindClusters(pbmc, resolution = r, verbose = FALSE)
  tibble(resolution = r, clusters = nlevels(Idents(cl)))
})
res_table
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()
pbmc
head(pbmc[[]], 3)
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
nrow(markers)
markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 4, with_ties = FALSE) |>
  summarise(top_markers = paste(gene, collapse = ", "))
canonical <- c("IL7R", "CCR7", "S100A4", "CD14", "LYZ", "MS4A1", "CD8A",
               "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP")
DotPlot(pbmc, features = canonical) + RotatedAxis()
marker_sets <- list(
  "Naive CD4 T" = c("IL7R", "CCR7", "LEF1"),
  "Memory CD4 T" = c("IL7R", "S100A4", "CD3E"),
  "CD14+ Mono" = c("CD14", "LYZ"),
  "B" = c("MS4A1", "CD79A"),
  "CD8 T" = c("CD8A", "CD8B"),
  "FCGR3A+ Mono" = c("FCGR3A", "MS4A7"),
  "NK" = c("GNLY", "NKG7", "KLRD1"),
  "DC" = c("FCER1A", "CST3"),
  "Platelet" = c("PPBP", "PF4")
)
avg <- AverageExpression(pbmc, features = unique(unlist(marker_sets)), layer =
                           "scale.data",
                         group.by = "seurat_clusters")$RNA
scores <- sapply(marker_sets, \(g) colMeans(avg[intersect(g, rownames(avg)), ,
                                                drop = FALSE]))
round(scores, 2)
assigned <- colnames(scores)[apply(scores, 1, which.max)]
names(assigned) <- rownames(scores)
assigned
names(assigned) <- sub("^g", "", names(assigned))
pbmc$cell_type <- unname(assigned[as.character(pbmc$seurat_clusters)])
Idents(pbmc) <- "cell_type"
DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
composition <- pbmc[[]] |>
  as_tibble(rownames = "barcode") |>
  count(cell_type, name = "cells") |>
  mutate(percent = round(100 * cells / sum(cells), 1)) |>
  arrange(desc(cells))
composition
ggplot(composition, aes(x = reorder(cell_type, cells), y = cells)) +
  geom_col() + coord_flip() +
  labs(x = NULL, y = "Cells", title = "Cell-type composition") + theme_bw()
markers_named <- markers |>
  as_tibble() |>
  mutate(cluster_cell_type = unname(assigned[as.character(cluster)])) |>
  relocate(cluster_cell_type, .after = cluster)
markers_named |> group_by(cluster_cell_type) |> slice_max(avg_log2FC, n = 2,
                                                          with_ties = FALSE) |>
  select(cluster_cell_type, gene, avg_log2FC, pct.1, pct.2, p_val_adj)
write_tsv(markers_named, "results/pbmc3k_markers.tsv")
saveRDS(pbmc, "results/pbmc3k_annotated.rds")
file.size("results/pbmc3k_annotated.rds") / 1e6
pbmc2 <- readRDS("results/pbmc3k_annotated.rds")
identical(Idents(pbmc2), Idents(pbmc))
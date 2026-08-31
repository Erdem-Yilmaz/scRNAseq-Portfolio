# Load libraries ----

library(Seurat)
library(tidyverse)
library(cowplot)
library(SingleCellExperiment)
library(scDblFinder)
library(viridis)
library(SingleR)
library(celldex)
library(BPCells)


# Import data ----

counts_1  <- Read10X(data.dir = "GSE318371_RAW/hd1")
counts_2  <- Read10X(data.dir = "GSE318371_RAW/hd2")
counts_3  <- Read10X(data.dir = "GSE318371_RAW/hd3")
counts_4  <- Read10X(data.dir = "GSE318371_RAW/hd4")
counts_5  <- Read10X(data.dir = "GSE318371_RAW/hd5")
counts_6  <- Read10X(data.dir = "GSE318371_RAW/hd6")
counts_7  <- Read10X(data.dir = "GSE318371_RAW/hd7")
counts_8  <- Read10X(data.dir = "GSE318371_RAW/hd8")
counts_9  <- Read10X(data.dir = "GSE318371_RAW/hd9")
counts_10 <- Read10X(data.dir = "GSE318371_RAW/hd10")
counts_11 <- Read10X(data.dir = "GSE318371_RAW/hd11")
counts_12 <- Read10X(data.dir = "GSE318371_RAW/hd12")


# Check dimensions

dim(counts_1)
dim(counts_2)
dim(counts_3)
dim(counts_4)
dim(counts_5)
dim(counts_6)
dim(counts_7)
dim(counts_8)
dim(counts_9)
dim(counts_10)
dim(counts_11)
dim(counts_12)


# Create Seurat objects for each healthy donor ----

seu_1 <- CreateSeuratObject(
  counts = counts_1,
  min.features = 200
)

seu_2 <- CreateSeuratObject(
  counts = counts_2,
  min.features = 200
)

seu_3 <- CreateSeuratObject(
  counts = counts_3,
  min.features = 200
)

seu_4 <- CreateSeuratObject(
  counts = counts_4,
  min.features = 200
)

seu_5 <- CreateSeuratObject(
  counts = counts_5,
  min.features = 200
)

seu_6 <- CreateSeuratObject(
  counts = counts_6,
  min.features = 200
)

seu_7 <- CreateSeuratObject(
  counts = counts_7,
  min.features = 200
)

seu_8 <- CreateSeuratObject(
  counts = counts_8,
  min.features = 200
)

seu_9 <- CreateSeuratObject(
  counts = counts_9,
  min.features = 200
)

seu_10 <- CreateSeuratObject(
  counts = counts_10,
  min.features = 200
)

seu_11 <- CreateSeuratObject(
  counts = counts_11,
  min.features = 200
)

seu_12 <- CreateSeuratObject(
  counts = counts_12,
  min.features = 200
)

# Assign donor IDs

seu_1$donor_id <- "hd1"
seu_2$donor_id <- "hd2"
seu_3$donor_id <- "hd3"
seu_4$donor_id <- "hd4"
seu_5$donor_id <- "hd5"
seu_6$donor_id <- "hd6"
seu_7$donor_id <- "hd7"
seu_8$donor_id <- "hd8"
seu_9$donor_id <- "hd9"
seu_10$donor_id <- "hd10"
seu_11$donor_id <- "hd11"
seu_12$donor_id <- "hd12"


# Define donor IDs and donor colors

donor_levels <- c(
  "hd1",
  "hd2",
  "hd3",
  "hd4",
  "hd5",
  "hd6",
  "hd7",
  "hd8",
  "hd9",
  "hd10",
  "hd11",
  "hd12"
)

donor_colors <- viridis::viridis(
  12,
  option = "viridis"
)

names(donor_colors) <- donor_levels


# Merge all healthy donors into one Seurat object ----

seurat <- merge(
  x = seu_1,
  y = c(
    seu_2,
    seu_3,
    seu_4,
    seu_5,
    seu_6,
    seu_7,
    seu_8,
    seu_9,
    seu_10,
    seu_11,
    seu_12
  ),
  add.cell.ids = donor_levels
)

# Set donor ID as an ordered factor

seurat$donor_id <- factor(
  seurat$donor_id,
  levels = donor_levels
)

# Check the merged Seurat object
seurat
Layers(seurat[["RNA"]])
table(seurat$donor_id)
dim(seurat)

# Save initial merged Seurat object

saveRDS(
  seurat,
  file = "objects/healthy_seurat_01_merged.rds"
)


# QC metrics: nFeature_RNA distribution ----

p_nFeature_before <- VlnPlot(
  seurat,
  features = "nFeature_RNA",
  group.by = "donor_id",
  pt.size = 0
) +
  scale_fill_manual(
    values = donor_colors
  ) +
  theme_classic() +
  labs(
    title = "Number of Detected Genes per Sample Before Filtering",
    x = "Sample",
    y = "Number of Detected Genes"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_nFeature_before

# QC metrics: nCount_RNA distribution ----

p_nCount_before <- VlnPlot(
  seurat,
  features = "nCount_RNA",
  group.by = "donor_id",
  pt.size = 0
) +
  scale_fill_manual(
    values = donor_colors
  ) +
  theme_classic() +
  labs(
    title = "RNA Counts per Sample Before Filtering",
    x = "Sample",
    y = "Total RNA Counts"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_nCount_before


# QC metrics: Mitochondrial RNA percentage ----

seurat[["percent.mt"]] <- PercentageFeatureSet(
  seurat,
  pattern = "^MT-"
)

p_percent_mt_before <- VlnPlot(
  seurat,
  features = "percent.mt",
  group.by = "donor_id",
  pt.size = 0
) +
  scale_fill_manual(
    values = donor_colors
  ) +
  theme_classic() +
  labs(
    title = "Mitochondrial RNA Percentage per Sample Before Filtering",
    x = "Sample",
    y = "Mitochondrial RNA (%)"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_percent_mt_before


# QC metrics: Scatter nCount_RNA versus nFeature_RNA before filtering ----

p_counts_features_before <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "donor_id",
  raster = TRUE
) +
  scale_color_manual(
    values = donor_colors
  ) +
  theme_classic() +
  labs(
    title = "RNA Counts Versus Detected Genes Before Filtering",
    x = "Detected Genes",
    y = "Total RNA",
    color = "Sample"
  )

p_counts_features_before


# QC metrics: Scatter nCount_RNA versus percent.mt before filtering ----

p_counts_mt_before <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  group.by = "donor_id",
  raster = TRUE
) +
  scale_color_manual(
    values = donor_colors
  ) +
  theme_classic() +
  labs(
    title = "RNA Counts Versus Mitochondrial RNA Before Filtering",
    x = "Total RNA",
    y = "Mitochondrial RNA (%)",
    color = "Sample"
  )

p_counts_mt_before


# QC metrics: Scatter nFeature_RNA versus percent.mt before filtering ----

p_features_mt_before <- FeatureScatter(
  seurat,
  feature1 = "nFeature_RNA",
  feature2 = "percent.mt",
  group.by = "donor_id",
  raster = TRUE
) +
  scale_color_manual(
    values = donor_colors
  ) +
  theme_classic() +
  labs(
    title = "Detected Genes Versus Mitochondrial RNA Before Filtering",
    x = "Detected Genes",
    y = "Mitochondrial RNA (%)",
    color = "Sample"
  )

p_features_mt_before


# Combine QC plots before filtering ----

qc_plots_aligned <- align_plots(
  p_nFeature_before,
  p_nCount_before,
  p_percent_mt_before,
  p_counts_features_before,
  p_counts_mt_before,
  p_features_mt_before,
  align = "hv",
  axis = "tblr"
)

# Combine aligned QC plots

qc_before <- plot_grid(
  plotlist = qc_plots_aligned,
  ncol = 2,
  nrow = 3,
  labels = c(
    "A", "B", "C",
    "D", "E", "F"
  )
)

# Add figure title

qc_before <- ggdraw() +
  draw_label(
    "QC Metrics Before Filtering",
    fontface = "bold",
    size = 18,
    x = 0.5,
    y = 0.98
  ) +
  draw_plot(
    qc_before,
    x = 0,
    y = 0,
    width = 1,
    height = 0.94
  )

qc_before

# Save QC figure

ggsave(
  "figures/Figure1_QC_Metrics_Before_Filtering.svg",
  plot = qc_before,
  width = 12,
  height = 16,
  units = "in"
)







# Inspect QC thresholds ----

qc_quantiles_nFeature <- quantile(
  seurat$nFeature_RNA,
  probs = c(
    0.01, 0.025, 0.05, 0.10,
    0.50, 0.90, 0.95, 0.975, 0.99
  )
)

qc_quantiles_nCount <- quantile(
  seurat$nCount_RNA,
  probs = c(
    0.01, 0.025, 0.05, 0.10,
    0.50, 0.90, 0.95, 0.975, 0.99
  )
)

qc_quantiles_percent_mt <- quantile(
  seurat$percent.mt,
  probs = c(
    0.50, 0.75, 0.90,
    0.95, 0.975, 0.99
  )
)

qc_quantiles_nFeature
qc_quantiles_nCount
qc_quantiles_percent_mt

# Filter cells based on QC metrics ----

seurat <- subset(
  seurat,
  subset =
    nFeature_RNA >= 500 &
    nCount_RNA >= 800 &
    percent.mt < 10
)

# Check remaining cells after QC filtering

dim(seurat)
table(seurat$donor_id)

# QC summaries after filtering ----

qc_summary_after_donor <- seurat@meta.data %>%
  group_by(donor_id) %>%
  summarise(
    cells = n(),
    median_nFeature_RNA = median(nFeature_RNA),
    median_nCount_RNA = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  )


# Plots after filtering ----

p_nFeature_after <- VlnPlot(
  seurat,
  features = "nFeature_RNA",
  group.by = "donor_id",
  pt.size = 0
) +
  scale_fill_manual(
    values = donor_colors
  ) +
  theme_classic() +
  labs(
    title = "Detected Genes After Filtering",
    x = "Sample",
    y = "nFeature_RNA"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_nCount_after <- VlnPlot(
  seurat,
  features = "nCount_RNA",
  group.by = "donor_id",
  pt.size = 0
) +
  scale_fill_manual(
    values = donor_colors
  ) +
  theme_classic() +
  labs(
    title = "RNA Counts After Filtering",
    x = "Sample",
    y = "nCount_RNA"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_percent_mt_after <- VlnPlot(
  seurat,
  features = "percent.mt",
  group.by = "donor_id",
  pt.size = 0
) +
  scale_fill_manual(
    values = donor_colors
  ) +
  theme_classic() +
  labs(
    title = "Mitochondrial RNA After Filtering",
    x = "Sample",
    y = "percent.mt"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# Scatter plots after filtering

p_counts_features_after <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "donor_id",
  raster = TRUE
) +
  scale_color_manual(
    values =donor_colors
  ) +
  theme_classic() +
  labs(
    title = "RNA Counts Versus Detected Genes After Filtering",
    x = "Detected Genes",
    y = "Total RNA",
    color = "Sample"
  )

p_counts_mt_after <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  group.by = "donor_id",
  raster = TRUE
) +
  scale_color_manual(
    values = donor_colors
  ) +
  theme_classic() +
  labs(
    title = "RNA Counts Versus Mitochondrial RNA After Filtering",
    x = "Total RNA",
    y = "Mitochondrial RNA (%)",
    color = "Sample"
  )

p_features_mt_after <- FeatureScatter(
  seurat,
  feature1 = "nFeature_RNA",
  feature2 = "percent.mt",
  group.by = "donor_id",
  raster = TRUE
) +
  scale_color_manual(
    values = donor_colors
  ) +
  theme_classic() +
  labs(
    title = "Detected Genes Versus Mitochondrial RNA After Filtering",
    x = "Detected Genes",
    y = "Mitochondrial RNA (%)",
    color = "Sample"
  )

# Align QC plots after filtering

qc_plots_after_aligned <- align_plots(
  p_nFeature_after,
  p_nCount_after,
  p_percent_mt_after,
  p_counts_features_after,
  p_counts_mt_after,
  p_features_mt_after,
  align = "hv",
  axis = "tblr"
)

# Combine QC plots after filtering

qc_after <- plot_grid(
  plotlist = qc_plots_after_aligned,
  ncol = 2,
  nrow = 3,
  labels = c(
    "A", "B", "C",
    "D", "E", "F"
  )
)

# Add figure title

qc_after <- ggdraw() +
  draw_label(
    "QC Metrics After Filtering",
    fontface = "bold",
    size = 18,
    x = 0.5,
    y = 0.98
  ) +
  draw_plot(
    qc_after,
    x = 0,
    y = 0,
    width = 1,
    height = 0.94
  )

qc_after

# Save QC figure

ggsave(
  "figures/Figure2_QC_Metrics_After_Filtering.svg",
  plot = qc_after,
  width = 12,
  height = 16,
  units = "in"
)



































# Save QC figure as SVG
ggplot2::ggsave(
  filename = "figures/Figure2_Healthy_QC_Metrics_After_Filtering.svg",
  plot = qc_after,
  width = 12,
  height = 16,
  units = "in"
)

# Clean up QC plot objects

rm(
  p_nFeature_before,
  p_nCount_before,
  p_percent_mt_before,
  p_counts_features_before,
  p_counts_mt_before,
  p_features_mt_before,
  qc_plots_before_aligned,
  qc_before
)

rm (qc_quantiles_nCount, qc_quantiles_nFeature, qc_quantiles_percent_mt)


rm(
  p_nFeature_after,
  p_nCount_after,
  p_percent_mt_after,
  p_counts_features_after,
  p_counts_mt_after,
  p_features_mt_after,
  qc_plots_after_aligned,
  qc_after
)

gc()


# Save filtered Seurat object

saveRDS(
  seurat,
  file = "objects/healthy_seurat_02_after_QC_filtering.rds"
)

# Join RNA layers for doublet detection ----

seurat[["RNA"]] <- JoinLayers(
  seurat[["RNA"]]
)
# Doublet detection with scDblFinder ----

sce <- as.SingleCellExperiment(
  seurat
)

# Run doublet detection

set.seed(1234)

sce <- scDblFinder(
  sce,
  samples = sce$donor_id,
  multiSampleMode = "split"
)

# Transfer doublet results to Seurat

seurat$scDblFinder.class <- sce$scDblFinder.class
seurat$scDblFinder.score <- sce$scDblFinder.score


# Check overall doublet classification

table(
  seurat$scDblFinder.class
)

# Check doublet classification by sample

table(
  seurat$donor_id,
  seurat$scDblFinder.class
)


# Check doublet classification by donor

table(
  seurat$donor_id,
  seurat$scDblFinder.class
)


# Doublet percentage by sample ----

doublet_percent_sample <- prop.table(
  table(
    seurat$donor_id,
    seurat$scDblFinder.class
  ),
  margin = 1
) * 100

doublet_percent_sample <- as.data.frame(
  doublet_percent_sample
)

colnames(doublet_percent_sample) <- c(
  "donor_id",
  "classification",
  "percent"
)

doublet_percent_sample$sample_id <- factor(
  doublet_percent_sample$donor_id,
  levels = c(
    "HD_1", "HD_2", "HD_3", "HD_4",
    "HD_5", "HD_6", "HD_7", "HD_8",
    "HD_9", "HD_10", "HD_11", "HD_12"
  )
)

# Save doublet summaries

write_csv(
  doublet_percent_sample,
  "summary_data/Doublet_Percentage_by_Sample.csv"
)

# Plot doublet percentages by sample

p_doublet_percent <- ggplot(
  doublet_percent_sample[
    doublet_percent_sample$classification == "doublet",
  ],
  aes(
    x = donor_id,
    y = percent,
    fill = donor_id
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = donor_colors
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +
  theme_classic() +
  labs(
    title = "Doublet Percentage by Sample",
    x = "Sample",
    y = "Doublets (%)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_doublet_percent


# Save doublet percentage figure

ggsave(
  "figures/Figure3_Doublet_Percentage_by_Sample.svg",
  plot = p_doublet_percent,
  width = 8,
  height = 5,
  units = "in"
)


# Save Seurat object after doublet detection

saveRDS(
  seurat,
  "objects/seurat_03_healthy_after_doublet_detection.rds"
)


# Remove predicted doublets ----

seurat <- subset(
  seurat,
  subset = scDblFinder.class == "singlet"
)


# Check remaining cells

dim(seurat)

# Check remaining cells by sample

table(
  seurat$donor_id
)

# Check that no doublets remain

table(
  seurat$scDblFinder.class
)


# Save filtered singlet-only Seurat object

saveRDS(
  seurat,
  "objects/seurat_04_healthy_after_doublet_removal.rds"
)


# Normalize data ----

seurat <- NormalizeData(
  seurat,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)


# Save normalized Seurat object

saveRDS(
  seurat,
  "objects/seurat_05_healthy_normalized.rds"
)

# Identify highly variable genes ----

seurat <- FindVariableFeatures(
  seurat,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)


# Check how many variable features were selected

length(
  VariableFeatures(seurat)
)


# Display the first few selected variable features

head(
  VariableFeatures(seurat)
)


# Save selected variable features

variable_features <- data.frame(
  gene = VariableFeatures(seurat)
)

write_csv(
  variable_features,
  "summary_data/healthy_variable_features_2000.csv"
)


# Scale the 2,000 highly variable genes ----

seurat <- ScaleData(
  seurat,
  assay = "RNA",
  features = VariableFeatures(seurat),
  verbose = TRUE
)


# Confirm that scale.data was created

Layers(
  seurat[["RNA"]]
)

saveRDS(
  seurat,
  "objects/seurat_06_healthy_scaled.rds"
)


# Run PCA using ----

seurat <- RunPCA(
  seurat,
  assay = "RNA",
  features = VariableFeatures(seurat),
  npcs = 50,
  verbose = TRUE
)


# Check that the PCA reduction was created

seurat[["pca"]]


# Check the dimensions of the PCA embedding

dim(
  Embeddings(
    seurat,
    "pca"
  )
)


# Show the first few cells and PCs

head(
  Embeddings(
    seurat,
    "pca"
  )
)


# Save PCA checkpoint

saveRDS(
  seurat,
  "objects/seurat_07_healthy_PCA.rds"
)


# PCA elbow plot ----

pca_elbow <- ElbowPlot(
  seurat,
  reduction = "pca",
  ndims = 50
) +
  theme_classic() +
  ggtitle(
    "PCA Elbow Plot"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5
    )
  )

# Display PCA elbow plot

pca_elbow


# Save PCA elbow plot

ggsave(
  "figures/Figure4_PCA_Elbow_Plot.svg",
  plot = pca_elbow,
  width = 8,
  height = 6,
  units = "in"
)


# Build the nearest-neighbor graph ----

seurat <- FindNeighbors(
  seurat,
  reduction = "pca",
  dims = 1:25,
  verbose = TRUE
)


# Check that the neighbor graph was created

names(
  seurat@graphs
)


# Cluster cells ----

seurat <- FindClusters(
  seurat,
  resolution = 0.4,
  verbose = TRUE
)

# Check how many cells are assigned to each cluster

table(
  Idents(seurat)
)


# Count the number of clusters

length(
  unique(
    Idents(seurat)
  )
)


# Save cluster sizes

cluster_sizes <- as.data.frame(
  table(
    seurat$seurat_clusters
  )
)

colnames(cluster_sizes) <- c(
  "cluster",
  "cells"
)

write_csv(
  cluster_sizes,
  "summary_data/healthy_cluster_sizes_resolution_0.4.csv"
)


# Run UMAP ----

seurat <- RunUMAP(
  seurat,
  reduction = "pca",
  dims = 1:25,
  verbose = TRUE
)


# Plot the UMAP colored by the current clusters

p_umap_clusters <- DimPlot(
  seurat,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  theme_classic() +
  ggtitle(
    "Healthy PBMC Clustering"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5
    )
  )


# Display clustered UMAP

p_umap_clusters


# Save clustered UMAP

ggsave(
  "figures/Figure5_Healthy_PBMC_Clustering_UMAP.svg",
  plot = p_umap_clusters,
  width = 10,
  height = 8,
  units = "in"
)


# Save clustering checkpoint

saveRDS(
  seurat,
  "objects/seurat_08_healthy_clustered_UMAP.rds"
)


# BPCells → Cluster Markers → SingleR Annotation ----

# Check the RNA assay class

class(
  seurat[["RNA"]]
)


# Check the current RNA layers

Layers(
  seurat,
  assay = "RNA"
)


# Check the class of the existing counts layer

class(
  seurat[["RNA"]]$counts
)


# Confirm the counts dimensions

dim(
  seurat[["RNA"]]$counts
)


# Convert the existing sparse counts matrix to BPCells ----
# Create BPCells parent directory

dir.create(
  "bpcells",
  showWarnings = FALSE
)

# Write RNA counts to BPCells

BPCells::write_matrix_dir(
  mat = seurat[["RNA"]]$counts,
  dir = "bpcells/"
)


# Open the BPCells matrix

rna_counts_bp <- BPCells::open_matrix_dir(
  dir = "bpcells/healthy_RNA_counts"
)

# Confirm BPCells matrix class

class(
  rna_counts_bp
)

# Confirm BPCells matrix dimensions

dim(
  rna_counts_bp
)

# Confirm rownames are preserved

head(
  rownames(rna_counts_bp)
)

# Confirm cell names are preserved

head(
  colnames(rna_counts_bp)
)


# Replace the in-memory counts layer with the BPCells-backed counts matrix

seurat[["RNA"]]$counts <- rna_counts_bp


# Confirm the counts layer class after conversion

class(
  seurat[["RNA"]]$counts
)

# Save

saveRDS(seurat, "objects/seurat_09_bpcells.rds")

# Find markers for every cluster ----

cluster_markers <- FindAllMarkers(
  seurat,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)


# Check the marker table

head(
  cluster_markers
)


# Check how many marker genes were detected for each cluster

table(
  cluster_markers$cluster
)


# Save complete cluster marker table

write_csv(
  cluster_markers,
  "summary_data/healthy_cluster_markers_all.csv"
)


# Select top 20 marker genes per cluster ----

top20_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup()


# Display the top markers for every cluster

top20_markers %>%
  select(
    cluster,
    gene,
    avg_log2FC,
    pct.1,
    pct.2
  ) %>%
  print(
    n = Inf
  )


# Save top 20 markers

write_csv(
  top20_markers,
  "summary_data/healthy_top20_markers_per_cluster.csv"
)


# Confirm donor and sample metadata

colnames(
  seurat@meta.data
)

table(
  seurat$donor_id
)

table(
  seurat$sample_id
)

table(
  seurat$donor_id,
  seurat$sample_id
)


# Save donor-by-sample structure

donor_sample_structure <- as.data.frame(
  table(
    donor_id = seurat$donor_id,
    sample_id = seurat$sample_id
  )
) %>%
  filter(
    Freq > 0
  )

write_csv(
  donor_sample_structure,
  "summary_data/healthy_donor_by_sample_structure.csv"
)


# SingleR ----

packageVersion(
  "SingleR"
)

packageVersion(
  "celldex"
)


# Check the class of the BPCells counts layer

class(
  LayerData(
    seurat,
    assay = "RNA",
    layer = "counts"
  )
)


# Confirm clusters are still present

table(
  Idents(seurat)
)


# Load Human Primary Cell Atlas reference

hpca_ref <- HumanPrimaryCellAtlasData()


# Check reference dimensions

dim(
  hpca_ref
)


# Check available reference cell types

table(
  hpca_ref$label.main
)


# Convert Seurat object to SingleCellExperiment

sce <- as.SingleCellExperiment(
  seurat,
  assay = "RNA"
)


# Identify genes shared with HPCA reference

common_genes <- intersect(
  rownames(sce),
  rownames(hpca_ref)
)


# Report number of shared genes

length(
  common_genes
)


# Save common gene information

common_gene_summary <- data.frame(
  common_genes = common_genes
)

write_csv(
  common_gene_summary,
  "summary_data/healthy_SingleR_common_genes.csv"
)


# Extract normalized expression for SingleR

singler_expr <- as(
  assay(
    sce,
    "logcounts"
  )[common_genes, ],
  "dgCMatrix"
)


# SingleR broad annotation ----

# Broad annotation is performed first to establish major cellular lineages.

singler_clusters_broad <- SingleR(
  test = singler_expr,
  ref = hpca_ref,
  labels = hpca_ref$label.main,
  clusters = sce$seurat_clusters
)


# Extract cluster IDs

cluster_ids <- rownames(
  singler_clusters_broad
)


# Extract broad predictions

broad_labels <- singler_clusters_broad$labels


# Create broad annotation table

broad_annotation <- data.frame(
  cluster = cluster_ids,
  broad = broad_labels,
  stringsAsFactors = FALSE,
  row.names = NULL
)


# Display broad annotation table

broad_annotation


# Confirm number of annotated clusters

nrow(
  broad_annotation
)


# Display number of clusters per broad cell type

table(
  broad_labels
)


# Save broad SingleR annotation

write_csv(
  broad_annotation,
  "summary_data/healthy_SingleR_broad_cluster_annotations.csv"
)


# Identify clusters belonging to each broad lineage ----

t_clusters <- broad_annotation$cluster[
  grepl(
    "T_cell|T cell",
    broad_annotation$broad,
    ignore.case = TRUE
  )
]

nk_clusters <- broad_annotation$cluster[
  grepl(
    "NK",
    broad_annotation$broad,
    ignore.case = TRUE
  )
]

monocyte_clusters <- broad_annotation$cluster[
  grepl(
    "Monocyte",
    broad_annotation$broad,
    ignore.case = TRUE
  )
]

b_clusters <- broad_annotation$cluster[
  grepl(
    "B_cell|B cell",
    broad_annotation$broad,
    ignore.case = TRUE
  )
]

platelet_clusters <- broad_annotation$cluster[
  grepl(
    "Platelet",
    broad_annotation$broad,
    ignore.case = TRUE
  )
]

cmp_clusters <- broad_annotation$cluster[
  grepl(
    "CMP",
    broad_annotation$broad,
    ignore.case = TRUE
  )
]


# Check clusters assigned to each broad lineage

t_clusters
nk_clusters
monocyte_clusters
b_clusters
platelet_clusters
cmp_clusters


# Save broad lineage assignments

broad_lineage_clusters <- bind_rows(
  data.frame(
    cluster = t_clusters,
    lineage = "T cells"
  ),
  data.frame(
    cluster = nk_clusters,
    lineage = "NK cells"
  ),
  data.frame(
    cluster = monocyte_clusters,
    lineage = "Monocytes"
  ),
  data.frame(
    cluster = b_clusters,
    lineage = "B cells"
  ),
  data.frame(
    cluster = cmp_clusters,
    lineage = "CMP"
  )
)

write_csv(
  broad_lineage_clusters,
  "summary_data/healthy_SingleR_broad_lineage_cluster_assignments.csv"
)


# Create cell indices for each broad population

t_cells <- colnames(sce)[
  sce$seurat_clusters %in% t_clusters
]

nk_cells <- colnames(sce)[
  sce$seurat_clusters %in% nk_clusters
]

monocyte_cells <- colnames(sce)[
  sce$seurat_clusters %in% monocyte_clusters
]

b_cells <- colnames(sce)[
  sce$seurat_clusters %in% b_clusters
]

platelet_cells <- colnames(sce)[
  sce$seurat_clusters %in% platelet_clusters
]

cmp_cells <- colnames(sce)[
  sce$seurat_clusters %in% cmp_clusters
]


# Create fine-level HPCA reference subsets ----

t_reference <- grepl(
  "^T_cell",
  hpca_ref$label.fine
)

nk_reference <- grepl(
  "^NK_cell",
  hpca_ref$label.fine
)

monocyte_reference <- grepl(
  "^Monocyte",
  hpca_ref$label.fine
)

b_reference <- grepl(
  "^B_cell",
  hpca_ref$label.fine
)

platelet_reference <- hpca_ref$label.fine == "Platelets"

cmp_reference <- hpca_ref$label.fine == "CMP"


# Fine SingleR — T cells

singler_fine_T <- SingleR(
  test = singler_expr[
    ,
    t_cells,
    drop = FALSE
  ],
  ref = hpca_ref[
    ,
    t_reference
  ],
  labels = hpca_ref$label.fine[
    t_reference
  ],
  clusters = sce$seurat_clusters[
    colnames(sce) %in% t_cells
  ]
)


# Fine SingleR — NK cells

singler_fine_NK <- SingleR(
  test = singler_expr[
    ,
    nk_cells,
    drop = FALSE
  ],
  ref = hpca_ref[
    ,
    nk_reference
  ],
  labels = hpca_ref$label.fine[
    nk_reference
  ],
  clusters = sce$seurat_clusters[
    colnames(sce) %in% nk_cells
  ]
)


# Fine SingleR — monocytes

singler_fine_monocyte <- SingleR(
  test = singler_expr[
    ,
    monocyte_cells,
    drop = FALSE
  ],
  ref = hpca_ref[
    ,
    monocyte_reference
  ],
  labels = hpca_ref$label.fine[
    monocyte_reference
  ],
  clusters = sce$seurat_clusters[
    colnames(sce) %in% monocyte_cells
  ]
)


# Fine SingleR — B cells

singler_fine_B <- SingleR(
  test = singler_expr[
    ,
    b_cells,
    drop = FALSE
  ],
  ref = hpca_ref[
    ,
    b_reference
  ],
  labels = hpca_ref$label.fine[
    b_reference
  ],
  clusters = sce$seurat_clusters[
    colnames(sce) %in% b_cells
  ]
)


# Fine SingleR — platelets

singler_fine_platelets <- SingleR(
  test = singler_expr[
    ,
    platelet_cells,
    drop = FALSE
  ],
  ref = hpca_ref[
    ,
    platelet_reference
  ],
  labels = hpca_ref$label.fine[
    platelet_reference
  ],
  clusters = sce$seurat_clusters[
    colnames(sce) %in% platelet_cells
  ]
)


# Fine SingleR — CMP

singler_fine_CMP <- SingleR(
  test = singler_expr[
    ,
    cmp_cells,
    drop = FALSE
  ],
  ref = hpca_ref[
    ,
    cmp_reference
  ],
  labels = hpca_ref$label.fine[
    cmp_reference
  ],
  clusters = sce$seurat_clusters[
    colnames(sce) %in% cmp_cells
  ]
)


# Extract fine labels

fine_T <- singler_fine_T$labels

fine_NK <- singler_fine_NK$labels

fine_monocyte <- singler_fine_monocyte$labels

fine_B <- singler_fine_B$labels

fine_platelets <- singler_fine_platelets$labels

fine_CMP <- singler_fine_CMP$labels


# Create final hierarchical labels

fine_labels_hierarchical <- broad_labels


# Replace T-cell clusters

fine_labels_hierarchical[
  match(
    rownames(singler_fine_T),
    cluster_ids
  )
] <- fine_T


# Replace NK-cell clusters

fine_labels_hierarchical[
  match(
    rownames(singler_fine_NK),
    cluster_ids
  )
] <- fine_NK


# Replace monocyte clusters

fine_labels_hierarchical[
  match(
    rownames(singler_fine_monocyte),
    cluster_ids
  )
] <- fine_monocyte


# Replace B-cell clusters

fine_labels_hierarchical[
  match(
    rownames(singler_fine_B),
    cluster_ids
  )
] <- fine_B


# Replace platelet clusters

fine_labels_hierarchical[
  match(
    rownames(singler_fine_platelets),
    cluster_ids
  )
] <- fine_platelets


# Replace CMP clusters

fine_labels_hierarchical[
  match(
    rownames(singler_fine_CMP),
    cluster_ids
  )
] <- fine_CMP


# Final hierarchical annotation table ----

hierarchical_annotation <- data.frame(
  cluster = cluster_ids,
  broad = broad_labels,
  fine = unname(
    fine_labels_hierarchical
  ),
  stringsAsFactors = FALSE,
  row.names = NULL
)


# Display final hierarchical annotation

hierarchical_annotation


# Save hierarchical annotation table

write_csv(
  hierarchical_annotation,
  "summary_data/healthy_SingleR_hierarchical_cluster_annotations.csv"
)


# Add SingleR annotations to Seurat metadata ----

seurat$SingleR_broad <- broad_labels[
  match(
    as.character(
      seurat$seurat_clusters
    ),
    cluster_ids
  )
]

seurat$SingleR_fine <- fine_labels_hierarchical[
  match(
    as.character(
      seurat$seurat_clusters
    ),
    cluster_ids
  )
]


# Verify broad annotations

table(
  seurat$SingleR_broad
)


# Verify fine annotations

table(
  seurat$SingleR_fine
)


# Verify final cluster-level annotation table

hierarchical_annotation


# Create clean display labels ----

seurat$SingleR_fine_clean <- dplyr::recode(
  seurat$SingleR_fine,
  "T_cell:CD4+_Naive" = "CD4 Naive T",
  "T_cell:CD4+_central_memory" = "CD4 Central Memory T",
  "T_cell:CD4+_effector_memory" = "CD4 Effector Memory T",
  "T_cell:CD8+" = "CD8 T",
  "T_cell:CD8+_effector_memory_RA" = "CD8 TEMRA",
  "B_cell:Naive" = "Naive B",
  "B_cell:Memory" = "Memory B",
  "B_cell:Plasma_cell" = "Plasma B",
  "Monocyte:CD14+" = "CD14+ Mono",
  "Monocyte:CD16+" = "CD16+ Mono",
  "Monocyte:CD16-" = "CD16- Mono",
  "NK_cell" = "NK",
  "Platelets" = "Platelets",
  "CMP" = "CMP"
)


# Check clean labels

table(
  seurat$SingleR_fine_clean
)


# Save annotation summary by cell

annotation_summary <- seurat@meta.data %>%
  select(
    donor_id,
    sample_id,
    orig.ident,
    seurat_clusters,
    SingleR_broad,
    SingleR_fine,
    SingleR_fine_clean
  )

write_csv(
  annotation_summary,
  "summary_data/healthy_SingleR_cell_level_annotations.csv"
)


# Plot broad SingleR annotation ----

p_singler_broad <- DimPlot(
  seurat,
  reduction = "umap",
  group.by = "SingleR_broad",
  label = TRUE,
  repel = TRUE
) +
  theme_classic() +
  ggplot2::ggtitle(
    "SingleR — Broad Annotation"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5
    )
  )


# Display broad annotation UMAP

p_singler_broad


# Plot clean hierarchical fine annotation

p_singler_fine_clean <- DimPlot(
  seurat,
  reduction = "umap",
  group.by = "SingleR_fine_clean",
  label = TRUE,
  repel = TRUE
) +
  theme_classic() +
  ggplot2::ggtitle(
    "SingleR — Hierarchical Fine Annotation"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5
    )
  )


# Display clean fine annotation UMAP

p_singler_fine_clean


# Save SingleR UMAPs

ggsave(
  "figures/Figure6_SingleR_Broad_Annotation_UMAP.svg",
  plot = p_singler_broad,
  width = 10,
  height = 8,
  units = "in"
)

ggsave(
  "figures/Figure7_SingleR_Fine_Annotation_UMAP.svg",
  plot = p_singler_fine_clean,
  width = 10,
  height = 8,
  units = "in"
)


# Save annotated Seurat checkpoint

saveRDS(
  seurat,
  "objects/seurat_09_healthy_SingleR_annotated.rds"
)


# Combine broad and fine SingleR UMAPs

p_singler_combined <- cowplot::plot_grid(
  p_singler_broad,
  p_singler_fine_clean,
  ncol = 2,
  labels = c(
    "A",
    "B"
  )
)


# Display combined figure

p_singler_combined


# Save combined figure

ggsave(
  "figures/Figure8_SingleR_Broad_Fine_UMAP.svg",
  plot = p_singler_combined,
  width = 16,
  height = 7,
  units = "in"
)


# Canonical lineage markers ----

lineage_markers <- c(
  "CD3D",
  "CD3E",
  "TRBC1",
  "TRBC2",
  "NKG7",
  "GNLY",
  "CD79A",
  "MS4A1",
  "CD74",
  "HLA-DRA",
  "LYZ",
  "LST1",
  "S100A8",
  "S100A9",
  "FCN1",
  "PPBP",
  "PF4"
)


# Keep only markers present in the Seurat object

lineage_markers_present <- lineage_markers[
  lineage_markers %in% rownames(seurat)
]


# Save canonical markers used

write_csv(
  data.frame(
    gene = lineage_markers_present
  ),
  "summary_data/healthy_canonical_lineage_markers_present.csv"
)


# Plot canonical lineage markers across clusters

p_lineage_markers <- DotPlot(
  seurat,
  features = lineage_markers_present,
  group.by = "seurat_clusters"
) +
  Seurat::RotatedAxis() +
  ggplot2::ggtitle(
    "Canonical Lineage Markers by Cluster"
  )


# Display marker validation plot

p_lineage_markers


# Save canonical lineage marker validation plot ----

ggplot2::ggsave(
  "figures/Figure9_Canonical_Lineage_Markers.svg",
  plot = p_lineage_markers,
  width = 14,
  height = 7,
  units = "in"
)


# Broad lineage markers ----

# CD3D = T cells
# NKG7 = NK cells
# MS4A1 = B cells
# LYZ = Monocytes
# PPBP = Platelets
# CD34 = CMP

broad_markers <- c(
  "CD3D",
  "NKG7",
  "MS4A1",
  "LYZ",
  "PPBP",
  "CD34"
)


# Keep only broad markers present in dataset

broad_markers_present <- broad_markers[
  broad_markers %in% rownames(seurat)
]


# Create broad lineage marker UMAP

p_umap_broad_markers <- FeaturePlot(
  seurat,
  features = broad_markers_present,
  reduction = "umap",
  ncol = 3,
  min.cutoff = "q10",
  max.cutoff = "q90",
  cols = c(
    "grey90",
    "red"
  )
)


# Remove axis labels

p_umap_broad_markers <- p_umap_broad_markers +
  ggplot2::labs(
    x = NULL,
    y = NULL
  )


# Change existing facet labels ----

p_umap_broad_markers$facet$params$labeller <- ggplot2::as_labeller(
  c(
    CD3D = "T cells — CD3D",
    NKG7 = "NK cells — NKG7",
    MS4A1 = "B cells — MS4A1",
    LYZ = "Monocytes — LYZ",
    PPBP = "Platelets — PPBP",
    CD34 = "CMP — CD34"
  )
)


# Display marker UMAP ----

p_umap_broad_markers


# Save broad lineage marker UMAP

ggsave(
  "figures/Figure10_Broad_Lineage_Marker_UMAP.svg",
  plot = p_umap_broad_markers,
  width = 12,
  height = 8,
  units = "in"
)


# Combine SingleR broad annotation and marker UMAP

p_broad_combined <- cowplot::plot_grid(
  p_singler_broad,
  p_umap_broad_markers,
  ncol = 1,
  labels = c(
    "A",
    "B"
  )
)


# Display combined broad annotation figure

p_broad_combined


# Save combined broad annotation figure

ggsave(
  "figures/Figure11_Broad_Annotation_Validation.svg",
  plot = p_broad_combined,
  width = 10,
  height = 14,
  units = "in"
)

# Save final annotated Seurat object

saveRDS(
  seurat,
  "objects/seurat_10_Healthy_SingleR_annotated_final.rds"
)

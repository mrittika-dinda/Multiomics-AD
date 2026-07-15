#!/usr/bin/env Rscript
# 01_snatac_preprocessing.R
# Single-nucleus ATAC-seq Preprocessing and Clustering
# Dataset: GSE174367 (32 samples: 16 AD, 16 Control)
# Author: Mrittika Dinda

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomeInfoDb)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 01_snatac_preprocessing.R <peak_matrix.h5> <metadata.csv>")
}

peak_file <- args[1]
meta_file <- args[2]

# Load data
message("Loading snATAC-seq data...")
counts <- Read10X_h5(peak_file)
metadata <- read.csv(meta_file, row.names = 1)

# Create ChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "hg38",
  min.cells = 10,
  min.features = 200
)

snatac <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)
message(sprintf("Loaded %d cells with %d peaks", ncol(snatac), nrow(snatac)))

# QC
message("Calculating QC metrics...")
snatac <- NucleosomeSignal(object = snatac)
snatac <- TSSEnrichment(object = snatac, fast = FALSE)
snatac$blacklist_ratio <- FractionCountsInRegion(
  object = snatac, assay = "peaks", regions = blacklist_hg38)
snatac$pct_reads_in_peaks <- snatac$peaks_count / snatac$nCount_peaks * 100

# Filter
snatac <- subset(x = snatac,
  subset = nCount_peaks > 1000 & nCount_peaks < 50000 &
    pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 &
    nucleosome_signal < 4 & TSS.enrichment > 2)

message(sprintf("Cells after QC: %d", ncol(snatac)))

# Normalization & Dimensionality Reduction
snatac <- RunTFIDF(snatac)
snatac <- FindTopFeatures(snatac, min.cutoff = "q5")
snatac <- RunSVD(snatac)
snatac <- RunUMAP(snatac, reduction = "lsi", dims = 2:30)

# Clustering
snatac <- FindNeighbors(snatac, reduction = "lsi", dims = 2:30)
snatac <- FindClusters(snatac, algorithm = 3, resolution = 0.5)

# Visualization
dir.create("results/atacseq", recursive = TRUE, showWarnings = FALSE)

DimPlot(snatac, reduction = "umap", group.by = "seurat_clusters") +
  theme_minimal() + labs(title = "snATAC-seq Clusters")
ggsave("results/atacseq/umap_clusters.png", width = 8, height = 7, dpi = 300)

DimPlot(snatac, reduction = "umap", group.by = "condition") +
  scale_color_manual(values = c("Control" = "#377EB8", "AD" = "#E41A1C")) +
  theme_minimal() + labs(title = "snATAC-seq: AD vs Control")
ggsave("results/atacseq/umap_condition.png", width = 8, height = 7, dpi = 300)

# Save
saveRDS(snatac, "results/atacseq/snatac_processed.rds")
message("snATAC-seq preprocessing complete!")

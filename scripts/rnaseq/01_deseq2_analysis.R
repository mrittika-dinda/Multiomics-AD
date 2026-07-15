#!/usr/bin/env Rscript
# 01_deseq2_analysis.R
# Differential Gene Expression Analysis using DESeq2
# Dataset: GSE125583 (210 samples: 110 AD, 100 Control)
# Author: Mrittika Dinda

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
})

# Parse command line args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 01_deseq2_analysis.R <metadata_file> [count_matrix_file]")
}

metadata_file <- args[1]
count_file <- ifelse(length(args) >= 2, args[2], "data/rnaseq/count_matrix.csv")

# Load data
message("Loading sample metadata...")
metadata <- read.csv(metadata_file, row.names = 1)
metadata$condition <- factor(metadata$condition, levels = c("Control", "AD"))

message("Loading count matrix...")
counts <- read.csv(count_file, row.names = 1)
counts <- counts[, rownames(metadata)]

message(sprintf("Samples: %d | Genes: %d", ncol(counts), nrow(counts)))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = metadata,
  design = ~ condition
)

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run DESeq2
message("Running DESeq2...")
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "AD", "Control"))

# Thresholding
res_df <- as.data.frame(res) %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange >= 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange <= -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

sig_genes <- res_df %>% filter(sig != "Not Significant")
message(sprintf("Significant DEGs: %d (Up: %d, Down: %d)", 
                nrow(sig_genes), 
                sum(sig_genes$sig == "Upregulated"),
                sum(sig_genes$sig == "Downregulated")))

# Save results
dir.create("results/rnaseq", recursive = TRUE, showWarnings = FALSE)
write.csv(res_df, "results/rnaseq/deseq2_results_all.csv")
write.csv(sig_genes, "results/rnaseq/degs_significant.csv")

# Volcano plot
volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Upregulated" = "#E41A1C", 
                                  "Downregulated" = "#377EB8", 
                                  "Not Significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal(base_size = 14) +
  labs(title = "Differential Gene Expression: AD vs Control",
       subtitle = sprintf("GSE125583 | %d DEGs", nrow(sig_genes)),
       x = expression(Log[2]~Fold~Change),
       y = expression(-Log[10]~P-value),
       color = "Regulation")

ggsave("results/rnaseq/volcano_plot.png", volcano, width = 10, height = 8, dpi = 300)

# Heatmap of top DEGs
top_genes <- sig_genes %>% arrange(padj) %>% head(50) %>% rownames()
mat <- counts(dds, normalized = TRUE)[top_genes, ]
mat_z <- t(scale(t(mat)))

annotation_col <- data.frame(Condition = metadata$condition, row.names = rownames(metadata))

pheatmap(mat_z, annotation_col = annotation_col,
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdYlBu")))(100),
         cluster_cols = TRUE, cluster_rows = TRUE,
         show_rownames = TRUE, show_colnames = FALSE,
         fontsize_row = 6,
         main = "Top 50 DEGs: AD vs Control",
         filename = "results/rnaseq/heatmap_top50_degs.png",
         width = 12, height = 10, dpi = 300)

print("RNA-seq analysis complete!")

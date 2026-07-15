#!/usr/bin/env Rscript
# 01_deg_dmr_integration.R
# Two-Layer Integration: DEGs + Promoter DMRs
# Author: Mrittika Dinda

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(clusterProfiler)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 01_deg_dmr_integration.R <deg_results.csv> <dmr_promoter_results.csv>")
}

# Load DEGs
degs <- read.csv(args[1], row.names = 1)
sig_degs <- degs %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1) %>%
  mutate(expression_direction = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"))

# Convert Ensembl to gene symbols
ensembl_ids <- rownames(sig_degs)
ensembl_ids_clean <- gsub("\\.[0-9]+$", "", ensembl_ids)
gene_symbols <- bitr(ensembl_ids_clean, fromType = "ENSEMBL", 
                     toType = "SYMBOL", OrgDb = org.Hs.eg.db)
sig_degs$ENSEMBL <- ensembl_ids_clean
sig_degs <- merge(sig_degs, gene_symbols, by = "ENSEMBL", all.x = TRUE)

# Load DMRs
dmrs <- read.csv(args[2]) %>%
  mutate(methylation_direction = case_when(
    methylation_diff > 0.1 ~ "Hypermethylated",
    methylation_diff < -0.1 ~ "Hypomethylated",
    TRUE ~ "Not_Significant"
  )) %>% filter(methylation_direction != "Not_Significant")

# Integrate
integrated <- merge(sig_degs, dmrs, by = "SYMBOL", all = FALSE) %>%
  mutate(regulatory_pattern = case_when(
    expression_direction == "Upregulated" & methylation_direction == "Hypomethylated" ~ "Hypo-Up",
    expression_direction == "Downregulated" & methylation_direction == "Hypermethylated" ~ "Hyper-Down",
    expression_direction == "Upregulated" & methylation_direction == "Hypermethylated" ~ "Hyper-Up",
    expression_direction == "Downregulated" & methylation_direction == "Hypomethylated" ~ "Hypo-Down",
    TRUE ~ "Discordant"
  ))

message("Integration Results:")
print(integrated %>% count(regulatory_pattern))

# Save
dir.create("results/integration", recursive = TRUE, showWarnings = FALSE)
write.csv(integrated, "results/integration/deg_dmr_integrated.csv", row.names = FALSE)

# Quadrant plot
plot_data <- integrated %>%
  mutate(color_group = case_when(
    regulatory_pattern == "Hypo-Up" ~ "Hypo-Up",
    regulatory_pattern == "Hyper-Down" ~ "Hyper-Down",
    TRUE ~ "Other"
  ))

quadrant_plot <- ggplot(plot_data, aes(x = methylation_diff, y = log2FoldChange, color = color_group)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Hypo-Up" = "#377EB8", "Hyper-Down" = "#E41A1C", "Other" = "grey60")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "Integrated DEG-DMR Epigenetic Regulation",
       x = "DNA Methylation Difference (delta beta)",
       y = "Gene Expression Log2 Fold Change",
       color = "Regulatory Pattern")

ggsave("results/integration/quadrant_plot_deg_dmr.png", quadrant_plot, width = 10, height = 8, dpi = 300)

message("Two-layer integration complete!")

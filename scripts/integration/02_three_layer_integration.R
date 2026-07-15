#!/usr/bin/env Rscript
# 02_three_layer_integration.R
# Three-Layer Integration: DEGs + DMRs + DARs
# Author: Mrittika Dinda

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 02_three_layer_integration.R <deg.csv> <dmr_promoter.csv> <dar_promoter.csv>")
}

# Load data
sig_degs <- read.csv(args[1], row.names = 1) %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1) %>%
  mutate(expression_direction = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"))

# Convert Ensembl to SYMBOL if needed
if (!"SYMBOL" %in% colnames(sig_degs)) {
  library(org.Hs.eg.db)
  library(clusterProfiler)
  ensembl_ids <- rownames(sig_degs)
  ensembl_ids_clean <- gsub("\\.[0-9]+$", "", ensembl_ids)
  gene_symbols <- bitr(ensembl_ids_clean, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  sig_degs$ENSEMBL <- ensembl_ids_clean
  sig_degs <- merge(sig_degs, gene_symbols, by = "ENSEMBL", all.x = TRUE)
}

deg_genes <- unique(na.omit(sig_degs$SYMBOL))
dmr_genes <- unique(na.omit(read.csv(args[2])$SYMBOL))
dar_genes <- unique(na.omit(read.csv(args[3])$SYMBOL))

# Three-way intersection
integrated_genes <- intersect(intersect(deg_genes, dmr_genes), dar_genes)
message(sprintf("Three-layer integrated genes: %d", length(integrated_genes)))

# Classify
classification <- data.frame(SYMBOL = integrated_genes, stringsAsFactors = FALSE)
classification <- merge(classification, sig_degs %>% dplyr::select(SYMBOL, log2FoldChange, expression_direction), by = "SYMBOL", all.x = TRUE)
classification <- merge(classification, read.csv(args[2]) %>% dplyr::select(SYMBOL, methylation_diff, methylation_direction), by = "SYMBOL", all.x = TRUE)
classification <- merge(classification, read.csv(args[3]) %>% dplyr::select(SYMBOL, avg_log2FC, chromatin_state), by = "SYMBOL", all.x = TRUE)

classification <- classification %>%
  mutate(
    regulatory_class = case_when(
      expression_direction == "Upregulated" & methylation_direction == "Hypomethylated" & chromatin_state == "Open" ~ "Activated",
      expression_direction == "Downregulated" & methylation_direction == "Hypermethylated" & chromatin_state == "Closed" ~ "Repressed",
      TRUE ~ "Mixed/Discordant"
    ),
    regulatory_signature = paste(
      ifelse(chromatin_state == "Open", "Open chromatin", "Closed chromatin"),
      ifelse(methylation_direction == "Hypomethylated", "Hypomethylation", "Hypermethylation"),
      ifelse(expression_direction == "Upregulated", "Upregulation", "Downregulation"),
      sep = " + "
    )
  )

message("Regulatory Classification:")
print(classification %>% count(regulatory_class))

# Save
dir.create("results/integration", recursive = TRUE, showWarnings = FALSE)
write.csv(classification, "results/integration/three_layer_integrated_genes.csv", row.names = FALSE)

# Print results
message("\n=== High-Confidence Regulatory Genes ===")
for (cls in unique(classification$regulatory_class)) {
  cls_genes <- classification %>% filter(regulatory_class == cls)
  message(sprintf("\n%s (%d genes):", cls, nrow(cls_genes)))
  for (i in 1:nrow(cls_genes)) {
    message(sprintf("  - %s: %s", cls_genes$SYMBOL[i], cls_genes$regulatory_signature[i]))
  }
}

message("\nThree-layer integration complete!")

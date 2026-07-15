#!/usr/bin/env Rscript
# 02_differential_accessibility.R
# Differential Chromatin Accessibility Analysis
# Author: Mrittika Dinda

suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 02_differential_accessibility.R <snatac_processed.rds>")
}

snatac <- readRDS(args[1])
Idents(snatac) <- "condition"

# Differential accessibility
message("Finding DARs...")
dars <- FindMarkers(object = snatac, ident.1 = "AD", ident.2 = "Control",
                    min.pct = 0.05, test.use = "LR", latent.vars = "nCount_peaks")

dars <- dars %>%
  mutate(direction = case_when(
    avg_log2FC > 0 & p_val_adj < 0.05 ~ "Open_Chromatin",
    avg_log2FC < 0 & p_val_adj < 0.05 ~ "Closed_Chromatin",
    TRUE ~ "Not_Significant"
  ))

sig_dars <- dars %>% filter(direction != "Not_Significant")
message(sprintf("Significant DARs: %d (Open: %d, Closed: %d)",
                nrow(sig_dars),
                sum(sig_dars$direction == "Open_Chromatin"),
                sum(sig_dars$direction == "Closed_Chromatin")))

# Annotate peaks
peak_coords <- StringToGRanges(rownames(sig_dars), sep = c(":", "-"))
peak_anno <- annotatePeak(peak_coords, tssRegion = c(-3000, 3000),
                          TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                          annoDb = "org.Hs.eg.db")

promoter_dars <- peak_anno@anno[
  grepl("Promoter", peak_anno@anno$annotation, ignore.case = TRUE)]

# Save
dir.create("results/atacseq", recursive = TRUE, showWarnings = FALSE)
write.csv(dars, "results/atacseq/dars_all.csv")
write.csv(sig_dars, "results/atacseq/dars_significant.csv")
write.csv(as.data.frame(promoter_dars), "results/atacseq/dars_promoter_associated.csv")

# Volcano plot
dar_volcano <- ggplot(dars, aes(x = avg_log2FC, y = -log10(p_val), color = direction)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Open_Chromatin" = "#E41A1C",
                                  "Closed_Chromatin" = "#377EB8",
                                  "Not_Significant" = "grey70")) +
  theme_minimal(base_size = 14) +
  labs(title = "Differential Chromatin Accessibility: AD vs Control",
       x = "Average Log2 Fold Change",
       y = "-Log10 P-value",
       color = "Chromatin State")

ggsave("results/atacseq/dar_volcano_plot.png", dar_volcano, width = 10, height = 8, dpi = 300)

message("Differential accessibility analysis complete!")

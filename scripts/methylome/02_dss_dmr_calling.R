#!/usr/bin/env Rscript
# 02_dss_dmr_calling.R
# Differential Methylation Analysis using DSS
# Dataset: GSE227194 (10 samples: 5 AD, 5 Control)
# Author: Mrittika Dinda

suppressPackageStartupMessages({
  library(DSS)
  library(annotatr)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 02_dss_dmr_calling.R <sample_info.csv>")
}

sample_info <- read.csv(args[1])

# Read methylation data
message("Reading methylation data...")
bs_obj <- makeBSseqData(
  lapply(sample_info$methylation_file, function(f) {
    dat <- data.table::fread(f, select = c(1, 2, 5, 6))
    colnames(dat) <- c("chr", "pos", "N", "X")
    dat$chr <- gsub("chr", "", dat$chr)
    dat
  }),
  sampleNames = sample_info$sample_id
)

# Differential methylation
message("Running DSS...")
dml_test <- DMLtest(bs_obj, 
                    group1 = sample_info$sample_id[sample_info$condition == "Control"],
                    group2 = sample_info$sample_id[sample_info$condition == "AD"])

message("Calling DMRs...")
dmrs <- callDMR(dml_test, p.threshold = 0.05, minCG = 3, minLen = 50, dis.merge = 100)

# Filter significant DMRs
dmrs_filtered <- dmrs %>%
  as.data.frame() %>%
  mutate(
    methylation_diff = diff.Methy,
    direction = case_when(
      methylation_diff > 0.1 ~ "Hypermethylated",
      methylation_diff < -0.1 ~ "Hypomethylated",
      TRUE ~ "Not Significant"
    )
  ) %>%
  filter(direction != "Not Significant")

message(sprintf("Significant DMRs: %d (Hyper: %d, Hypo: %d)",
                nrow(dmrs_filtered),
                sum(dmrs_filtered$direction == "Hypermethylated"),
                sum(dmrs_filtered$direction == "Hypomethylated")))

# Save results
dir.create("results/methylome", recursive = TRUE, showWarnings = FALSE)
write.csv(dmrs_filtered, "results/methylome/dmrs_significant.csv", row.names = FALSE)

# Visualization
meth_dist <- ggplot(dmrs_filtered, aes(x = methylation_diff, fill = direction)) +
  geom_histogram(bins = 100, alpha = 0.8) +
  scale_fill_manual(values = c("Hypermethylated" = "#E41A1C", "Hypomethylated" = "#377EB8")) +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution of Differential Methylation",
       subtitle = sprintf("GSE227194 | %d DMRs", nrow(dmrs_filtered)),
       x = "Methylation Difference (delta beta)",
       y = "Count",
       fill = "Direction") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

ggsave("results/methylome/methylation_distribution.png", meth_dist, width = 10, height = 7, dpi = 300)

message("Methylome analysis complete!")

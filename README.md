# Multiomics-AD
Multi-omics analysis of differential gene  expression regulation in human brain during  Alzheimer’s Disease (AD) : RNA-seq, methylome, and snATAC-seq integration

# Alzheimer's Disease Multi-omics Analysis Pipeline

[![Python](https://img.shields.io/badge/Python-3.8+-3776AB?logo=python)](https://python.org)
[![R](https://img.shields.io/badge/R-4.3+-276DC3?logo=r)](https://r-project.org)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **End-to-end integrative analysis of transcriptomic, epigenomic, and chromatin accessibility data to identify gene expression regulation mechanisms in Alzheimer's Disease prefrontal cortex.**

---

## Overview

This repository contains the complete analysis pipeline from my **Master's dissertation** at **IBAB Bangalore**, supervised by **Dr. Kshitish K Acharya** and **Dr. Shyam Sundar Rajagopalan**.

**Biological Question:** Are transcriptional changes in Alzheimer's Disease driven by coordinated alterations in DNA methylation and chromatin accessibility at gene promoters?

**Key Finding:** AD-associated transcriptional dysregulation is characterized by **widespread promoter hypomethylation and epigenetic erosion**, challenging the classical paradigm of hypermethylation-driven gene silencing. Only **4 high-confidence genes** showed coordinated three-layer regulation.

---

## Datasets

| Omics Layer | GEO ID | Samples | Tissue | Platform |
|-------------|--------|---------|--------|----------|
| **RNA-seq** | [GSE125583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125583) | 210 (110 AD, 100 Control) | Prefrontal Cortex | Bulk RNA-seq |
| **Methylome** | [GSE227194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227194) | 10 (5 AD, 5 Control) | Prefrontal Cortex | RRBS |
| **snATAC-seq** | [GSE174367](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367) | 32 (16 AD, 16 Control) | Prefrontal Cortex | 10x Genomics |

---

## Pipeline Architecture
Raw Data -> QC -> Alignment -> Quantification -> Differential Analysis -> Integration -> Visualization
plain

### Layer 1: Transcriptomics (RNA-seq)
- **Tool:** DESeq2
- **Output:** 577 DEGs (366 upregulated, 211 downregulated)

### Layer 2: Epigenomics (RRBS Methylome)
- **Tool:** Bismark + DSS
- **Output:** 170,921 significant DMRs

### Layer 3: Chromatin Accessibility (snATAC-seq)
- **Tool:** Cell Ranger ATAC -> Seurat/Signac
- **Output:** Differentially accessible regions (DARs)

### Integration
- **Two-layer:** DEG intersect DMR (promoter +/-3kb TSS) -> 39 genes
- **Three-layer:** DEG intersect DMR intersect DAR -> **4 high-confidence genes**

---

## Repository Structure

| Path | Description |
|------|-------------|
| `scripts/rnaseq/` | RNA-seq differential expression (DESeq2) |
| `scripts/methylome/` | RRBS methylation analysis (Bismark, DSS) |
| `scripts/atacseq/` | snATAC-seq analysis (Seurat, Signac) |
| `scripts/integration/` | Multi-omics integration scripts |
| `configs/` | SLURM configs, sample metadata templates |
| `docs/` | Methodology and results documentation |
| `figures/` | Example visualizations |
| `data/` | Data download instructions |

---

## Quick Start

### Prerequisites

**R Packages:**
```r
install.packages(c("DESeq2", "ggplot2", "pheatmap", "dplyr", "Signac", "Seurat"))
BiocManager::install(c("DSS", "Bismark", "annotatr", "rGREAT", "ChIPseeker", 
                       "GenomicRanges", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene"))
Running the Pipeline
bash
# 1. RNA-seq Differential Expression
Rscript scripts/rnaseq/01_deseq2_analysis.R --metadata configs/rnaseq_metadata.csv

# 2. Methylome Analysis
bash scripts/methylome/01_bismark_pipeline.sh
Rscript scripts/methylome/02_dss_dmr_calling.R

# 3. snATAC-seq Analysis
Rscript scripts/atacseq/01_snatac_preprocessing.R
Rscript scripts/atacseq/02_differential_accessibility.R

# 4. Multi-omics Integration
Rscript scripts/integration/01_deg_dmr_integration.R
Rscript scripts/integration/02_three_layer_integration.R
Key Results
Differential Gene Expression
577 significant DEGs identified
Top upregulated: GAD1/2, NPAS4, BDNF, SNAP25 (compensatory synaptic remodeling)
Top downregulated: MMP9, CXCR4, FOXO4 (impaired amyloid clearance, immune dysregulation)
DNA Methylation
170,921 significant DMRs (|delta beta| > 0.1)
54.7% hypomethylated vs 45.3% hypermethylated -> epigenetic erosion
Three-Layer Integration (DEG + DMR + DAR)
Table
Gene	Regulation	Biological Role
PCDHGC5	Activated (Open + Hypo + Up)	Synaptic adhesion
LINC01962	Activated (Open + Hypo + Up)	lncRNA, chromatin regulation
NR4A1-AS	Repressed (Closed + Hyper + Down)	Neuronal survival
LOC100132249	Mixed/Discordant	Alternative regulation

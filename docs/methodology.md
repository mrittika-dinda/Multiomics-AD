# Methodology

## Study Design

This study investigates gene expression regulation in Alzheimer's Disease (AD) through integrative multi-omics analysis of the human prefrontal cortex.

### Datasets

| Omics | GEO ID | Samples | Platform |
|-------|--------|---------|----------|
| RNA-seq | GSE125583 | 210 (110 AD, 100 Control) | Bulk RNA-seq |
| Methylome | GSE227194 | 10 (5 AD, 5 Control) | RRBS |
| snATAC-seq | GSE174367 | 32 (16 AD, 16 Control) | 10x Genomics |

### Brain Region: Prefrontal Cortex

Selected because it is one of the earliest and most severely affected regions in AD with available multi-omics datasets.

## Pipeline Overview

### 1. RNA-seq Analysis
- **QC:** FastQC + MultiQC
- **Trimming:** Trimmomatic
- **Alignment:** HISAT2 to GRCh38
- **Quantification:** featureCounts (GENCODE v43)
- **DE Analysis:** DESeq2 (padj &lt; 0.05, |log2FC| &gt;= 1)

### 2. RRBS Methylome Analysis
- **QC:** FastQC + MultiQC
- **Trimming:** Trim Galore (RRBS mode)
- **Alignment:** Bismark (Bowtie2) to bisulfite-converted GRCh38
- **Deduplication:** deduplicate_bismark
- **DMR Calling:** DSS (|delta beta| &gt; 0.1, FDR &lt; 0.05)
- **Annotation:** annotatr + TxDb.Hsapiens.UCSC.hg38.knownGene

### 3. snATAC-seq Analysis
- **Preprocessing:** Cell Ranger ATAC
- **QC:** Nucleosome signal, TSS enrichment, FRiP
- **Clustering:** Seurat/Signac (LSI + UMAP)
- **Differential Accessibility:** FindMarkers (LR test, padj &lt; 0.05)
- **Annotation:** ChIPseeker (promoter +/-3kb TSS)

### 4. Multi-Omics Integration

**Two-Layer (DEG + DMR):**
- Promoter regions: +/-3kb of TSS
- Gene symbol harmonization: Ensembl -&gt; HGNC
- Integration: DEG intersect promoter-DMR

**Three-Layer (DEG + DMR + DAR):**
- Classification:
  - **Activated:** Open + Hypomethylated + Upregulated
  - **Repressed:** Closed + Hypermethylated + Downregulated
  - **Mixed/Discordant:** All other combinations

## Key Findings

1. **577 DEGs** identified (366 up, 211 down)
2. **170,921 significant DMRs** (54.7% hypomethylated)
3. **39 genes** in two-layer integration
4. **4 high-confidence genes** in three-layer integration

#!/bin/bash
# 01_bismark_pipeline.sh
# RRBS Methylation Analysis Pipeline using Bismark
# Dataset: GSE227194 (10 samples: 5 AD, 5 Control)
# Author: Mrittika Dinda

set -e

# Configuration
GENOME_DIR="/path/to/GRCh38_bismark"
RAW_DIR="data/methylome/raw"
TRIM_DIR="results/methylome/trimmed"
ALIGN_DIR="results/methylome/aligned"
METHYL_DIR="results/methylome/methylation"
QC_DIR="results/methylome/qc"

mkdir -p "$TRIM_DIR" "$ALIGN_DIR" "$METHYL_DIR" "$QC_DIR"
THREADS=8

# Step 1: Quality Control
echo "[1/6] Running FastQC..."
fastqc -t $THREADS -o "$QC_DIR" ${RAW_DIR}/*.fastq.gz

# Step 2: Adapter Trimming
echo "[2/6] Trimming adapters with Trim Galore..."
for fq in ${RAW_DIR}/*_1.fastq.gz; do
    base=$(basename "$fq" _1.fastq.gz)
    if [ -f "${RAW_DIR}/${base}_2.fastq.gz" ]; then
        trim_galore --paired --rrbs --fastqc --cores $THREADS -o "$TRIM_DIR" \
            "${RAW_DIR}/${base}_1.fastq.gz" "${RAW_DIR}/${base}_2.fastq.gz"
    else
        trim_galore --rrbs --fastqc --cores $THREADS -o "$TRIM_DIR" "$fq"
    fi
done

# Step 3: MultiQC
echo "[3/6] MultiQC..."
multiqc "$QC_DIR" -o "$QC_DIR"

# Step 4: Alignment
echo "[4/6] Aligning with Bismark..."
for trim_fq in ${TRIM_DIR}/*_trimmed.fq.gz; do
    bismark --genome "$GENOME_DIR" --parallel $THREADS --output_dir "$ALIGN_DIR" "$trim_fq"
done

# Step 5: Deduplication
echo "[5/6] Removing PCR duplicates..."
for bam in ${ALIGN_DIR}/*.bam; do
    deduplicate_bismark --bam "$bam"
done

# Step 6: Methylation Extraction
echo "[6/6] Extracting methylation calls..."
for dedup_bam in ${ALIGN_DIR}/*_deduplicated.bam; do
    bismark_methylation_extractor --bedGraph --counts --comprehensive \
        --multicore $THREADS --output_dir "$METHYL_DIR" "$dedup_bam"
done

echo "[Done] Bismark pipeline complete!"

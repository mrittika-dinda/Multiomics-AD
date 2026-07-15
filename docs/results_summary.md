# Results Summary

## RNA-seq Differential Expression

- **Total DEGs: 577** (366 upregulated, 211 downregulated)

### Top Upregulated Genes
| Gene | log2FC | Function |
|------|--------|----------|
| PVALB | +1.66 | Ca2+ buffering in interneurons |
| GAD2 | +1.58 | GABA synthesis |
| NPAS4 | +1.55 | Activity-dependent plasticity |
| GAD1 | +1.37 | GABA synthesis |
| BDNF | +1.16 | Synaptic plasticity |
| SNAP25 | +1.03 | Synaptic vesicle docking |

### Top Downregulated Genes
| Gene | log2FC | Function |
|------|--------|----------|
| TREML4 | -1.77 | Innate immune activation |
| MMP9 | -1.47 | ECM remodeling, A-beta clearance |
| CXCR4 | -1.20 | Immune cell migration |
| FOXO4 | -1.00 | Oxidative stress response |

## DNA Methylation

- **Total DMRs: 170,921** (|delta beta| &gt; 0.1)
- **Hypermethylated:** 77,500 (45.3%)
- **Hypomethylated:** 93,421 (54.7%)

**Key Finding:** Genome-wide hypomethylation bias suggests **epigenetic erosion** rather than selective gene silencing.

## Two-Layer Integration (DEG + DMR)

| Pattern | Count |
|---------|-------|
| Hypo-Up | 18 |
| Hyper-Down | 5 |
| Hypo-Down | 16 |
| Hyper-Up | 0 |

## Three-Layer Integration (DEG + DMR + DAR)

| Gene | Class | Signature | Biological Role |
|------|-------|-----------|-----------------|
| PCDHGC5 | Activated | Open + Hypo + Up | Synaptic adhesion |
| LINC01962 | Activated | Open + Hypo + Up | lncRNA, chromatin regulation |
| NR4A1-AS | Repressed | Closed + Hyper + Down | Neuronal survival |
| LOC100132249 | Mixed | Discordant | Alternative regulation |

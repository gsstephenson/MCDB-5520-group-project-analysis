# AlphaGenome Variant Analysis

This directory contains AlphaGenome deep learning model predictions for FCGR2B and GBP1 gene variants.

## Files

- `variant_analysis_FCGR2B_GBP1.ipynb` - Main analysis notebook with predictions for 6 variants across 2 genes
- `variant_summary.csv` - Variant metadata (positions, annotations, effects)
- `visualization_modality_tour.ipynb` - Reference example from AlphaGenome documentation
- `variant_analysis_results/` - All generated plots (14 total)

## Analysis Overview

### Genes Analyzed
- **GBP1** (Guanylate Binding Protein 1) - Immune response, GTPase activity
- **FCGR2B** (Fc Fragment Of IgG Receptor IIb) - Immune regulation

### Variants Analyzed (6 total)
**GBP1:**
1. chr1:89013380:G>A (p.Arg225Trp) - Missense in GTPase domain
2. chr1:89014075:T>C (c.625+8A>G) - Splice region variant

**FCGR2B:**
1. chr1:161506414:C>T (p.Gln63*) - **STOP-GAINED** (high impact)
2. chr1:161509955:A>G (p.His167Arg) - Missense in IgG binding domain

### Predictions Generated
- **RNA-seq expression** (bone marrow, liver, B cells, T cells)
- **CAGE** (transcription start sites)
- **Splice sites & usage**
- **Chromatin accessibility** (DNase, ATAC)

### Output Structure
```
variant_analysis_results/
├── GBP1/
│   ├── expression/ (baseline)
│   ├── variants/ (4 plots: 2 variants × 2 zoom levels)
│   └── chromatin/ (accessibility)
└── FCGR2B/
    ├── expression/ (baseline)
    ├── variants/ (4 plots: 2 variants × 2 zoom levels)
    └── chromatin/ (accessibility)
```

Each variant has **UNZOOMED** (1MB gene context) and **ZOOMED** (20KB variant-centered) views.

## Data Dependencies
- Uses `../data/` directory (shared with main analysis)
- Loads GENCODE v46 annotations from AlphaGenome cloud storage
- Requires ALPHA_GENOME_KEY environment variable

## Key Features
- REF/ALT overlay visualization (grey vs red)
- Variant annotation markers (orange vertical lines)
- Strand-specific filtering for negative strand genes
- High-resolution outputs (300 DPI)

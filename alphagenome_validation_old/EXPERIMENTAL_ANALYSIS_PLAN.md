# Experimental Variant AlphaGenome Analysis Plan

## Overview
Comprehensive multi-scale functional validation of 15 variants from our Nextflow trio analysis (HG002, HG003, HG004) using Google DeepMind's AlphaGenome AI.

## Selected Variants

### ITPKB (Inositol-trisphosphate 3-kinase B) - 5 variants
- chr1:228,742,035 A→G (intronic, heterozygous) - HG002
- chr1:228,738,074 A→G (intronic, heterozygous) - HG003  
- chr1:228,742,035 A→G (intronic, heterozygous) - HG003
- chr1:228,645,553 G→C (intergenic, homozygous_alt) - HG002
- chr1:228,645,553 G→C (intergenic, homozygous_alt) - HG004

### FCGR2B (Fc fragment of IgG receptor IIb) - 5 variants
- chr1:161,506,414 C→T (5'UTR/stop_gained, homozygous_alt) - HG002 ⭐
- chr1:161,506,414 C→T (5'UTR/stop_gained, heterozygous) - HG004 ⭐
- chr1:161,506,414 C→T (5'UTR/stop_gained, heterozygous) - HG003 ⭐
- chr1:161,509,955 A→G (missense p.His167Arg, heterozygous) - HG003 ⭐
- chr1:161,506,415 A→G (missense/5'UTR p.Gln63Arg, homozygous_alt) - HG002 ⭐

### GBP1 (Guanylate binding protein 1) - 5 variants
- chr1:89,014,075 T→C (splice_region c.625+8A>G, homozygous_alt) - HG002 ⭐
- chr1:89,014,075 T→C (splice_region c.625+8A>G, homozygous_alt) - HG003 ⭐
- chr1:89,014,075 T→C (splice_region c.625+8A>G, heterozygous) - HG004 ⭐
- chr1:89,013,380 G→A (missense p.Arg225Trp, heterozygous) - HG002 ⭐
- chr1:89,013,391 C→T (missense p.Arg221Gln, heterozygous) - HG002 ⭐

⭐ = High-impact variant (missense, splice, stop-gained)

## Analysis Parameters (API-Validated)

### Window Sizes
- **16 KB** (16,384 bp) - Fine detail, immediate variant context
- **512 KB** (524,288 bp) - Local gene regulatory landscape  
- **1 MB** (1,048,576 bp) - Broad regional chromatin context

### Prediction Types (6 modalities)
1. **RNA_SEQ** - Gene expression levels
2. **SPLICE_JUNCTIONS** - Splice junction usage patterns
3. **ATAC** - Chromatin accessibility
4. **CAGE** - Transcription start site usage
5. **CHIP_HISTONE** - Histone modifications (H3K27AC, H3K4ME3, H3K36ME3, etc.)
6. **CHIP_TF** - Transcription factor binding sites

### Analysis Scope
- **Total variants**: 15
- **Scales per variant**: 3 (16KB, 512KB, 1MB)
- **Prediction types**: 6
- **Total API calls**: 15 × 3 × 6 = **270 predictions**
- **Estimated runtime**: ~2.5 hours (0.5s per call + API latency)

## Why These Variants?

### Selection Criteria
1. **Functional impact** - Missense, splice region, 5'UTR regulatory
2. **Quality** - PASS filter in VCF
3. **Zygosity** - Heterozygous or homozygous alternate (real variants)
4. **Gene relevance** - ITPKB, FCGR2B, GBP1 are well-characterized genes
5. **Variant diversity** - Mix of coding, regulatory, and splice variants

### Expected Insights
- **Gene expression changes** - How variants affect transcription
- **Splicing alterations** - Impact on splice site usage
- **Regulatory effects** - Changes in chromatin state and TF binding
- **Multi-scale context** - How effects propagate across genomic scales

## Output Structure
```
experimental_alphagenome_results/
├── predictions/          # Raw .npz prediction arrays
│   ├── ITPKB_HG002_chr1_pos228742035_16KB_RNA_SEQ.npz
│   ├── ITPKB_HG002_chr1_pos228742035_512KB_RNA_SEQ.npz
│   └── ... (270 total files)
├── statistics/           # Per-variant JSON results
│   ├── ITPKB_HG002_chr1_pos228742035_result.json
│   └── ... (15 files)
└── analysis_summary.json # Complete analysis metadata
```

## Next Steps After Analysis
1. **Visualization** - Generate 15 plots per variant (heatmaps, tracks, comparisons)
2. **Statistical analysis** - Quantify effect sizes, compute integrated scores
3. **Cross-variant comparison** - Identify patterns across genes
4. **Report generation** - HTML interactive reports with findings

## Running the Analysis

### Prerequisites
✅ API key validated (ALPHA_GENOME_KEY)  
✅ Interval sizes verified (16KB, 512KB, 1MB)  
✅ Output types confirmed (6 available)  
✅ 15 variants selected from trio VCF data  

### Execution
```bash
cd /mnt/work_1/gest9386/CU_Boulder/MCDB-4520/MCDB-5520-group-project-analysis/alphagenome_validation/scripts
./run_experimental.sh
```

### Monitoring
```bash
# Watch progress
tail -f experimental_analysis.log

# Check if running
ps aux | grep run_experimental_analysis
```

## Validation Approach

This analysis validates **our experimental Nextflow pipeline results** by:
1. Confirming variant functional impacts with AI predictions
2. Identifying potential false positives (benign variants)
3. Discovering high-impact variants for follow-up
4. Understanding multi-scale regulatory context

This is much more appropriate than analyzing arbitrary HGVS coordinates from literature, as we're validating **our own experimental data** from the trio analysis we performed.

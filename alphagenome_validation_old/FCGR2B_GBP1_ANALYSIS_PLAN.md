# FCGR2B & GBP1 Deep Dive Analysis

**Date:** November 26, 2025  
**Focus:** Immune genes REF vs ALT comparison across all AlphaGenome tracks

---

## Objective

Compare reference (REF) and alternate (ALT) alleles for all variants in **FCGR2B** and **GBP1** genes across all available AlphaGenome prediction tracks to identify which biological processes are most affected by these variants.

---

## Gene Background

### FCGR2B (Fc Fragment of IgG Receptor IIb)
- **Location:** chr1:161,505,000-161,521,000
- **Function:** Low-affinity IgG receptor, immune system inhibitor
- **Disease relevance:** Autoimmune disorders, immune response regulation
- **Variants:** 5 (3× 5'UTR, 1× missense, 1× 5'UTR)

### GBP1 (Guanylate Binding Protein 1)
- **Location:** chr1:89,013,000-89,015,000
- **Function:** Interferon-induced GTPase, immune defense
- **Disease relevance:** Viral/bacterial resistance, inflammatory responses
- **Variants:** 5 (3× splice region, 2× missense)

---

## Variants Summary

### FCGR2B Variants (n=5)

| Position | REF→ALT | Type | Effect | Sample | Genotype |
|----------|---------|------|--------|--------|----------|
| 161,506,414 | C→T | 5'UTR | c.-347C>T | HG002 | homozygous_alt |
| 161,506,414 | C→T | 5'UTR | c.-347C>T | HG004 | heterozygous |
| 161,506,414 | C→T | 5'UTR | c.-347C>T | HG003 | heterozygous |
| 161,509,955 | A→G | Missense | p.His28Arg | HG003 | heterozygous |
| 161,506,415 | A→G | 5'UTR | c.-346A>G | HG002 | homozygous_alt |

**Key observation:** Position 161,506,414 appears in all three samples (trio), position 161,506,415 is adjacent (+1 bp).

### GBP1 Variants (n=5)

| Position | REF→ALT | Type | Effect | Sample | Genotype |
|----------|---------|------|--------|--------|----------|
| 89,014,075 | T→C | Splice region | c.625+8A>G | HG002 | homozygous_alt |
| 89,014,075 | T→C | Splice region | c.625+8A>G | HG003 | homozygous_alt |
| 89,014,075 | T→C | Splice region | c.625+8A>G | HG004 | heterozygous |
| 89,013,380 | G→A | Missense | p.Arg225Trp | HG002 | heterozygous |
| 89,013,391 | C→T | Missense | p.Arg221Gln | HG002 | heterozygous |

**Key observation:** Position 89,014,075 appears in all three samples (trio), two additional missense variants at positions 11bp apart.

---

## AlphaGenome Analysis Plan

### Tracks to Analyze (6 total)

1. **RNA_SEQ** - Gene expression levels
   - *Process:* Transcription output
   - *Metric:* Expression fold-change REF vs ALT

2. **SPLICE_JUNCTIONS** - Splice junction usage
   - *Process:* RNA splicing patterns
   - *Metric:* Junction usage differences

3. **ATAC** - Chromatin accessibility (ATAC-seq)
   - *Process:* DNA accessibility
   - *Metric:* Accessibility score changes

4. **CAGE** - Transcription start sites
   - *Process:* Transcription initiation
   - *Metric:* TSS usage differences

5. **CHIP_HISTONE** - Histone modifications
   - *Process:* Epigenetic regulation
   - *Metric:* H3K27AC, H3K4ME3, H3K36ME3 changes

6. **CHIP_TF** - Transcription factor binding
   - *Process:* TF recruitment
   - *Metric:* Binding site affinity changes

### Analysis Strategy

For each of 10 variants:
1. Generate REF prediction (1 Mbp window)
2. Generate ALT prediction (same window with variant applied)
3. Calculate delta (|ALT - REF|)
4. Quantify:
   - Maximum delta (peak effect)
   - Mean delta (average effect)
   - Spatial distribution of effects

### Window Size
- **1 Mbp (1,048,576 bp)** - Captures:
  - Gene body
  - Regulatory elements (enhancers, silencers)
  - Long-range chromatin interactions
  - Adjacent gene effects

---

## Expected Outcomes

### Hypothesis 1: FCGR2B 5'UTR variants affect transcription initiation
- **Tracks most affected:** CAGE, CHIP_HISTONE (H3K4ME3)
- **Process:** Transcription start site usage, promoter activity
- **Prediction:** ALT alleles reduce TSS usage → lower expression

### Hypothesis 2: FCGR2B missense (p.His28Arg) affects gene expression
- **Tracks most affected:** RNA_SEQ, possibly SPLICE_JUNCTIONS
- **Process:** mRNA stability, nonsense-mediated decay
- **Prediction:** Altered protein sequence → mRNA instability

### Hypothesis 3: GBP1 splice region variants disrupt splicing
- **Tracks most affected:** SPLICE_JUNCTIONS (primary)
- **Process:** Exon inclusion/exclusion
- **Prediction:** ALT alleles cause aberrant splicing → isoform shifts

### Hypothesis 4: GBP1 missense variants affect chromatin state
- **Tracks most affected:** CHIP_HISTONE, ATAC
- **Process:** Secondary chromatin effects from protein dysfunction
- **Prediction:** Missense → protein misfolding → altered feedback regulation

---

## Deliverables

1. **Jupyter Notebook** - `FCGR2B_GBP1_RefAlt_Analysis.ipynb`
   - Interactive analysis with visualizations
   - Step-by-step REF vs ALT comparisons

2. **Summary Tables**
   - `FCGR2B_GBP1_ref_alt_summary.csv` - All delta metrics
   - `FCGR2B_GBP1_detailed_results.json` - Full prediction arrays

3. **Visualizations**
   - Heatmap: Variant × Track impact matrix
   - Gene comparison: FCGR2B vs GBP1 across tracks
   - Process ranking: Which biological processes most affected

4. **Key Findings Report**
   - Top 10 highest impact changes
   - Process-level impact summary
   - Biological interpretation

---

## Computational Requirements

- **API calls:** 120 (10 variants × 6 tracks × 2 alleles)
- **Estimated time:** 30-60 minutes (depending on API response times)
- **Storage:** ~500 MB for prediction arrays
- **Memory:** ~8 GB RAM for analysis

---

## Next Steps

1. Run notebook to generate all REF vs ALT predictions
2. Analyze delta metrics to identify highest impact tracks
3. Focus detailed visualization on top 3 most affected processes
4. Generate publication-quality figures
5. Write biological interpretation for each gene

---

## Status

- [ ] Notebook created
- [ ] AlphaGenome client initialized
- [ ] Predictions generated (0/120)
- [ ] Summary statistics calculated
- [ ] Visualizations created
- [ ] Key findings documented
- [ ] Results exported

---

## Notes

- This analysis complements previous work on de novo variants
- Focus is specifically on immune function genes
- REF vs ALT approach isolates variant-specific effects
- Results will inform experimental validation priorities

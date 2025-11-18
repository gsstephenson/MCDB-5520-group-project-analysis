# MCDB-5520 Trio Variant Analysis: De Novo Variants + AI Functional Validation

**Analysis by:** George Stephenson  
**Date:** November 2025

---

## Overview

Analyzed chromosome 1 variants from GIAB trio (HG002/child, HG003/father, HG004/mother) to identify de novo variants and validate their functional impact using AlphaGenome AI. Found 70 de novo candidates; AI analysis revealed 5 ClinVar-benign variants that are likely pathogenic.

**Critical Context:** The GIAB Ashkenazi Jewish trio are **healthy reference individuals** used for sequencing validation - they exhibit no clinical disease despite carrying potentially pathogenic variants. This provides a unique opportunity to study dormant or incompletely penetrant disease risk variants.

### Key Results

- **70 de novo candidates** (chr1) - present in child, absent from both parents
- **55/70 (79%)** validated with AlphaGenome 
- **6 likely pathogenic** variants identified in healthy individuals
- **5 ClinVar misclassifications** detected (benign variants with strong functional impact)

**Major Finding:** AlphaGenome identified variants in serious disease genes (TARDBP/ALS, CASQ2/cardiac death, LAMC1-2/cancer, TGFB2/tumor suppressor) marked "benign" in ClinVar but showing significant splicing disruption (8-22%). These individuals are healthy, suggesting these variants represent **dormant disease risk alleles** with incomplete penetrance, requiring additional genetic/environmental factors for disease manifestation.

---

## Analysis Summary for Write-Up

### 1. De Novo Variant Discovery (BCFtools)

**Method:** BCFtools `isec -n=1` to find variants in child only  
**Dataset:** ClinVar-annotated chr1 variants (~3,600 per sample from ~530,000 total)

**Results:**
- **70 de novo candidates** in HG002 (child)
- **3 pathogenic/conflicting**: DBT, SPTA1, FH (all non-coding)
- **476** father-specific, **537** mother-specific variants

**Important:** These are candidates (not confirmed de novo) - could be false negatives in parents due to low coverage.

**Files:** `bcf_analysis/0000.vcf`, `analysis_results/de_novo_candidate_variants.csv`

---

### 2. ClinVar Annotations (3 Pathogenic/Conflicting)

All 3 are non-coding with conflicting classifications:

| Gene | Position | Type | Disease | ClinVar Status |
|------|----------|------|---------|----------------|
| **DBT** | 100189351 | 3'UTR del | Maple syrup urine disease | Conflicting (2 path/3 uncertain) |
| **SPTA1** | 158668075 | Intronic | Hereditary spherocytosis | Conflicting (3 uncertain/1 benign) |
| **FH** | 241500602 | Intronic | HLRCC cancer syndrome | Conflicting (majority benign) |

**Challenge:** Non-coding + conflicting = need functional validation

**Files:** `analysis_results/variant_summary_statistics.csv`, `bcf_analysis/INTERPRETATION.md`

---

### 3. AlphaGenome AI Validation (Google DeepMind)

**Method:** Multimodal predictions (RNA-seq, splicing, chromatin) on 1 Mbp windows  
**Success:** 55/70 variants (79%)  
**Thresholds:** Splice Δ >0.08 OR chromatin >0.5 = likely pathogenic

**Results:**
- 6 likely pathogenic (11%)
- 11 likely benign (20%)
- 38 uncertain (69%)

**Critical Discovery - ClinVar Discrepancies:**

5 "Benign" variants show strong functional impact:

| Gene | Disease | Splice Disruption | Clinical Risk |
|------|---------|-------------------|---------------|
| **TARDBP** | ALS/FTD | 12.7% | Neurodegeneration |
| **CASQ2** | Cardiac arrhythmia | 13.4% | Sudden death |
| **LAMC1** | Cancer metastasis | 14.4% | Tumor progression |
| **LAMC2** | Cancer metastasis | 11.5% | Tumor progression |
| **TGFB2** | Tumor suppressor | **22.1%** (highest) | Cancer risk |

Plus: **FH** (conflicting→likely pathogenic) - HLRCC cancer syndrome, 11.5% splice impact

**Interpretation:** All 6 are intronic/splice-region variants showing splicing disruption missed by standard tools. ClinVar likely lacked functional RNA data.

**Files:** `alphagenome_validation/outputs/validation_summary.csv` (summary table), `*_results.json` (55 detailed files)

---

## Project Structure

```
MCDB-5520-group-project-analysis/
│
├── README.md                              # This file
│
├── bcf_analysis/                          # BCFtools de novo identification
│   ├── 0000.vcf                           # 70 de novo variants ★
│   ├── COMPREHENSIVE_ANALYSIS_SUMMARY.md  # Methodology ★
│   └── INTERPRETATION.md                  # Clinical interpretation ★
│
├── alphagenome_validation/outputs/        # AlphaGenome results
│   ├── validation_summary.csv             # Main results table ★★
│   └── *_results.json                     # 55 detailed predictions
│
├── analysis_results/                      # Summary statistics
│   ├── de_novo_candidate_variants.csv     # All 70 variants annotated ★
│   └── variant_summary_statistics.csv     # Trio counts ★
│
├── analysis.ipynb                         # Exploratory analysis notebook
│
└── data/                                  # Raw VCF files (not in repo)

★★ = Critical for write-up
```

---

## Methods Summary

### BCFtools De Novo Detection
```bash
bcftools isec -n=1 HG002_clinvar.vcf HG003_clinvar.vcf HG004_clinvar.vcf
```
- Identifies variants in child only (not in either parent)
- Used ClinVar subset (~3,600/530,000 variants) for clinical relevance
- Chr1 only (computationally tractable, gene-rich)
- **Caveat:** Candidates only - could be false negatives in parents

### AlphaGenome Validation
- **Predictions:** RNA-seq (expression), splice junctions, ATAC-seq (chromatin)
- **Window:** 1 Mbp centered on variant
- **Classification:** Splice Δ >0.08 OR chromatin >0.5 = likely pathogenic
- **Script:** `validate_variants.py`

---

---

## Key Findings for Discussion

### 1. ClinVar Misclassifications (5 variants)
AI revealed "benign" variants with strong functional impact:
- **TARDBP** (ALS gene): 12.7% splicing disruption
- **CASQ2** (cardiac death): 13.4% splicing disruption  
- **LAMC1/LAMC2** (cancer): 11-14% splicing disruption
- **TGFB2** (tumor suppressor): 22.1% splicing disruption (highest)

**Why missed:** ClinVar lacks RNA functional data; traditional tools classify intronic variants as "modifier/low impact"

### 2. Non-Coding Variants Matter
All 6 likely pathogenic variants are intronic/splice-region:
- Primary mechanism: Splicing disruption (8-22% change)
- Secondary: Regulatory element changes
- Challenge: Standard tools underestimate non-coding impact

### 3. Resolving Conflicting Classifications
- **FH** (cancer gene): ClinVar majority "benign" → AlphaGenome "likely pathogenic" (11.5% splice disruption)
- **DBT/SPTA1**: AlphaGenome supports benign/uncertain (minimal functional impact)

---

## For Write-Up

### Methods
- BCFtools: `isec -n=1` on ClinVar-annotated chr1 variants to find de novo candidates
- AlphaGenome: Multimodal AI predictions (RNA/splicing/chromatin), 1 Mbp windows
- Classification thresholds: Splice Δ >0.08 OR chromatin >0.5 = likely pathogenic

### Results
- 70 de novo candidates (3 with pathogenic/conflicting ClinVar annotations)
- 55/70 successfully validated: 6 likely pathogenic, 11 benign, 38 uncertain
- 5 ClinVar-benign variants show strong functional impact (potential misclassifications)

### Discussion Points
- **ClinVar limitations:** Database lacks RNA functional data for non-coding variants
- **AI value:** AlphaGenome captures splicing/regulatory effects missed by standard annotation tools
- **Non-coding impact:** All 6 likely pathogenic variants are intronic/splice-region (not protein-changing)
- **Clinical concern:** Variants in CASQ2 (cardiac death risk), TARDBP (ALS), FH (cancer syndrome)
- **Healthy phenotype paradox:** Study subjects are healthy despite pathogenic predictions → incomplete penetrance

### Files to Reference
- `validation_summary.csv` - Main results table
- `de_novo_candidate_variants.csv` - All 70 variants
- `COMPREHENSIVE_ANALYSIS_SUMMARY.md` - Methods details
- Individual JSON files - Detailed evidence

---

## Reproducibility

```bash
# 1. De novo identification
cd bcf_analysis && ./run_bcftools_denovo_analysis.sh

# 2. AlphaGenome validation  
export ALPHA_GENOME_KEY="your_key"
cd alphagenome_validation/scripts && python3 validate_variants.py

# 3. Explore (optional)
jupyter notebook analysis.ipynb
```

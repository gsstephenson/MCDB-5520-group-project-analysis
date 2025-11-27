#!/usr/bin/env python3
"""
Experimental Variant Definitions from Trio Nextflow Analysis

Selected variants from ITPKB, FCGR2B, and GBP1 genes for comprehensive
AlphaGenome validation. These variants are from our experimental trio
analysis (HG002, HG003, HG004) chr1 data.

Generated: select_experimental_variants
Total variants: 15
"""

# No imports needed - this is just data definitions

# AlphaGenome prediction types (validated with API)
# Note: CHIP_HISTONE includes H3K27AC, H3K4ME3, H3K36ME3, etc.
ALL_OUTPUTS = [
    'RNA_SEQ',           # Gene expression
    'SPLICE_JUNCTIONS',  # Splice junction usage
    'ATAC',              # Chromatin accessibility
    'CAGE',              # Transcription start sites
    'CHIP_HISTONE',      # Histone modifications (includes H3K27AC, H3K4ME3, H3K36ME3)
    'CHIP_TF'            # Transcription factor binding
]

# Multi-scale analysis windows (validated with AlphaGenome API)
INTERVAL_SIZES = [16384, 524288, 1048576]  # 16KB, 512KB, 1MB (API-supported sizes)

# Experimental variants from trio analysis
EXPERIMENTAL_VARIANTS = [
    {
        'name': 'ITPKB_HG002_chr1_pos228742035',
        'gene': 'ITPKB',
        'sample': 'HG002_chr1',
        'chromosome': 'chr1',
        'position': 228742035,
        'ref': 'A',
        'alt': 'G',
        'type': 'intronic_variant',
        'genotype': 'heterozygous',
        'effect': 'intron_variant',
        'hgvs_c': 'c.322-1250A>G',
        'hgvs_p': '',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'ITPKB_HG003_chr1_pos228738074',
        'gene': 'ITPKB',
        'sample': 'HG003_chr1',
        'chromosome': 'chr1',
        'position': 228738074,
        'ref': 'A',
        'alt': 'G',
        'type': 'intronic_variant',
        'genotype': 'heterozygous',
        'effect': 'intron_variant',
        'hgvs_c': 'c.321+343A>G',
        'hgvs_p': '',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'ITPKB_HG003_chr1_pos228742035',
        'gene': 'ITPKB',
        'sample': 'HG003_chr1',
        'chromosome': 'chr1',
        'position': 228742035,
        'ref': 'A',
        'alt': 'G',
        'type': 'intronic_variant',
        'genotype': 'heterozygous',
        'effect': 'intron_variant',
        'hgvs_c': 'c.322-1250A>G',
        'hgvs_p': '',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'ITPKB_HG002_chr1_pos228645553',
        'gene': 'ITPKB',
        'sample': 'HG002_chr1',
        'chromosome': 'chr1',
        'position': 228645553,
        'ref': 'G',
        'alt': 'C',
        'type': 'other',
        'genotype': 'homozygous_alt',
        'effect': 'upstream_gene_variant',
        'hgvs_c': 'n.-3867C>G',
        'hgvs_p': '',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'ITPKB_HG004_chr1_pos228645553',
        'gene': 'ITPKB',
        'sample': 'HG004_chr1',
        'chromosome': 'chr1',
        'position': 228645553,
        'ref': 'G',
        'alt': 'C',
        'type': 'other',
        'genotype': 'homozygous_alt',
        'effect': 'upstream_gene_variant',
        'hgvs_c': 'n.-3867C>G',
        'hgvs_p': '',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'FCGR2B_HG002_chr1_pos161506414',
        'gene': 'FCGR2B',
        'sample': 'HG002_chr1',
        'chromosome': 'chr1',
        'position': 161506414,
        'ref': 'C',
        'alt': 'T',
        'type': '5_prime_UTR_variant',
        'genotype': 'homozygous_alt',
        'effect': '5_prime_UTR_variant',
        'hgvs_c': 'c.-347C>T',
        'hgvs_p': '',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'FCGR2B_HG004_chr1_pos161506414',
        'gene': 'FCGR2B',
        'sample': 'HG004_chr1',
        'chromosome': 'chr1',
        'position': 161506414,
        'ref': 'C',
        'alt': 'T',
        'type': '5_prime_UTR_variant',
        'genotype': 'heterozygous',
        'effect': '5_prime_UTR_variant',
        'hgvs_c': 'c.-347C>T',
        'hgvs_p': '',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'FCGR2B_HG003_chr1_pos161506414',
        'gene': 'FCGR2B',
        'sample': 'HG003_chr1',
        'chromosome': 'chr1',
        'position': 161506414,
        'ref': 'C',
        'alt': 'T',
        'type': '5_prime_UTR_variant',
        'genotype': 'heterozygous',
        'effect': '5_prime_UTR_variant',
        'hgvs_c': 'c.-347C>T',
        'hgvs_p': '',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'FCGR2B_HG003_chr1_pos161509955',
        'gene': 'FCGR2B',
        'sample': 'HG003_chr1',
        'chromosome': 'chr1',
        'position': 161509955,
        'ref': 'A',
        'alt': 'G',
        'type': 'missense_variant',
        'genotype': 'heterozygous',
        'effect': 'missense_variant',
        'hgvs_c': 'c.83A>G',
        'hgvs_p': 'p.His28Arg',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'FCGR2B_HG002_chr1_pos161506415',
        'gene': 'FCGR2B',
        'sample': 'HG002_chr1',
        'chromosome': 'chr1',
        'position': 161506415,
        'ref': 'A',
        'alt': 'G',
        'type': '5_prime_UTR_variant',
        'genotype': 'homozygous_alt',
        'effect': '5_prime_UTR_variant',
        'hgvs_c': 'c.-346A>G',
        'hgvs_p': '',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'GBP1_HG002_chr1_pos89014075',
        'gene': 'GBP1',
        'sample': 'HG002_chr1',
        'chromosome': 'chr1',
        'position': 89014075,
        'ref': 'T',
        'alt': 'C',
        'type': 'splice_variant',
        'genotype': 'homozygous_alt',
        'effect': 'splice_region_variant&intron_variant',
        'hgvs_c': 'c.625+8A>G',
        'hgvs_p': '',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'GBP1_HG003_chr1_pos89014075',
        'gene': 'GBP1',
        'sample': 'HG003_chr1',
        'chromosome': 'chr1',
        'position': 89014075,
        'ref': 'T',
        'alt': 'C',
        'type': 'splice_variant',
        'genotype': 'homozygous_alt',
        'effect': 'splice_region_variant&intron_variant',
        'hgvs_c': 'c.625+8A>G',
        'hgvs_p': '',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'GBP1_HG004_chr1_pos89014075',
        'gene': 'GBP1',
        'sample': 'HG004_chr1',
        'chromosome': 'chr1',
        'position': 89014075,
        'ref': 'T',
        'alt': 'C',
        'type': 'splice_variant',
        'genotype': 'heterozygous',
        'effect': 'splice_region_variant&intron_variant',
        'hgvs_c': 'c.625+8A>G',
        'hgvs_p': '',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'GBP1_HG002_chr1_pos89013380',
        'gene': 'GBP1',
        'sample': 'HG002_chr1',
        'chromosome': 'chr1',
        'position': 89013380,
        'ref': 'G',
        'alt': 'A',
        'type': 'missense_variant',
        'genotype': 'heterozygous',
        'effect': 'missense_variant',
        'hgvs_c': 'c.673C>T',
        'hgvs_p': 'p.Arg225Trp',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    },
    {
        'name': 'GBP1_HG002_chr1_pos89013391',
        'gene': 'GBP1',
        'sample': 'HG002_chr1',
        'chromosome': 'chr1',
        'position': 89013391,
        'ref': 'C',
        'alt': 'T',
        'type': 'missense_variant',
        'genotype': 'heterozygous',
        'effect': 'missense_variant',
        'hgvs_c': 'c.662G>A',
        'hgvs_p': 'p.Arg221Gln',
        'quality_filter': 'PASS',
        'outputs': ALL_OUTPUTS,
        'interval_sizes': INTERVAL_SIZES,
    }
]

# Summary statistics
GENE_SUMMARY = [{'gene': 'ITPKB', 'total': 1186, 'interesting': 577, 'selected': 5}, {'gene': 'FCGR2B', 'total': 323, 'interesting': 323, 'selected': 5}, {'gene': 'GBP1', 'total': 185, 'interesting': 149, 'selected': 5}]

if __name__ == '__main__':
    print("Experimental Variants for AlphaGenome Analysis")
    print("=" * 80)
    print(f"Total variants: {len(EXPERIMENTAL_VARIANTS)}")
    print("\nBy gene:")
    for summary in GENE_SUMMARY:
        print(f"  {summary['gene']:10s}: {summary['selected']}/{summary['interesting']} selected "
              f"({summary['total']} total in region)")
    
    print("\nVariants:")
    for v in EXPERIMENTAL_VARIANTS:
        print(f"  {v['name']:40s} {v['chromosome']}:{v['position']:,} "
              f"{v['ref']}→{v['alt']} [{v['type']}]")

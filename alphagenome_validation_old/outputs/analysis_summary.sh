#!/bin/bash

echo "=========================================================================="
echo "ALPHAGENOME VALIDATION ANALYSIS: 70 DE NOVO VARIANTS"
echo "=========================================================================="
echo ""

# Count classifications
echo "1. CLASSIFICATION SUMMARY:"
echo "────────────────────────────────────────────────────────────────────────"
echo -n "Total variants analyzed: "
tail -n +2 validation_summary.csv | wc -l

echo ""
echo "AlphaGenome Classifications:"
tail -n +2 validation_summary.csv | cut -d',' -f6 | sort | uniq -c | sort -rn

echo ""
echo "ClinVar Classifications:"
tail -n +2 validation_summary.csv | cut -d',' -f5 | sort | uniq -c | sort -rn

echo ""
echo "2. KEY DISCREPANCIES (ClinVar Benign vs AlphaGenome Pathogenic):"
echo "────────────────────────────────────────────────────────────────────────"
grep "Benign,Likely_pathogenic" validation_summary.csv | cut -d',' -f1,2,5,6,8 | head -10

echo ""
echo "3. HIGH CONFIDENCE PATHOGENIC FINDINGS:"
echo "────────────────────────────────────────────────────────────────────────"
grep "Likely_pathogenic,Moderate" validation_summary.csv | cut -d',' -f1,2,8 | head -10

echo ""
echo "4. SPLICE JUNCTION IMPACT:"
echo "────────────────────────────────────────────────────────────────────────"
grep "splicing: likely_pathogenic" validation_summary.csv | cut -d',' -f1,2,4 | head -10

echo ""
echo "5. CONFIDENT BENIGN CALLS:"
echo "────────────────────────────────────────────────────────────────────────"
grep "Likely_benign,Moderate" validation_summary.csv | cut -d',' -f1,2,4,8 | head -10

echo ""
echo "=========================================================================="


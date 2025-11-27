#!/bin/bash
#
# Run AlphaGenome Analysis on Experimental Trio Variants
#
# This script runs comprehensive multi-scale AlphaGenome validation on
# 15 variants selected from our experimental Nextflow trio analysis.
#

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "=================================================================="
echo "AlphaGenome Analysis: Experimental Trio Variants"
echo "=================================================================="
echo ""
echo "Variants: 15 (from ITPKB, FCGR2B, GBP1 genes)"
echo "Scales: 3 (100kb, 1MB, 2MB)"
echo "Prediction types: 8 (all modalities)"
echo "Total API calls: 360"
echo "Estimated time: ~3 hours"
echo ""
echo "Output: experimental_alphagenome_results/"
echo ""
echo "=================================================================="
echo ""

# Check API key
if [ -z "$ALPHA_GENOME_KEY" ]; then
    echo "ERROR: ALPHA_GENOME_KEY environment variable not set"
    exit 1
fi

# Check conda environment
if ! conda env list | grep -q alphagenome-env; then
    echo "ERROR: alphagenome-env conda environment not found"
    exit 1
fi

# Run analysis in background with logging
echo "Starting analysis..."
echo "Log file: experimental_analysis.log"
echo ""

nohup conda run -n alphagenome-env python run_experimental_analysis.py \
    > experimental_analysis.log 2>&1 &

ANALYSIS_PID=$!
echo "Analysis running in background (PID: $ANALYSIS_PID)"
echo ""
echo "Monitor progress:"
echo "  tail -f experimental_analysis.log"
echo ""
echo "Check if still running:"
echo "  ps -p $ANALYSIS_PID"
echo ""

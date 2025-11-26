#!/usr/bin/env python3
"""
Run Comprehensive AlphaGenome Analysis on Experimental Variants

This script runs the multi-scale, multi-modal AlphaGenome validation
on variants selected from our experimental trio Nextflow analysis.

Usage:
    python run_experimental_analysis.py
    
Or with custom output directory:
    python run_experimental_analysis.py --output custom_output_dir
"""

import os
import sys
import json
import time
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any

# Add AlphaGenome to path
sys.path.insert(0, '/opt/alphagenome')

try:
    from alphagenome.data import genome
    from alphagenome.models import dna_client
except ImportError as e:
    print(f"ERROR: Failed to import AlphaGenome modules: {e}")
    sys.exit(1)

# Import experimental variants
from experimental_variants_data import EXPERIMENTAL_VARIANTS, ALL_OUTPUTS, GENE_SUMMARY


class ExperimentalVariantValidator:
    """Validate experimental variants with AlphaGenome multi-scale analysis."""
    
    def __init__(self, api_key: str, output_dir: Path):
        """Initialize validator with API key and output directory."""
        self.api_key = api_key
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create subdirectories
        (self.output_dir / 'predictions').mkdir(exist_ok=True)
        (self.output_dir / 'statistics').mkdir(exist_ok=True)
        
        print(f"Initializing AlphaGenome client...")
        try:
            self.model = dna_client.create(api_key)
            print("✓ AlphaGenome client created successfully")
        except Exception as e:
            print(f"✗ Failed to create AlphaGenome client: {e}")
            raise
    
    def validate_variant_multiscale(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validate a single variant at multiple scales with all prediction types.
        
        Args:
            variant: Variant dictionary with genomic coordinates and metadata
            
        Returns:
            Dictionary with results for all scales and prediction types
        """
        print(f"\n{'='*80}")
        print(f"Validating: {variant['name']}")
        print(f"  Gene: {variant['gene']}")
        print(f"  Position: {variant['chromosome']}:{variant['position']:,}")
        print(f"  Change: {variant['ref']} → {variant['alt']}")
        print(f"  Type: {variant['type']}")
        print(f"  Genotype: {variant['genotype']}")
        print(f"  HGVS: {variant['hgvs_c']}")
        print(f"{'='*80}")
        
        results = {
            'variant_info': variant,
            'scales': {},
            'timestamp': datetime.now().isoformat()
        }
        
        # Analyze at each scale
        for interval_size in variant['interval_sizes']:
            scale_name = f"{interval_size//1000}kb" if interval_size < 1000000 else f"{interval_size//1048576}MB"
            print(f"\nAnalyzing at {scale_name} scale (window: {interval_size:,} bp)...")
            
            try:
                scale_result = self._analyze_at_scale(variant, interval_size)
                results['scales'][scale_name] = scale_result
                print(f"  ✓ Completed {len(scale_result)} prediction types")
            except Exception as e:
                print(f"  ✗ Error at {scale_name} scale: {e}")
                results['scales'][scale_name] = {'error': str(e)}
        
        return results
    
    def _analyze_at_scale(self, variant: Dict[str, Any], interval_size: int) -> Dict[str, Any]:
        """Analyze variant at a specific scale with all prediction types."""
        scale_results = {}
        
        # Create interval
        interval = genome.Interval(
            chromosome=variant['chromosome'],
            start=variant['position'] - (interval_size // 2),
            end=variant['position'] + (interval_size // 2)
        )
        
        # Create variant object
        ag_variant = genome.Variant(
            chromosome=variant['chromosome'],
            position=variant['position'],
            reference_bases=variant['ref'],
            alternate_bases=variant['alt']
        )
        
        for output_type_name in variant['outputs']:
            print(f"    - {output_type_name}...", end=' ', flush=True)
            
            try:
                # Get OutputType enum
                output_type = getattr(dna_client.OutputType, output_type_name)
                
                # Get prediction using correct API
                prediction = self.model.predict_variant(
                    interval=interval,
                    variant=ag_variant,
                    requested_outputs=[output_type],
                    ontology_terms=[]
                )
                
                # Save and analyze prediction
                pred_result = self._save_predictions(
                    prediction, 
                    output_type_name, 
                    variant, 
                    interval_size
                )
                
                scale_results[output_type_name] = pred_result
                print("✓")
                
                # Rate limiting
                time.sleep(0.5)
                
            except Exception as e:
                print(f"✗ ({e})")
                scale_results[output_type] = {'error': str(e)}
        
        return scale_results
    
    def _save_predictions(self, prediction, output_type: str, variant: Dict[str, Any], 
                         interval_size: int) -> Dict[str, Any]:
        """Save raw predictions and calculate statistics."""
        scale_name = f"{interval_size//1000}kb" if interval_size < 1000000 else f"{interval_size//1048576}MB"
        
        # Get prediction data
        pred_obj = getattr(prediction, output_type.lower(), None)
        if pred_obj is None:
            return {'error': f'No {output_type} prediction available'}
        
        # Extract values
        if hasattr(pred_obj, 'values'):
            values = pred_obj.values
        elif hasattr(pred_obj, 'value'):
            values = pred_obj.value
        else:
            return {'error': 'Unknown prediction format'}
        
        # Save raw prediction
        filename = f"{variant['name']}_{scale_name}_{output_type}.npz"
        filepath = self.output_dir / 'predictions' / filename
        np.savez_compressed(filepath, values=values, variant=variant, output_type=output_type)
        
        # Calculate statistics
        stats = {
            'shape': values.shape if hasattr(values, 'shape') else (len(values),),
            'mean': float(np.mean(values)),
            'std': float(np.std(values)),
            'min': float(np.min(values)),
            'max': float(np.max(values)),
            'median': float(np.median(values)),
            'prediction_file': str(filepath.name)
        }
        
        return stats
    
    def run_analysis(self, variants: List[Dict[str, Any]]) -> None:
        """Run comprehensive analysis on all variants."""
        print("="*80)
        print("COMPREHENSIVE EXPERIMENTAL VARIANT ANALYSIS")
        print("="*80)
        print(f"\nOutput directory: {self.output_dir}")
        print(f"Variants to analyze: {len(variants)}")
        
        # Calculate total API calls
        total_calls = sum(
            len(v['outputs']) * len(v['interval_sizes']) 
            for v in variants
        )
        print(f"Total API calls: {total_calls}")
        print(f"Estimated time: {total_calls * 0.5 / 60:.1f} minutes")
        
        # Analyze each variant
        all_results = []
        start_time = time.time()
        
        for i, variant in enumerate(variants, 1):
            print(f"\n{'#'*80}")
            print(f"# VARIANT {i}/{len(variants)}")
            print(f"{'#'*80}")
            
            try:
                result = self.validate_variant_multiscale(variant)
                all_results.append(result)
                
                # Save individual result
                result_file = self.output_dir / 'statistics' / f"{variant['name']}_result.json"
                with open(result_file, 'w') as f:
                    json.dump(result, f, indent=2)
                
            except Exception as e:
                print(f"\n✗ Failed to analyze variant: {e}")
                all_results.append({
                    'variant_info': variant,
                    'error': str(e)
                })
        
        # Save summary
        elapsed = time.time() - start_time
        summary = {
            'variants_analyzed': len(variants),
            'total_api_calls': total_calls,
            'elapsed_time_seconds': elapsed,
            'timestamp': datetime.now().isoformat(),
            'gene_summary': GENE_SUMMARY,
            'results': all_results
        }
        
        summary_file = self.output_dir / 'analysis_summary.json'
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        print("\n" + "="*80)
        print("ANALYSIS COMPLETE")
        print("="*80)
        print(f"Total time: {elapsed/60:.1f} minutes")
        print(f"Summary saved: {summary_file}")
        print(f"Results directory: {self.output_dir}")


def main():
    """Main execution function."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Run AlphaGenome analysis on experimental variants')
    parser.add_argument('--output', default='experimental_alphagenome_results',
                       help='Output directory for results')
    args = parser.parse_args()
    
    # Get API key
    api_key = os.environ.get('ALPHA_GENOME_KEY')
    if not api_key:
        print("ERROR: ALPHA_GENOME_KEY environment variable not set")
        sys.exit(1)
    
    # Initialize validator
    output_dir = Path(args.output)
    validator = ExperimentalVariantValidator(api_key, output_dir)
    
    # Run analysis
    validator.run_analysis(EXPERIMENTAL_VARIANTS)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
AlphaGenome Validation Script for De Novo Pathogenic Variants

This script validates 3 pathogenic de novo variants using AlphaGenome's
multimodal genomic predictions.

Variants:
1. DBT (chr1:100189351) - 3'UTR deletion, Maple syrup urine disease
2. SPTA1 (chr1:158668075) - Intronic deletion, Hereditary spherocytosis
3. FH (chr1:241500602) - Intronic microsatellite, Fumarase deficiency/HLRCC

Author: MCDB-4520 Bioinformatics Project
Date: December 2024
"""

import os
import sys
import json
import argparse
import time
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Any

# Add AlphaGenome to path
sys.path.insert(0, '/opt/alphagenome')

try:
    from alphagenome.data import genome
    from alphagenome.models import dna_client
    import numpy as np
except ImportError as e:
    print(f"ERROR: Failed to import AlphaGenome modules: {e}")
    print("Please ensure /opt/alphagenome is accessible and PYTHONPATH is set.")
    sys.exit(1)

# Import all 70 variants from auto-generated data file
from all_variants_data import ALL_VARIANTS

# Use imported variants
VARIANTS = ALL_VARIANTS


class AlphaGenomeValidator:
    """Validates variants using AlphaGenome API."""
    
    def __init__(self, api_key: str, output_dir: Path):
        """
        Initialize validator.
        
        Args:
            api_key: AlphaGenome API key
            output_dir: Directory to save results
        """
        self.api_key = api_key
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"Initializing AlphaGenome client...")
        try:
            self.model = dna_client.create(api_key)
            print("✓ AlphaGenome client created successfully")
        except Exception as e:
            print(f"✗ Failed to create AlphaGenome client: {e}")
            raise
    
    def validate_variant(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validate a single variant using AlphaGenome.
        
        Args:
            variant: Dictionary containing variant information
            
        Returns:
            Dictionary with validation results
        """
        print(f"\n{'='*80}")
        print(f"Processing: {variant['name']}")
        print(f"Gene: {variant['gene']}")
        print(f"Position: {variant['chromosome']}:{variant['position']}")
        print(f"Change: {variant['ref']} → {variant['alt']}")
        print(f"Strategy: {variant['analysis_strategy']}")
        print(f"{'='*80}\n")
        
        # Define genomic interval
        interval_start = variant['position'] - (variant['interval_size'] // 2)
        interval_end = variant['position'] + (variant['interval_size'] // 2)
        
        interval = genome.Interval(
            chromosome=variant['chromosome'],
            start=interval_start,
            end=interval_end
        )
        
        print(f"Genomic interval: {variant['chromosome']}:{interval_start}-{interval_end}")
        print(f"Interval size: {interval_end - interval_start:,} bp")
        
        # Create variant object
        ag_variant = genome.Variant(
            chromosome=variant['chromosome'],
            position=variant['position'],
            reference_bases=variant['ref'],
            alternate_bases=variant['alt']
        )
        
        # Create serializable variant info (convert OutputType enums to strings)
        variant_info_serializable = variant.copy()
        variant_info_serializable['outputs'] = [ot.name for ot in variant['outputs']]
        
        results = {
            'variant_info': variant_info_serializable,
            'interval': {
                'chromosome': variant['chromosome'],
                'start': interval_start,
                'end': interval_end,
                'size': interval_end - interval_start
            },
            'predictions': {},
            'timestamp': datetime.now().isoformat(),
            'status': 'pending'
        }
        
        # Make API calls for each output type
        for output_type in variant['outputs']:
            output_name = output_type.name
            print(f"\nRequesting {output_name} predictions...")
            
            try:
                # Call AlphaGenome API with retry logic
                outputs = self._call_api_with_retry(
                    interval=interval,
                    variant=ag_variant,
                    output_type=output_type
                )
                
                # Process and store results
                prediction_data = self._process_predictions(
                    outputs, output_type, variant['position']
                )
                results['predictions'][output_name] = prediction_data
                
                print(f"✓ {output_name} predictions completed")
                
            except Exception as e:
                print(f"✗ Failed to get {output_name} predictions: {e}")
                results['predictions'][output_name] = {
                    'error': str(e),
                    'status': 'failed'
                }
        
        # Determine overall status
        if all(p.get('status') != 'failed' for p in results['predictions'].values()):
            results['status'] = 'success'
        else:
            results['status'] = 'partial_failure'
        
        # Calculate pathogenicity score
        results['pathogenicity_assessment'] = self._assess_pathogenicity(
            results['predictions'], variant
        )
        
        return results
    
    def _call_api_with_retry(
        self, 
        interval: genome.Interval,
        variant: genome.Variant,
        output_type: dna_client.OutputType,
        max_retries: int = 3
    ) -> Any:
        """
        Call AlphaGenome API with exponential backoff retry.
        
        Args:
            interval: Genomic interval
            variant: Variant to predict
            output_type: Type of prediction output
            max_retries: Maximum number of retry attempts
            
        Returns:
            AlphaGenome prediction outputs
        """
        for attempt in range(max_retries):
            try:
                outputs = self.model.predict_variant(
                    interval=interval,
                    variant=variant,
                    requested_outputs=[output_type],
                    ontology_terms=[]  # Empty list for general predictions
                )
                return outputs
                
            except Exception as e:
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    print(f"  Retry {attempt + 1}/{max_retries} after {wait_time}s... ({e})")
                    time.sleep(wait_time)
                else:
                    raise
    
    def _process_predictions(
        self,
        outputs: Any,
        output_type: dna_client.OutputType,
        variant_position: int
    ) -> Dict[str, Any]:
        """
        Process AlphaGenome prediction outputs.
        
        Args:
            outputs: Raw AlphaGenome outputs
            output_type: Type of output
            variant_position: Position of variant
            
        Returns:
            Processed prediction data
        """
        result = {
            'output_type': output_type.name,
            'status': 'success'
        }
        
        if output_type == dna_client.OutputType.RNA_SEQ:
            # Extract RNA-seq data
            ref_data = outputs.reference.rna_seq
            alt_data = outputs.alternate.rna_seq
            
            # Get values around variant position
            variant_idx = variant_position - ref_data.interval.start
            window = 50  # ±50bp window
            start_idx = max(0, variant_idx - window)
            end_idx = min(len(ref_data.values), variant_idx + window)
            
            ref_values = ref_data.values[start_idx:end_idx]
            alt_values = alt_data.values[start_idx:end_idx]
            
            # Calculate statistics
            ref_mean = float(np.mean(ref_values))
            alt_mean = float(np.mean(alt_values))
            fold_change = alt_mean / ref_mean if ref_mean != 0 else 0
            
            result.update({
                'reference_mean': ref_mean,
                'alternate_mean': alt_mean,
                'fold_change': fold_change,
                'log2_fold_change': float(np.log2(fold_change)) if fold_change > 0 else None,
                'absolute_difference': float(abs(alt_mean - ref_mean)),
                'effect_direction': 'increased' if alt_mean > ref_mean else 'decreased',
                'window_size': end_idx - start_idx
            })
            
        elif output_type == dna_client.OutputType.SPLICE_JUNCTIONS:
            # Extract splice junction data
            ref_junctions = outputs.reference.splice_junctions
            alt_junctions = outputs.alternate.splice_junctions
            
            # Get values arrays (2D: positions x junction types)
            ref_values = ref_junctions.values
            alt_values = alt_junctions.values
            
            # Calculate absolute differences
            diff = np.abs(alt_values - ref_values)
            
            # Global statistics
            max_change = float(np.max(diff))
            mean_change = float(np.mean(diff))
            
            # Try to find changes near the variant position
            # Map genomic position to array index
            interval_start = ref_junctions.interval.start
            interval_end = ref_junctions.interval.end
            array_length = ref_values.shape[0]
            
            # Calculate position in array
            # The array represents the interval with some resolution
            bp_per_position = (interval_end - interval_start) / array_length
            variant_array_idx = int((variant_position - interval_start) / bp_per_position)
            
            # Extract local region around variant (if in bounds)
            if 0 <= variant_array_idx < array_length:
                window = max(1, min(50, array_length // 20))  # At least 1, adaptive window
                start_idx = max(0, variant_array_idx - window)
                end_idx = min(array_length, variant_array_idx + window + 1)
                
                # Ensure we have at least some data
                if end_idx > start_idx and start_idx < array_length:
                    local_diff = diff[start_idx:end_idx, :]
                    if local_diff.size > 0:
                        local_max_change = float(np.max(local_diff))
                        local_mean_change = float(np.mean(local_diff))
                    else:
                        local_max_change = max_change
                        local_mean_change = mean_change
                else:
                    local_max_change = max_change
                    local_mean_change = mean_change
            else:
                local_max_change = max_change
                local_mean_change = mean_change
            
            result.update({
                'max_splice_change': max_change,
                'mean_splice_change': mean_change,
                'local_max_splice_change': local_max_change,
                'local_mean_splice_change': local_mean_change,
                'junction_array_shape': list(ref_values.shape),
                'variant_array_index': variant_array_idx if 0 <= variant_array_idx < array_length else None,
                'splicing_affected': local_max_change > 0.1
            })
            
        elif output_type == dna_client.OutputType.ATAC:
            # Extract ATAC (chromatin accessibility) data
            ref_atac = outputs.reference.atac
            alt_atac = outputs.alternate.atac
            
            # Get values around variant
            variant_idx = variant_position - ref_atac.interval.start
            window = 50
            start_idx = max(0, variant_idx - window)
            end_idx = min(len(ref_atac.values), variant_idx + window)
            
            ref_access = ref_atac.values[start_idx:end_idx]
            alt_access = alt_atac.values[start_idx:end_idx]
            
            mean_change = float(np.mean(np.abs(alt_access - ref_access)))
            max_change = float(np.max(np.abs(alt_access - ref_access)))
            
            result.update({
                'mean_accessibility_change': mean_change,
                'max_accessibility_change': max_change,
                'accessibility_affected': max_change > 0.2,
                'window_size': end_idx - start_idx
            })
        
        return result
    
    def _assess_pathogenicity(
        self,
        predictions: Dict[str, Any],
        variant: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Assess variant pathogenicity based on AlphaGenome predictions.
        
        Args:
            predictions: Dictionary of prediction results
            variant: Variant information
            
        Returns:
            Pathogenicity assessment
        """
        assessment = {
            'variant': variant['name'],
            'gene': variant['gene'],
            'evidence': []
        }
        
        # RNA-seq evidence
        if 'RNA_SEQ' in predictions and predictions['RNA_SEQ'].get('status') == 'success':
            rna = predictions['RNA_SEQ']
            fold_change = rna.get('fold_change', 1.0)
            
            if abs(np.log2(fold_change)) > 1.0:  # >2x change
                assessment['evidence'].append({
                    'type': 'gene_expression',
                    'effect': 'large',
                    'value': fold_change,
                    'interpretation': f"{rna['effect_direction']} expression by {fold_change:.2f}x - strong functional impact"
                })
            elif abs(np.log2(fold_change)) > 0.585:  # >1.5x change
                assessment['evidence'].append({
                    'type': 'gene_expression',
                    'effect': 'moderate',
                    'value': fold_change,
                    'interpretation': f"{rna['effect_direction']} expression by {fold_change:.2f}x - moderate impact"
                })
            else:
                assessment['evidence'].append({
                    'type': 'gene_expression',
                    'effect': 'minimal',
                    'value': fold_change,
                    'interpretation': f"Minimal expression change ({fold_change:.2f}x) - likely benign"
                })
        
        # Splicing evidence
        if 'SPLICE_JUNCTIONS' in predictions and predictions['SPLICE_JUNCTIONS'].get('status') == 'success':
            splice = predictions['SPLICE_JUNCTIONS']
            # Use local change (near variant) if available, otherwise global max
            local_change = splice.get('local_max_splice_change', splice.get('max_splice_change', 0.0))
            global_change = splice.get('max_splice_change', 0.0)
            
            if local_change > 0.1:
                assessment['evidence'].append({
                    'type': 'splicing',
                    'effect': 'likely_pathogenic',
                    'value': local_change,
                    'interpretation': f"Local splice junction change of {local_change:.3f} (global max: {global_change:.3f}) - likely affects splicing"
                })
            elif local_change > 0.05:
                assessment['evidence'].append({
                    'type': 'splicing',
                    'effect': 'uncertain',
                    'value': local_change,
                    'interpretation': f"Moderate splice change ({local_change:.3f}) near variant - uncertain impact"
                })
            else:
                assessment['evidence'].append({
                    'type': 'splicing',
                    'effect': 'minimal',
                    'value': local_change,
                    'interpretation': f"Minimal splice change ({local_change:.3f}) near variant - unlikely to affect splicing"
                })
        
        # Chromatin accessibility evidence
        if 'ATAC' in predictions and predictions['ATAC'].get('status') == 'success':
            atac = predictions['ATAC']
            max_change = atac.get('max_accessibility_change', 0.0)
            
            if max_change > 0.2:
                assessment['evidence'].append({
                    'type': 'chromatin_accessibility',
                    'effect': 'moderate',
                    'value': max_change,
                    'interpretation': f"Accessibility change of {max_change:.3f} - may affect regulatory function"
                })
            else:
                assessment['evidence'].append({
                    'type': 'chromatin_accessibility',
                    'effect': 'minimal',
                    'value': max_change,
                    'interpretation': f"Minimal accessibility change ({max_change:.3f}) - limited regulatory impact"
                })
        
        # Overall classification
        large_effects = sum(1 for e in assessment['evidence'] if e['effect'] in ['large', 'likely_pathogenic'])
        moderate_effects = sum(1 for e in assessment['evidence'] if e['effect'] in ['moderate', 'uncertain'])
        minimal_effects = sum(1 for e in assessment['evidence'] if e['effect'] == 'minimal')
        
        if large_effects >= 1:
            assessment['classification'] = 'Likely_pathogenic'
            assessment['confidence'] = 'High' if large_effects > 1 else 'Moderate'
        elif moderate_effects >= 1:
            assessment['classification'] = 'Uncertain_significance'
            assessment['confidence'] = 'Low'
        else:
            assessment['classification'] = 'Likely_benign'
            assessment['confidence'] = 'Moderate'
        
        # Compare with ClinVar
        assessment['clinvar_comparison'] = {
            'clinvar_classification': variant['clinvar_classification'],
            'clinvar_votes': variant.get('clnsigconf', 'N/A'),
            'agreement': 'To be manually assessed'
        }
        
        return assessment
    
    def save_results(self, results: Dict[str, Any], variant_name: str):
        """Save results to JSON file."""
        output_file = self.output_dir / f"{variant_name}_results.json"
        
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\n✓ Results saved to: {output_file}")
    
    def generate_summary(self, all_results: List[Dict[str, Any]]):
        """Generate summary CSV of all validations."""
        summary_file = self.output_dir / 'validation_summary.csv'
        
        with open(summary_file, 'w') as f:
            # Header
            f.write("Variant,Gene,Position,Type,ClinVar_Classification,")
            f.write("AlphaGenome_Classification,Confidence,Key_Evidence\n")
            
            # Data rows
            for result in all_results:
                info = result['variant_info']
                assessment = result.get('pathogenicity_assessment', {})
                
                evidence_summary = '; '.join([
                    f"{e['type']}: {e['effect']}" 
                    for e in assessment.get('evidence', [])
                ])
                
                f.write(f"{info['name']},")
                f.write(f"{info['gene']},")
                f.write(f"{info['chromosome']}:{info['position']},")
                f.write(f"{info['type']},")
                f.write(f"{info['clinvar_classification']},")
                f.write(f"{assessment.get('classification', 'N/A')},")
                f.write(f"{assessment.get('confidence', 'N/A')},")
                f.write(f"\"{evidence_summary}\"\n")
        
        print(f"\n✓ Summary saved to: {summary_file}")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Validate de novo variants using AlphaGenome'
    )
    parser.add_argument(
        '--api-key',
        type=str,
        default=os.environ.get('ALPHA_GENOME_KEY'),
        help='AlphaGenome API key (or set ALPHA_GENOME_KEY env var)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('../outputs'),
        help='Output directory for results'
    )
    parser.add_argument(
        '--variant',
        type=str,
        choices=['DBT', 'SPTA1', 'FH', 'all'],
        default='all',
        help='Which variant(s) to validate'
    )
    
    args = parser.parse_args()
    
    # Validate API key
    if not args.api_key:
        print("ERROR: AlphaGenome API key required!")
        print("Set ALPHAGENOME_API_KEY environment variable or use --api-key")
        sys.exit(1)
    
    print("="*80)
    print("AlphaGenome Variant Validation Pipeline")
    print("="*80)
    print(f"API Key: {args.api_key[:10]}...")
    print(f"Output Directory: {args.output_dir}")
    print(f"Timestamp: {datetime.now().isoformat()}")
    print("="*80)
    
    # Initialize validator
    validator = AlphaGenomeValidator(args.api_key, args.output_dir)
    
    # Select variants to process
    if args.variant == 'all':
        variants_to_process = VARIANTS
    else:
        variants_to_process = [v for v in VARIANTS if v['gene'] == args.variant]
    
    print(f"\nProcessing {len(variants_to_process)} variant(s)...")
    
    # Process each variant
    all_results = []
    for variant in variants_to_process:
        try:
            results = validator.validate_variant(variant)
            validator.save_results(results, variant['name'])
            all_results.append(results)
            
            # Print summary
            assessment = results.get('pathogenicity_assessment', {})
            print(f"\n{'─'*80}")
            print(f"SUMMARY for {variant['name']}:")
            print(f"  ClinVar: {variant['clinvar_classification']}")
            print(f"  AlphaGenome: {assessment.get('classification', 'N/A')}")
            print(f"  Confidence: {assessment.get('confidence', 'N/A')}")
            print(f"  Evidence:")
            for evidence in assessment.get('evidence', []):
                print(f"    - {evidence['interpretation']}")
            print(f"{'─'*80}\n")
            
        except Exception as e:
            print(f"\n✗ Failed to process {variant['name']}: {e}")
            import traceback
            traceback.print_exc()
    
    # Generate overall summary
    if all_results:
        validator.generate_summary(all_results)
    
    print("\n" + "="*80)
    print("Validation Complete!")
    print(f"Results saved to: {args.output_dir}")
    print("="*80)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
AlphaGenome Visualization Script

Generates publication-quality visualizations of AlphaGenome variant predictions.

Author: MCDB-4520 Bioinformatics Project
Date: December 2024
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Dict, Any, List

# Add AlphaGenome to path
sys.path.insert(0, '/opt/alphagenome')

try:
    from alphagenome.data import genome
    from alphagenome.models import dna_client
    from alphagenome.visualization import plot_components
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
except ImportError as e:
    print(f"ERROR: Failed to import required modules: {e}")
    print("Required: alphagenome, matplotlib, numpy")
    sys.exit(1)


class AlphaGenomeVisualizer:
    """Creates visualizations from AlphaGenome validation results."""
    
    def __init__(self, results_dir: Path, output_dir: Path):
        """
        Initialize visualizer.
        
        Args:
            results_dir: Directory containing JSON results
            output_dir: Directory to save visualizations
        """
        self.results_dir = Path(results_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set publication-quality matplotlib parameters
        plt.rcParams.update({
            'figure.figsize': (12, 8),
            'figure.dpi': 300,
            'font.size': 10,
            'axes.labelsize': 12,
            'axes.titlesize': 14,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 10,
            'font.family': 'sans-serif',
            'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans']
        })
    
    def load_results(self, variant_name: str) -> Dict[str, Any]:
        """Load results JSON file."""
        results_file = self.results_dir / f"{variant_name}_results.json"
        
        if not results_file.exists():
            raise FileNotFoundError(f"Results file not found: {results_file}")
        
        with open(results_file, 'r') as f:
            return json.load(f)
    
    def visualize_gene_expression(self, results: Dict[str, Any]) -> Path:
        """
        Visualize gene expression changes (RNA-seq).
        
        Args:
            results: Validation results dictionary
            
        Returns:
            Path to saved figure
        """
        variant_info = results['variant_info']
        predictions = results['predictions']
        
        if 'RNA_SEQ' not in predictions:
            print(f"No RNA-seq data for {variant_info['name']}")
            return None
        
        rna_data = predictions['RNA_SEQ']
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Top panel: Expression comparison
        ref_mean = rna_data['reference_mean']
        alt_mean = rna_data['alternate_mean']
        fold_change = rna_data['fold_change']
        
        categories = ['Reference', 'Alternate']
        values = [ref_mean, alt_mean]
        colors = ['#2E86AB', '#A23B72']
        
        bars = ax1.bar(categories, values, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
        ax1.set_ylabel('Mean Expression Level', fontweight='bold')
        ax1.set_title(
            f'{variant_info["gene"]} Gene Expression: {variant_info["ref"]}→{variant_info["alt"]}\n'
            f'Position: {variant_info["chromosome"]}:{variant_info["position"]:,}',
            fontweight='bold'
        )
        ax1.grid(axis='y', alpha=0.3, linestyle='--')
        
        # Add value labels on bars
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax1.text(
                bar.get_x() + bar.get_width()/2., height,
                f'{value:.3f}',
                ha='center', va='bottom', fontweight='bold'
            )
        
        # Add fold change annotation
        ax1.text(
            0.5, max(values) * 0.9,
            f'Fold Change: {fold_change:.2f}x\n'
            f'Log₂ FC: {rna_data.get("log2_fold_change", 0):.2f}\n'
            f'Direction: {rna_data["effect_direction"].upper()}',
            ha='center',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
            fontsize=11,
            fontweight='bold'
        )
        
        # Bottom panel: Pathogenicity assessment
        assessment = results['pathogenicity_assessment']
        classification = assessment.get('classification', 'Unknown')
        confidence = assessment.get('confidence', 'Unknown')
        
        ax2.axis('off')
        
        # Classification box
        class_color = {
            'Likely_pathogenic': '#E63946',
            'Uncertain_significance': '#F77F00',
            'Likely_benign': '#06A77D'
        }.get(classification, '#CCCCCC')
        
        summary_text = (
            f"ALPHAGENOME CLASSIFICATION\n"
            f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"
            f"Classification: {classification.replace('_', ' ').upper()}\n"
            f"Confidence: {confidence.upper()}\n\n"
            f"ClinVar Classification:\n{variant_info['clinvar_classification']}\n\n"
            f"Evidence Summary:\n"
        )
        
        for evidence in assessment.get('evidence', []):
            summary_text += f"  • {evidence['interpretation']}\n"
        
        ax2.text(
            0.5, 0.5,
            summary_text,
            ha='center', va='center',
            fontsize=11,
            family='monospace',
            bbox=dict(
                boxstyle='round,pad=1',
                facecolor=class_color,
                alpha=0.3,
                edgecolor='black',
                linewidth=2
            )
        )
        
        plt.tight_layout()
        
        # Save figure
        output_file = self.output_dir / f"{variant_info['name']}_gene_expression.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Gene expression visualization saved: {output_file}")
        return output_file
    
    def visualize_splicing(self, results: Dict[str, Any]) -> Path:
        """
        Visualize splicing predictions.
        
        Args:
            results: Validation results dictionary
            
        Returns:
            Path to saved figure
        """
        variant_info = results['variant_info']
        predictions = results['predictions']
        
        if 'SPLICE_JUNCTIONS' not in predictions:
            print(f"No splicing data for {variant_info['name']}")
            return None
        
        splice_data = predictions['SPLICE_JUNCTIONS']
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Top panel: Splice site probability changes
        # Note: AlphaGenome returns combined splice junction data, not separate donor/acceptor
        categories = ['Global Max', 'Local Max (near variant)']
        values = [
            splice_data.get('max_splice_change', 0),
            splice_data.get('local_max_splice_change', 0)
        ]
        
        colors = ['#E63946' if v > 0.1 else '#F77F00' if v > 0.05 else '#06A77D' for v in values]
        
        bars = ax1.bar(categories, values, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
        ax1.set_ylabel('Max Probability Change (Δp)', fontweight='bold')
        ax1.set_title(
            f'{variant_info["gene"]} Splicing Impact: {variant_info["ref"]}→{variant_info["alt"]}\n'
            f'Position: {variant_info["chromosome"]}:{variant_info["position"]:,} '
            f'({variant_info["type"].replace("_", " ")})',
            fontweight='bold'
        )
        ax1.axhline(y=0.1, color='red', linestyle='--', alpha=0.5, linewidth=2, label='Likely pathogenic (>0.1)')
        ax1.axhline(y=0.05, color='orange', linestyle='--', alpha=0.5, linewidth=2, label='Uncertain (>0.05)')
        ax1.legend(loc='upper right')
        ax1.grid(axis='y', alpha=0.3, linestyle='--')
        ax1.set_ylim(0, max(values) * 1.2 if max(values) > 0 else 0.2)
        
        # Add value labels
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax1.text(
                bar.get_x() + bar.get_width()/2., height,
                f'{value:.4f}',
                ha='center', va='bottom', fontweight='bold'
            )
        
        # Bottom panel: Assessment summary
        assessment = results['pathogenicity_assessment']
        classification = assessment.get('classification', 'Unknown')
        confidence = assessment.get('confidence', 'Unknown')
        
        ax2.axis('off')
        
        class_color = {
            'Likely_pathogenic': '#E63946',
            'Uncertain_significance': '#F77F00',
            'Likely_benign': '#06A77D'
        }.get(classification, '#CCCCCC')
        
        splicing_affected = splice_data.get('splicing_affected', False)
        max_change = splice_data.get('max_splice_change', 0)
        
        summary_text = (
            f"SPLICING ANALYSIS SUMMARY\n"
            f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"
            f"Splicing Affected: {'YES' if splicing_affected else 'NO'}\n"
            f"Max Change: {max_change:.4f}\n"
            f"Classification: {classification.replace('_', ' ').upper()}\n"
            f"Confidence: {confidence.upper()}\n\n"
            f"Interpretation:\n"
        )
        
        for evidence in assessment.get('evidence', []):
            if evidence['type'] == 'splicing':
                summary_text += f"  {evidence['interpretation']}\n"
        
        summary_text += f"\nClinVar: {variant_info['clinvar_classification']}"
        
        ax2.text(
            0.5, 0.5,
            summary_text,
            ha='center', va='center',
            fontsize=11,
            family='monospace',
            bbox=dict(
                boxstyle='round,pad=1',
                facecolor=class_color,
                alpha=0.3,
                edgecolor='black',
                linewidth=2
            )
        )
        
        plt.tight_layout()
        
        # Save figure
        output_file = self.output_dir / f"{variant_info['name']}_splicing_predictions.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Splicing visualization saved: {output_file}")
        return output_file
    
    def visualize_chromatin_accessibility(self, results: Dict[str, Any]) -> Path:
        """
        Visualize chromatin accessibility changes.
        
        Args:
            results: Validation results dictionary
            
        Returns:
            Path to saved figure
        """
        variant_info = results['variant_info']
        predictions = results['predictions']
        
        if 'ATAC' not in predictions:
            print(f"No ATAC data for {variant_info['name']}")
            return None
        
        atac_data = predictions['ATAC']
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Chromatin accessibility changes
        mean_change = atac_data.get('mean_accessibility_change', 0)
        max_change = atac_data.get('max_accessibility_change', 0)
        
        categories = ['Mean Change', 'Max Change']
        values = [mean_change, max_change]
        
        colors = ['#E63946' if max_change > 0.2 else '#F77F00' if max_change > 0.1 else '#06A77D']
        
        bars = ax.bar(categories, values, color=colors * 2, alpha=0.7, edgecolor='black', linewidth=1.5)
        ax.set_ylabel('Accessibility Change (Δ)', fontweight='bold')
        ax.set_title(
            f'{variant_info["gene"]} Chromatin Accessibility: {variant_info["ref"]}→{variant_info["alt"]}\n'
            f'Position: {variant_info["chromosome"]}:{variant_info["position"]:,}',
            fontweight='bold'
        )
        ax.axhline(y=0.2, color='red', linestyle='--', alpha=0.5, linewidth=2, label='High impact (>0.2)')
        ax.axhline(y=0.1, color='orange', linestyle='--', alpha=0.5, linewidth=2, label='Moderate impact (>0.1)')
        ax.legend(loc='upper right')
        ax.grid(axis='y', alpha=0.3, linestyle='--')
        
        # Add value labels
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width()/2., height,
                f'{value:.4f}',
                ha='center', va='bottom', fontweight='bold'
            )
        
        # Add interpretation box
        accessibility_affected = atac_data.get('accessibility_affected', False)
        
        interp_text = (
            f"Accessibility Impact: {'YES' if accessibility_affected else 'NO'}\n"
            f"Interpretation: {'May affect transcription factor binding and gene regulation' if accessibility_affected else 'Limited regulatory impact expected'}"
        )
        
        ax.text(
            0.5, max(values) * 0.7,
            interp_text,
            ha='center',
            transform=ax.transData,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
            fontsize=10
        )
        
        plt.tight_layout()
        
        # Save figure
        output_file = self.output_dir / f"{variant_info['name']}_chromatin_accessibility.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Chromatin accessibility visualization saved: {output_file}")
        return output_file
    
    def create_summary_figure(self, all_results: List[Dict[str, Any]]) -> Path:
        """
        Create comprehensive summary figure for all variants.
        
        Args:
            all_results: List of all validation results
            
        Returns:
            Path to saved figure
        """
        fig, axes = plt.subplots(len(all_results), 1, figsize=(14, 5 * len(all_results)))
        
        if len(all_results) == 1:
            axes = [axes]
        
        for idx, results in enumerate(all_results):
            ax = axes[idx]
            variant_info = results['variant_info']
            assessment = results['pathogenicity_assessment']
            
            # Prepare data
            gene = variant_info['gene']
            classification = assessment.get('classification', 'Unknown')
            confidence = assessment.get('confidence', 'Unknown')
            
            # Classification color
            class_color = {
                'Likely_pathogenic': '#E63946',
                'Uncertain_significance': '#F77F00',
                'Likely_benign': '#06A77D'
            }.get(classification, '#CCCCCC')
            
            # Create summary visualization
            ax.axis('off')
            
            summary_text = (
                f"═══════════════════════════════════════════════════════════════════════════════\n"
                f"VARIANT: {variant_info['name']}\n"
                f"═══════════════════════════════════════════════════════════════════════════════\n\n"
                f"Gene: {gene}\n"
                f"Position: {variant_info['chromosome']}:{variant_info['position']:,}\n"
                f"Change: {variant_info['ref']} → {variant_info['alt']}\n"
                f"Type: {variant_info['type'].replace('_', ' ').title()}\n"
                f"Disease: {variant_info['disease'].replace('_', ' ')}\n\n"
                f"─────────────────────────────────────────────────────────────────────────────────\n"
                f"ALPHAGENOME CLASSIFICATION: {classification.replace('_', ' ').upper()}\n"
                f"Confidence: {confidence.upper()}\n"
                f"─────────────────────────────────────────────────────────────────────────────────\n\n"
                f"Evidence:\n"
            )
            
            for evidence in assessment.get('evidence', []):
                summary_text += f"  • {evidence['interpretation']}\n"
            
            summary_text += (
                f"\n─────────────────────────────────────────────────────────────────────────────────\n"
                f"ClinVar Classification: {variant_info.get('clinvar_classification', 'N/A')}\n"
                f"ClinVar Votes: {variant_info.get('clnsigconf', 'N/A')}\n"
                f"═══════════════════════════════════════════════════════════════════════════════\n"
            )
            
            ax.text(
                0.5, 0.5,
                summary_text,
                ha='center', va='center',
                fontsize=9,
                family='monospace',
                bbox=dict(
                    boxstyle='round,pad=1.5',
                    facecolor=class_color,
                    alpha=0.2,
                    edgecolor='black',
                    linewidth=3
                )
            )
        
        plt.tight_layout()
        
        # Save figure
        output_file = self.output_dir / 'all_variants_summary.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Summary figure saved: {output_file}")
        return output_file


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Visualize AlphaGenome variant validation results'
    )
    parser.add_argument(
        '--results-dir',
        type=Path,
        default=Path('../outputs'),
        help='Directory containing validation results JSON files'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('../visualizations'),
        help='Output directory for visualizations'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("AlphaGenome Visualization Pipeline")
    print("="*80)
    print(f"Results Directory: {args.results_dir}")
    print(f"Output Directory: {args.output_dir}")
    print("="*80)
    
    # Initialize visualizer
    visualizer = AlphaGenomeVisualizer(args.results_dir, args.output_dir)
    
    # Find all result files
    result_files = list(args.results_dir.glob('*_results.json'))
    
    if not result_files:
        print(f"\nNo result files found in {args.results_dir}")
        print("Please run validate_variants.py first.")
        sys.exit(1)
    
    print(f"\nFound {len(result_files)} result file(s)")
    
    all_results = []
    
    # Process each result file
    for result_file in result_files:
        variant_name = result_file.stem.replace('_results', '')
        print(f"\nProcessing: {variant_name}")
        
        try:
            results = visualizer.load_results(variant_name)
            all_results.append(results)
            
            # Generate visualizations based on available data
            predictions = results['predictions']
            
            if 'RNA_SEQ' in predictions:
                visualizer.visualize_gene_expression(results)
            
            if 'SPLICE_JUNCTIONS' in predictions:
                visualizer.visualize_splicing(results)
            
            if 'ATAC' in predictions:
                visualizer.visualize_chromatin_accessibility(results)
            
        except Exception as e:
            print(f"✗ Failed to visualize {variant_name}: {e}")
            import traceback
            traceback.print_exc()
    
    # Create summary figure
    if all_results:
        print("\nGenerating summary figure...")
        visualizer.create_summary_figure(all_results)
    
    print("\n" + "="*80)
    print("Visualization Complete!")
    print(f"Figures saved to: {args.output_dir}")
    print("="*80)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Batch Variant Effect Processor for VCF files

This script processes multiple VCF files and creates a combined report
of variant effects.
"""

import os
import argparse
import glob
import subprocess
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

def process_single_vcf(vcf_file, fasta_file, genbank_file, output_dir):
    """Process a single VCF file and generate the variant effect prediction."""
    base_name = os.path.basename(vcf_file).split('.')[0]
    output_file = os.path.join(output_dir, f"{base_name}_effects.tsv")
    
    # Call the variant effect predictor script
    cmd = ["python", "vep.py", vcf_file, fasta_file, genbank_file, "-o", output_file]    
    try:
        subprocess.run(cmd, check=True)
        return output_file
    except subprocess.CalledProcessError as e:
        print(f"Error processing {vcf_file}: {e}")
        return None

def combine_results(result_files, output_dir):
    """Combine results from multiple files into a single report."""
    if not result_files:
        print("No result files to combine.")
        return
    
    # Load all results
    all_results = []
    sample_names = []
    
    for result_file in result_files:
        if result_file is None:
            continue
            
        sample_name = os.path.basename(result_file).split('_effects')[0]
        sample_names.append(sample_name)
        
        try:
            df = pd.read_csv(result_file, sep='\t')
            df['sample'] = sample_name
            all_results.append(df)
        except Exception as e:
            print(f"Error loading {result_file}: {e}")
    
    if not all_results:
        print("No valid results to combine.")
        return
    
    # Combine all results
    combined_df = pd.concat(all_results, ignore_index=True)
    
    # Save combined results
    combined_file = os.path.join(output_dir, "combined_effects.tsv")
    combined_df.to_csv(combined_file, sep='\t', index=False)
    
    return combined_df, combined_file


def generate_summary(combined_df, output_dir):
    """Generate summary statistics and plots from the combined results."""
    if combined_df is None or combined_df.empty:
        print("No data to generate summary.")
        return
    
    # Create an output directory for plots
    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    # Summary statistics
    summary_file = os.path.join(output_dir, "variant_summary.txt")
    
    with open(summary_file, 'w') as f:
        # Overall statistics
        f.write("=== OVERALL VARIANT EFFECT SUMMARY ===\n")
        f.write(f"Total variants analyzed: {len(combined_df)}\n\n")
        
        effect_counts = combined_df['EFFECT'].value_counts()
        f.write("Variant effects distribution:\n")
        for effect, count in effect_counts.items():
            percent = count / len(combined_df) * 100
            f.write(f"  {effect}: {count} ({percent:.2f}%)\n")
        
        # Non-synonymous changes
        nonsyn_df = combined_df[combined_df['EFFECT'] == 'nonsynonymous']
        f.write(f"\nNon-synonymous mutations: {len(nonsyn_df)}\n")
        if not nonsyn_df.empty:
            f.write("Top affected genes:\n")
            gene_counts = nonsyn_df['GENE'].value_counts().head(10)
            for gene, count in gene_counts.items():
                if gene != '.':
                    f.write(f"  {gene}: {count}\n")
        
        # Stop gain/loss mutations
        stop_df = combined_df[combined_df['EFFECT'].isin(['stop_gain', 'stop_loss'])]
        f.write(f"\nStop gain/loss mutations: {len(stop_df)}\n")
        if not stop_df.empty:
            f.write("Details:\n")
            for _, row in stop_df.iterrows():
                f.write(f"  {row['sample']}: {row['GENE']} - {row['EFFECT']} at position {row['POS']} ({row['REF']}>{row['ALT']})\n")
        
        # Sample-specific statistics
        f.write("\n=== SAMPLE-SPECIFIC STATISTICS ===\n")
        sample_counts = combined_df['sample'].value_counts()
        for sample, count in sample_counts.items():
            f.write(f"\nSample: {sample}\n")
            f.write(f"Total variants: {count}\n")
            
            sample_df = combined_df[combined_df['sample'] == sample]
            sample_effects = sample_df['EFFECT'].value_counts()
            
            f.write("Effect distribution:\n")
            for effect, eff_count in sample_effects.items():
                percent = eff_count / count * 100
                f.write(f"  {effect}: {eff_count} ({percent:.2f}%)\n")
    
    # Generate plots
    try:
        # Plot 1: Distribution of variant effects across all samples
        plt.figure(figsize=(10, 6))
        sns.countplot(data=combined_df, x='EFFECT', order=effect_counts.index)
        plt.title('Distribution of Variant Effects')
        plt.xlabel('Effect Type')
        plt.ylabel('Count')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'effect_distribution.png'))
        plt.close()
        
        # Plot 2: Variant effects by sample
        plt.figure(figsize=(12, 8))
        effect_by_sample = combined_df.groupby(['sample', 'EFFECT']).size().unstack(fill_value=0)
        effect_by_sample.plot(kind='bar', stacked=True, figsize=(12, 8))
        plt.title('Variant Effects by Sample')
        plt.xlabel('Sample')
        plt.ylabel('Count')
        plt.legend(title='Effect Type')
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'effects_by_sample.png'))
        plt.close()
        
        # Plot 3: Top affected genes (for nonsynonymous mutations)
        if not nonsyn_df.empty:
            plt.figure(figsize=(10, 6))
            top_genes = nonsyn_df['GENE'].value_counts().head(15)
            top_genes = top_genes[top_genes.index != '.']  # Remove intergenic
            if not top_genes.empty:
                sns.barplot(x=top_genes.index, y=top_genes.values)
                plt.title('Top Genes with Nonsynonymous Mutations')
                plt.xlabel('Gene')
                plt.ylabel('Count')
                plt.xticks(rotation=90)
                plt.tight_layout()
                plt.savefig(os.path.join(plots_dir, 'top_affected_genes.png'))
                plt.close()
                
    except Exception as e:
        print(f"Error generating plots: {e}")
    
    print(f"Summary written to {summary_file}")
    print(f"Plots saved to {plots_dir}")

def main():
    parser = argparse.ArgumentParser(description="Process multiple VCF files for variant effect prediction.")
    parser.add_argument("--vcf_dir", required=True, help="Directory containing VCF files")
    parser.add_argument("--fasta", required=True, help="Reference genome FASTA file")
    parser.add_argument("--genbank", required=True, help="GenBank annotation file")
    parser.add_argument("--output_dir", default="variant_results", help="Output directory (default: variant_results)")
    parser.add_argument("--pattern", default="*.vcf", help="Pattern to match VCF files (default: *.vcf)")
    
    args = parser.parse_args()
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    # Find VCF files
    vcf_pattern = os.path.join(args.vcf_dir, args.pattern)
    vcf_files = glob.glob(vcf_pattern)
    
    if not vcf_files:
        print(f"No VCF files found matching pattern {vcf_pattern}")
        return
    
    print(f"Found {len(vcf_files)} VCF files to process.")
    
    # Process each VCF file
    result_files = []
    for vcf_file in vcf_files:
        print(f"Processing {vcf_file}...")
        result_file = process_single_vcf(vcf_file, args.fasta, args.genbank, args.output_dir)
        result_files.append(result_file)
    
    # Combine results and generate summary
    combined_df, combined_file = combine_results(result_files, args.output_dir)
    
    if combined_df is not None:
        print(f"Combined results written to {combined_file}")
        generate_summary(combined_df, args.output_dir)
        print("Processing complete!")
    else:
        print("Failed to generate combined results.")

if __name__ == "__main__":
    main()

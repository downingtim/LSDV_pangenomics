#!/usr/bin/env python3
"""
Clade Differentiation Analysis
This script analyzes variant effect data combined with metadata to identify
mutations that differentiate between clades.
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from scipy.stats import fisher_exact


def load_metadata(metadata_file):
    """Load sample metadata containing clade information."""
    try:
        metadata = pd.read_csv(metadata_file, sep='\t')
        # Ensure required columns exist
        if 'Sample' not in metadata.columns or 'Clade' not in metadata.columns:
            print("Error: Metadata file must contain 'Sample' and 'Clade' columns")
            return None
        
        # Set Sample as index for easier lookup
        metadata = metadata.set_index('Sample')
        return metadata
    except Exception as e:
        print(f"Error loading metadata file: {e}")
        return None


def load_variant_effects(effects_dir, suffix="_effects.tsv"):
    """
    Load variant effect data from multiple files.
    
    Args:
        effects_dir: Directory containing variant effect TSV files
        suffix: Suffix of effect files to identify them
        
    Returns:
        DataFrame with combined variant effects and sample information
    """
    all_data = []
    
    # List files in the directory
    for filename in os.listdir(effects_dir):
        if filename.endswith(suffix):
            file_path = os.path.join(effects_dir, filename)
            sample_name = filename.replace(suffix, "")
            
            try:
                df = pd.read_csv(file_path, sep='\t')
                df['Sample'] = sample_name
                all_data.append(df)
            except Exception as e:
                print(f"Error loading {file_path}: {e}")
    
    if not all_data:
        print(f"No variant effect files found in {effects_dir}")
        return None
    
    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)
    return combined_df


def merge_data(variant_data, metadata):
    """Merge variant effect data with metadata."""
    if variant_data is None or metadata is None:
        return None
    
    # Create a mapping of sample names to their clades
    sample_to_clade = metadata['Clade'].to_dict()
    
    # Apply the mapping to variant data
    variant_data['Clade'] = variant_data['Sample'].map(sample_to_clade)
    
    # Drop rows where clade information is missing
    variants_with_clade = variant_data.dropna(subset=['Clade'])
    
    if len(variants_with_clade) == 0:
        print("No variants with matching clade information found.")
        return None
    
    return variants_with_clade


def identify_clade_mutations(data, min_freq_diff=0.7, min_samples_per_clade=2):
    """
    Identify mutations that differentiate between clades.
    
    Args:
        data: DataFrame with variant and clade information
        min_freq_diff: Minimum frequency difference between clades to consider a mutation distinctive
        min_samples_per_clade: Minimum number of samples required in a clade for analysis
        
    Returns:
        DataFrame with clade-differentiating mutations
    """
    # Skip intergenic mutations - we're interested in coding changes
    coding_data = data[data['EFFECT'] != 'intergenic'].copy()
    
    # Get unique clades and check if we have enough data
    clades = data['Clade'].unique()
    if len(clades) < 2:
        print("Need at least two different clades for comparison.")
        return None
    
    # Count samples per clade
    clade_sample_counts = data.groupby('Clade')['Sample'].nunique()
    valid_clades = clade_sample_counts[clade_sample_counts >= min_samples_per_clade].index.tolist()
    
    if len(valid_clades) < 2:
        print(f"Need at least two clades with {min_samples_per_clade}+ samples each.")
        print("Clade sample counts:", clade_sample_counts.to_dict())
        return None
    
    # Create a unique mutation identifier
    coding_data['mutation_id'] = (
        coding_data['CHROM'] + '_' + 
        coding_data['POS'].astype(str) + '_' + 
        coding_data['REF'] + '>' + 
        coding_data['ALT']
    )
    
    # Add amino acid change for coding variants
    def format_aa_change(row):
        if row['REF_AA'] and row['ALT_AA'] and row['GENE'] != '.':
            return f"{row['GENE']}:{row['REF_AA']}{row['AA_POS']}{row['ALT_AA']}"
        return None
    
    coding_data['aa_change'] = coding_data.apply(format_aa_change, axis=1)
    
    # Calculate mutation frequency in each clade
    results = []
    
    # Get unique mutations
    unique_mutations = coding_data['mutation_id'].unique()
    
    for mutation in unique_mutations:
        mut_data = coding_data[coding_data['mutation_id'] == mutation]
        
        # Skip if only present in one sample overall
        if mut_data['Sample'].nunique() <= 1:
            continue
        
        # Get first row for basic mutation info
        mut_info = mut_data.iloc[0]
        
        # Calculate frequency in each clade
        clade_freqs = {}
        clade_counts = {}
        significance_tests = {}
        
        for clade in valid_clades:
            # Get samples in this clade
            clade_samples = data[data['Clade'] == clade]['Sample'].unique()
            
            # Count samples with this mutation in this clade
            samples_with_mut = mut_data[mut_data['Clade'] == clade]['Sample'].unique()
            
            # Calculate frequency
            freq = len(samples_with_mut) / len(clade_samples) if len(clade_samples) > 0 else 0
            clade_freqs[clade] = freq
            clade_counts[clade] = f"{len(samples_with_mut)}/{len(clade_samples)}"
            
            # For each other clade, calculate significance using Fisher's exact test
            for other_clade in valid_clades:
                if clade >= other_clade:  # Avoid duplicate tests and self-comparisons
                    continue
                
                # Create contingency table
                other_samples = data[data['Clade'] == other_clade]['Sample'].unique()
                other_with_mut = mut_data[mut_data['Clade'] == other_clade]['Sample'].unique()
                
                table = [
                    [len(samples_with_mut), len(clade_samples) - len(samples_with_mut)],
                    [len(other_with_mut), len(other_samples) - len(other_with_mut)]
                ]
                
                # Fisher's exact test
                try:
                    odds_ratio, p_value = fisher_exact(table)
                    significance_tests[f"{clade}_vs_{other_clade}"] = p_value
                except:
                    significance_tests[f"{clade}_vs_{other_clade}"] = 1.0
        
        # Look for frequency differences between clades
        max_diff = 0
        diff_clades = []
        
        for i, clade1 in enumerate(valid_clades):
            for clade2 in valid_clades[i+1:]:
                freq_diff = abs(clade_freqs[clade1] - clade_freqs[clade2])
                if freq_diff > max_diff:
                    max_diff = freq_diff
                    diff_clades = [clade1, clade2]
        
        if max_diff >= min_freq_diff:
            result = {
                'mutation_id': mutation,
                'CHROM': mut_info['CHROM'],
                'POS': mut_info['POS'],
                'REF': mut_info['REF'],
                'ALT': mut_info['ALT'],
                'EFFECT': mut_info['EFFECT'],
                'GENE': mut_info['GENE'],
                'PRODUCT': mut_info['PRODUCT'],
                'REF_AA': mut_info['REF_AA'],
                'ALT_AA': mut_info['ALT_AA'],
                'AA_POS': mut_info['AA_POS'],
                'aa_change': mut_info['aa_change'],
                'max_freq_diff': max_diff,
                'diff_clades': '_vs_'.join(diff_clades),
            }
            
            # Add frequency for each clade
            for clade in valid_clades:
                result[f'{clade}_freq'] = clade_freqs[clade]
                result[f'{clade}_count'] = clade_counts[clade]
            
            # Add p-values for clade comparisons
            for comparison, p_value in significance_tests.items():
                result[f'p_{comparison}'] = p_value
            
            results.append(result)
    
    if not results:
        print("No clade-differentiating mutations found with the specified criteria.")
        return None
    
    # Convert to DataFrame and sort by maximum frequency difference
    result_df = pd.DataFrame(results)
    result_df = result_df.sort_values('max_freq_diff', ascending=False)
    
    return result_df


def analyze_clade_patterns(data):
    """
    Analyze patterns of mutations across different clades.
    Returns summary statistics about clade-specific mutations.
    """
    if data is None or len(data) == 0:
        return None
    
    # Count mutations by effect type for each clade comparison
    effect_counts = data.groupby(['diff_clades', 'EFFECT']).size().unstack(fill_value=0)
    
    # Count genes with differentiating mutations
    gene_counts = data.groupby('diff_clades')['GENE'].nunique()
    
    # Identify top genes with differentiating mutations
    top_genes = data.groupby('GENE').size().sort_values(ascending=False).head(10)
    
    # Calculate average frequency difference by effect type
    avg_diff_by_effect = data.groupby('EFFECT')['max_freq_diff'].mean().sort_values(ascending=False)
    
    return {
        'effect_counts': effect_counts,
        'gene_counts': gene_counts,
        'top_genes': top_genes,
        'avg_diff_by_effect': avg_diff_by_effect
    }


def generate_visualizations(data, patterns, output_dir):
    """Generate visualizations for clade differentiation analysis."""
    if data is None or patterns is None:
        return
    
    # Create output directory
    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    # Plot 1: Distribution of mutation effects by clade comparison
    if 'effect_counts' in patterns and not patterns['effect_counts'].empty:
        plt.figure(figsize=(12, 8))
        patterns['effect_counts'].plot(kind='bar', stacked=True)
        plt.title('Distribution of Mutation Effects by Clade Comparison')
        plt.xlabel('Clade Comparison')
        plt.ylabel('Number of Mutations')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'effect_by_clade.png'))
        plt.close()
    
    # Plot 2: Top genes with differentiating mutations
    if 'top_genes' in patterns and not patterns['top_genes'].empty:
        plt.figure(figsize=(10, 6))
        patterns['top_genes'].plot(kind='bar')
        plt.title('Top Genes with Clade-Differentiating Mutations')
        plt.xlabel('Gene')
        plt.ylabel('Number of Differentiating Mutations')
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'top_differentiating_genes.png'))
        plt.close()
    
    # Plot 3: Heatmap of mutation frequencies by clade
    try:
        # Extract frequency columns
        freq_cols = [col for col in data.columns if col.endswith('_freq')]
        if freq_cols:
            # Prepare data for heatmap
            clade_names = [col.replace('_freq', '') for col in freq_cols]
            heatmap_data = data[freq_cols].copy()
            heatmap_data.columns = clade_names
            
            # Create mutation labels
            if 'aa_change' in data.columns:
                mutation_labels = data['aa_change'].fillna(data['mutation_id'])
            else:
                mutation_labels = data['mutation_id']
            
            # Top 20 mutations for readability
            top_mutations = data.head(20)
            top_heatmap_data = heatmap_data.iloc[:20]
            top_labels = mutation_labels.iloc[:20]
            
            plt.figure(figsize=(10, 12))
            sns.heatmap(top_heatmap_data, annot=True, cmap='YlOrRd', 
                      yticklabels=top_labels, fmt='.2f', cbar_kws={'label': 'Frequency'})
            plt.title('Mutation Frequency by Clade (Top 20 Differentiating Mutations)')
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, 'mutation_frequency_heatmap.png'))
            plt.close()
    except Exception as e:
        print(f"Error generating heatmap: {e}")


def main():
    parser = argparse.ArgumentParser(description="Analyze mutations that differentiate between clades.")
    parser.add_argument("--effects_dir", required=True, help="Directory containing variant effect TSV files")
    parser.add_argument("--metadata", required=True, help="CSV file with sample metadata including clade information")
    parser.add_argument("--output", default="clade_analysis", help="Output directory (default: clade_analysis)")
    parser.add_argument("--min_freq_diff", type=float, default=0.7, 
                        help="Minimum frequency difference to consider a mutation distinctive (default: 0.7)")
    parser.add_argument("--min_samples", type=int, default=2, 
                        help="Minimum samples per clade for analysis (default: 2)")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Load data
    print("Loading metadata...")
    metadata = load_metadata(args.metadata)
    
    print("Loading variant effect data...")
    variant_data = load_variant_effects(args.effects_dir)
    
    if metadata is None or variant_data is None:
        print("Failed to load required data. Exiting.")
        return
    
    # Merge data
    print("Merging metadata with variant data...")
    merged_data = merge_data(variant_data, metadata)
    
    if merged_data is None:
        print("Failed to merge data. Exiting.")
        return
    
    # Identify clade-differentiating mutations
    print("Identifying mutations that differentiate between clades...")
    clade_mutations = identify_clade_mutations(
        merged_data, 
        min_freq_diff=args.min_freq_diff,
        min_samples_per_clade=args.min_samples
    )
    
    if clade_mutations is None:
        print("No clade-differentiating mutations found.")
        return
    
    # Save results
    output_file = os.path.join(args.output, "clade_differentiating_mutations.tsv")
    clade_mutations.to_csv(output_file, sep='\t', index=False)
    print(f"Saved clade-differentiating mutations to {output_file}")
    
    # Analyze patterns
    print("Analyzing patterns of clade-specific mutations...")
    patterns = analyze_clade_patterns(clade_mutations)
    
    # Generate visualizations
    print("Generating visualizations...")
    generate_visualizations(clade_mutations, patterns, args.output)
    
    # Generate summary report
    summary_file = os.path.join(args.output, "clade_analysis_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("=== CLADE DIFFERENTIATION ANALYSIS SUMMARY ===\n\n")
        
        # Basic stats
        f.write(f"Total variant effect data loaded: {len(variant_data)} variants\n")
        f.write(f"Samples with clade information: {merged_data['Sample'].nunique()}\n")
        f.write(f"Clades in analysis: {', '.join(merged_data['Clade'].unique())}\n")
        f.write(f"Clade-differentiating mutations found: {len(clade_mutations)}\n\n")
        
        # Effect counts
        if 'effect_counts' in patterns and not patterns['effect_counts'].empty:
            f.write("Mutation effects by clade comparison:\n")
            f.write(patterns['effect_counts'].to_string() + "\n\n")
        
        # Gene counts
        if 'gene_counts' in patterns and not patterns['gene_counts'].empty:
            f.write("Number of genes with differentiating mutations by clade comparison:\n")
            f.write(patterns['gene_counts'].to_string() + "\n\n")
        
        # Top genes
        if 'top_genes' in patterns and not patterns['top_genes'].empty:
            f.write("Top genes with differentiating mutations:\n")
            for gene, count in patterns['top_genes'].items():
                gene_mutations = clade_mutations[clade_mutations['GENE'] == gene]
                f.write(f"  {gene} ({count} mutations)\n")
                
                # List top mutations for this gene
                for _, row in gene_mutations.iterrows():
                    mut_str = f"{row['REF_AA']}{row['AA_POS']}{row['ALT_AA']}" if row['aa_change'] else f"{row['REF']}>{row['ALT']} at {row['POS']}"
                    f.write(f"    - {mut_str} ({row['EFFECT']}), differentiates {row['diff_clades'].replace('_vs_', ' from ')}\n")
                f.write("\n")
        
        # Key differentiating mutations
        f.write("Top 10 clade-differentiating mutations:\n")
        for i, (_, row) in enumerate(clade_mutations.head(10).iterrows()):
            f.write(f"{i+1}. ")
            if row['aa_change']:
                f.write(f"{row['aa_change']} ({row['EFFECT']})\n")
            else:
                f.write(f"{row['mutation_id']} ({row['EFFECT']})\n")
                
            f.write(f"   Differentiates {row['diff_clades'].replace('_vs_', ' from ')}\n")
            
            # Add frequency information for each clade
            for clade in metadata['Clade'].unique():
                if f'{clade}_freq' in row:
                    f.write(f"   {clade}: {row[f'{clade}_freq']:.2f} ({row[f'{clade}_count']})\n")
            
            # Add p-value if available
            p_cols = [col for col in row.index if col.startswith('p_')]
            for p_col in p_cols:
                clades = p_col.replace('p_', '')
                f.write(f"   p-value ({clades}): {row[p_col]:.4f}\n")
            
            f.write("\n")
    
    print(f"Analysis complete! Summary written to {summary_file}")


if __name__ == "__main__":
    main()

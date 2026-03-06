#!/usr/bin/env python3
"""
Variant Effect Predictor for VCF files
This script predicts the effect of variants (especially SNPs) on protein-coding genes
using a reference genome and GenBank annotation file.
"""

import os
import re
import argparse
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import gzip


def parse_vcf(vcf_file):
    """Parse a VCF file and extract variants."""
    variants = []
    
    if vcf_file.endswith(".gz"):
        f = gzip.open(vcf_file, "rt")
    else:
        f = open(vcf_file)

            
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            
            # For simplicity, we'll only handle simple SNPs
            if len(ref) == 1 and len(alt) == 1 and ',' not in alt:
                variants.append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt
                })
    
    return variants


def parse_genbank(gb_file):
    """
    Parse a GenBank file and extract gene annotations.
    Returns a dictionary of genes with their positions and coding sequences.
    """
    genes = {}
    record = next(SeqIO.parse(gb_file, "genbank"))
    
    for feature in record.features:
        if feature.type == "CDS":
            # Extract gene information
            if 'gene' in feature.qualifiers:
                gene_name = feature.qualifiers['gene'][0]
            elif 'locus_tag' in feature.qualifiers:
                gene_name = feature.qualifiers['locus_tag'][0]
            else:
                # Skip features without gene names
                continue
                
            # Get the location of the gene
            start = int(feature.location.start) + 1  # Convert to 1-based
            end = int(feature.location.end)
            strand = 1 if feature.location.strand == 1 else -1
            
            # Get the coding sequence
            cds = feature.extract(record.seq)
            
            # Extract protein product information
            product = "Unknown protein"
            if 'product' in feature.qualifiers:
                product = feature.qualifiers['product'][0]
            
            # Get protein translation if available
            translation = None
            if 'translation' in feature.qualifiers:
                translation = feature.qualifiers['translation'][0]
            
            genes[gene_name] = {
                'start': start,
                'end': end,
                'strand': strand,
                'cds': str(cds),
                'product': product,
                'translation': translation
            }
    
    return genes, record.seq


def load_genome(fasta_file):
    """Load the reference genome from a FASTA file."""
    for record in SeqIO.parse(fasta_file, "fasta"):
        return record.id, str(record.seq)
    return None, None


def identify_variant_location(variant, genes):
    """Identify which gene (if any) a variant falls within."""
    chrom = variant['chrom']
    pos = variant['pos']
    
    for gene_name, gene_info in genes.items():
        if gene_info['start'] <= pos <= gene_info['end']:
            return gene_name, gene_info
    
    return None, None


def predict_effect(variant, gene_info, genome_seq):
    """
    Predict the effect of a variant on a gene.
    Returns the type of mutation (synonymous/nonsynonymous) and the amino acid change.
    """
    pos = variant['pos']
    ref_base = variant['ref']
    alt_base = variant['alt']
    
    # Make sure the variant position is within the gene
    if not (gene_info['start'] <= pos <= gene_info['end']):
        return "intergenic", None, None, None, None
    
    # Adjust position to be relative to gene start (0-based)
    gene_start = gene_info['start']
    gene_end = gene_info['end']
    strand = gene_info['strand']
    
    # Position in the gene (0-based)
    rel_pos = pos - gene_start
    
    # For genes on the minus strand, we need to complement the bases
    # and calculate the position differently
    if strand == -1:
        rel_pos = gene_end - pos
        ref_base = str(Seq(ref_base).complement())
        alt_base = str(Seq(alt_base).complement())
    
    # Calculate codon position and offset within codon
    codon_pos = rel_pos // 3
    offset = rel_pos % 3
    
    # Extract the codon from the CDS
    cds = gene_info['cds']
    codon_start = codon_pos * 3
    
    # Make sure the codon position is valid
    if codon_start + 2 >= len(cds):
        return "invalid_position", None, None, None, None
    
    ref_codon = cds[codon_start:codon_start+3]
    
    # Create the mutated codon
    alt_codon_list = list(ref_codon)
    alt_codon_list[offset] = alt_base
    alt_codon = ''.join(alt_codon_list)
    
    # Translate the codons
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    
    try:
        ref_aa = standard_table.forward_table.get(ref_codon, '*')
        alt_aa = standard_table.forward_table.get(alt_codon, '*')
    except:
        return "invalid_codon", ref_codon, alt_codon, None, None
    
    # Determine the effect
    if ref_aa == alt_aa:
        effect = "synonymous"
    elif alt_aa == '*':
        effect = "stop_gain"
    elif ref_aa == '*':
        effect = "stop_loss"
    else:
        effect = "nonsynonymous"
    
    return effect, ref_codon, alt_codon, ref_aa, alt_aa


def main(vcf_file, fasta_file, genbank_file, output_file):
    """Main function to run the variant effect prediction."""
    
    # Load reference genome
    genome_id, genome_seq = load_genome(fasta_file)
    
    # Load gene annotations
    genes, gb_seq = parse_genbank(genbank_file)
    
    # Parse VCF file
    variants = parse_vcf(vcf_file)
    
    # Initialize counters
    stats = Counter()
    
    # Process each variant
    results = []
    
    for variant in variants:
        # Find which gene (if any) the variant falls within
        gene_name, gene_info = identify_variant_location(variant, genes)
        
        if gene_name is None:
            # Variant is not in a gene
            effect = "intergenic"
            result = {
                'chrom': variant['chrom'],
                'pos': variant['pos'],
                'ref': variant['ref'],
                'alt': variant['alt'],
                'effect': effect,
                'gene': None,
                'product': None,
                'ref_codon': None,
                'alt_codon': None,
                'ref_aa': None,
                'alt_aa': None,
                'codon_pos': None,
                'aa_pos': None
            }
        else:
            # Variant is in a gene, predict its effect
            pos = variant['pos']
            rel_pos = pos - gene_info['start']
            codon_pos = rel_pos // 3
            aa_pos = codon_pos + 1  # 1-indexed for amino acid position
            
            effect, ref_codon, alt_codon, ref_aa, alt_aa = predict_effect(variant, gene_info, genome_seq)
            
            result = {
                'chrom': variant['chrom'],
                'pos': variant['pos'],
                'ref': variant['ref'],
                'alt': variant['alt'],
                'effect': effect,
                'gene': gene_name,
                'product': gene_info['product'],
                'ref_codon': ref_codon,
                'alt_codon': alt_codon,
                'ref_aa': ref_aa,
                'alt_aa': alt_aa,
                'codon_pos': codon_pos + 1,  # 1-indexed for codon position
                'aa_pos': aa_pos
            }
        
        # Update statistics
        stats[effect] += 1
        results.append(result)
    
    # Write results to output file
    with open(output_file, 'w') as f:
        # Write header
        f.write("CHROM\tPOS\tREF\tALT\tEFFECT\tGENE\tPRODUCT\tREF_CODON\tALT_CODON\tREF_AA\tALT_AA\tCODON_POS\tAA_POS\n")
        
        # Write variant effects
        for result in results:
            f.write(f"{result['chrom']}\t{result['pos']}\t{result['ref']}\t{result['alt']}\t{result['effect']}\t")
            f.write(f"{result['gene'] or '.'}\t{result['product'] or '.'}\t")
            f.write(f"{result['ref_codon'] or '.'}\t{result['alt_codon'] or '.'}\t")
            f.write(f"{result['ref_aa'] or '.'}\t{result['alt_aa'] or '.'}\t")
            f.write(f"{result['codon_pos'] or '.'}\t{result['aa_pos'] or '.'}\n")
    
    # Generate summary
    print(f"Variant Effect Predictor: Summary for {vcf_file}")
    print(f"Total variants processed: {len(variants)}")
    print("\nVariant effects:")
    for effect, count in stats.most_common():
        print(f"  {effect}: {count} ({count/len(variants)*100:.1f}%)")
    
    print(f"\nResults written to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict the effect of variants on genes.")
    parser.add_argument("vcf", help="Input VCF file")
    parser.add_argument("fasta", help="Reference genome FASTA file")
    parser.add_argument("genbank", help="GenBank annotation file")
    parser.add_argument("-o", "--output", default="variant_effects.tsv",
                        help="Output file name (default: variant_effects.tsv)")
    
    args = parser.parse_args()
    
    main(args.vcf, args.fasta, args.genbank, args.output)

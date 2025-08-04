# LSDV_pangenomics
LSDV pangenome variation graph analysis, read mapping &amp; interpretation

# Virus read mapping & variant calling Pipeline
A comprehensive bioinformatics pipeline for mapping virus (principally Lumpy Skin Disease Virus, LSDV) sequencing reads to reference genomes and calling variants using multiple approaches.

# Overview

This pipeline supports mapping reads to different reference genome sets (1, 3, 6, or 121 genomes) and performs variant calling using both traditional alignment methods and graph-based approaches. It integrates multiple tools including minimap2, vg (variation graph), freebayes, and bcftools to provide comprehensive variant analysis.

Features
•	Multi-reference support: Works with 1, 3, 6, or 121 genome references
•	Graph-based analysis: Includes vg (variation graph) support for single genome reference
•	Multiple variant callers: FreeBayes, BCFtools, and vg for comprehensive variant detection
•	Quality control: Integrated qualimap and coverage analysis
•	SLURM compatibility: Designed for HPC cluster execution
•	Comprehensive output: BAM files, VCF files, quality reports, and consensus sequences

Requirements
Software Dependencies
 •	Alignment tools: 
  o	minimap2
  o	vg (variation graph toolkit)
  o	samtools (version 1.9+)

 •	Variant calling: 
  o	freebayes
  o	bcftools
  o	bgzip
  o	tabix

 •	Quality control: 
  o	qualimap (version 2.2.1+)

 •	Preprocessing (assumed to be completed): 
  o	fastp (for read trimming)
  o	kraken (for validation)

# System Requirements
 •	Linux/Unix environment
 •	SLURM workload manager (for cluster execution)
 •	Recommended: 20+ CPU cores, 50GB+ RAM per job

# Input Files Structure

The pipeline expects the following directory structure:
/mnt/lustre/RDS-archive/downing/LSDV/
├── KRAKEN_VALID_FILES/     # Processed FASTQ files
└── FASTP_FILES/            # Alternative FASTQ location

# Installation
1.	Clone this repository:

   git clone https://github.com/downingtim/lsdv-pipeline.git
   cd lsdv-pipeline

2.	Make the script executable: chmod +x merged_lsdv_pipeline.sh

3.	Ensure all dependencies are installed and in your PATH.

# Reference Genome Setup

e.g.fFor 1-genome analysis (with vg support):

# Build variation graph (run once)
pggb -i 1genomes.fasta -m -S -o LSDV3 -t 36 -p 90 -s 1k -n 121
vg convert -t 12 -g LSDV1.gfa -v > LSDV1.vg
vg autoindex -p LSDV1 -g LSDV1.gfa -t 22
samtools faidx 1genomes.fasta
vg convert LSDV1.xg -p > LSDV1.pg

For multi-genome analysis:
Ensure your reference FASTA files are named correctly:
•	3 genomes: 3genomes.fasta
•	6 genomes: 6genomes.fasta
•	121 genomes: LSDV_genomes.fasta

#  Usage

Basic Usage
# Single sample with 1 genome (default, includes vg analysis)
./merged_lsdv_pipeline.sh SampleName

# Single sample with specific reference type
./merged_lsdv_pipeline.sh SampleName 3genome
./merged_lsdv_pipeline.sh SampleName 6genome
./merged_lsdv_pipeline.sh SampleName 121genome

# Batch Processing with SLURM

Create a sample list file (acc_list.txt) with one sample name per line, then:
# Process all samples with 1 genome
more acc_list.txt | perl -e 'while(<>){ chomp; print "sbatch -c 20 --mem=50G -J $_.job merged_lsdv_pipeline.sh $_\n"; }' | bash

# Process all samples with 6 genomes
more acc_list.txt | perl -e 'while(<>){ chomp; print "sbatch -c 20 --mem=50G -J $_.job merged_lsdv_pipeline.sh $_ 6genome\n"; }' | bash

Parameters
 •	sample_name: Base name of the sample (required)
 •	reference_type: Reference genome set to use (optional) 
  o	1genome (default): Single reference with vg support
  o	3genome: Three reference genomes
  o	6genome: Six reference genomes
  o	121genome: 121 reference genomes

Pipeline Workflow

 1. Setup and Configuration
  •	Validates input parameters
  •	Configures reference genomes and file paths
  •	Sets up conditional processing flags

 2. Read Mapping
  •	minimap2: Standard alignment for all reference types
  •	vg giraffe: Graph-based alignment with Giraffe
  •	Supports both paired-end and single-end reads

 3. BAM Processing
  •	SAM to BAM conversion
  •	Coordinate sorting
  •	Duplicate marking and removal
  •	Indexing

 4. Quality Control
  •	qualimap: Comprehensive BAM quality assessment
  •	samtools coverage: Coverage statistics
  •	samtools flagstat: Alignment statistics

 5. Variant Calling
  •	vg call: Graph-based variant calling with VG
  •	freebayes: Haplotype-based variant detection
  •	bcftools: Consensus variant calling

 6. Post-processing
  •	VCF normalization and compression
  •	Indexing with tabix
  •	Variant count summaries

 7. Consensus Generation
  •	Reference and alternate consensus sequences
  •	Uses bcftools consensus

#  Output Structure

The pipeline creates the following output directories:
├── GAM_FILES/              # vg alignment files  
├── GAM_FILES2/             # Secondary vg alignments
├── SAM_FILES/              # SAM alignment files
├── BAM_FILES/              # Primary BAM files
├── BAM_FILES2/             # Secondary BAM files  
├── VG_VCF_FILES/           # vg variant calls  
├── FB_VCF_FILES/           # FreeBayes variant calls
├── FB_VCF_FILES2/          # Secondary FreeBayes calls  
├── BCF_VCF_FILES/          # BCFtools variant calls
├── BCF_VCF_FILES2/         # Secondary BCFtools calls  
├── QUALIMAP/               # Quality control reports
├── COVERAGE/               # Coverage statistics
├── CONSENSUS/              # Consensus sequences
├── ERROR_FILES/            # Error logs
├── AUG/                    # Augmented graphs  
├── PACK_FILES/             # Packed coverage  
└── DEPTH_FILES/            # Depth information  

#  Key Output Files

For each sample, the pipeline generates:
 •	Alignments: BAM_FILES/{sample}.rmdup.sorted.bam
 •	Variants: 
  o	FB_VCF_FILES/{sample}.norm.fb.vcf.gz (FreeBayes)
  o	BCF_VCF_FILES/{sample}.norm.bcf.vcf.gz (BCFtools)
  o	VG_VCF_FILES/{sample}.norm.vg.vcf.gz (VG)
 •	Quality: QUALIMAP/{sample}/qualimapReport.html
 •	Coverage: COVERAGE/{sample}.coverage.txt
 •	Consensus: 
  o	CONSENSUS/{sample}.LA.fasta (reference allele)
  o	CONSENSUS/{sample}.LR.fasta (alternate allele)

# Monitoring and Troubleshooting
 Log Files
 •	Check ERROR_FILES/ for detailed error logs
 •	Each processing step generates specific error files
 •	SLURM job logs provide additional debugging information

# Common Issues
 1.	Missing input files: Ensure FASTQ files are in the expected directories
 2.	Reference genome errors: Verify reference files exist and are properly indexed
 3.	Memory issues: Increase SLURM memory allocation for large datasets
 4.	Permission errors: Ensure write permissions for output directories

Debugging
# Check if input files exist
ls -la /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/{sample}*

# Monitor running jobs
squeue -u $USER

# Check specific error logs

tail -f ERROR_FILES/{sample}.*.errors.txt

Performance Considerations
 •	CPU usage: Pipeline uses 10+ cores per job
 •	Memory: 50GB+ recommended for large reference sets
 •	Storage: Ensure sufficient disk space for intermediate files
 •	Network: Consider data locality for large-scale processing

Customization

The script can be customized by modifying:
 •	Reference paths: Update fasta and fasta2 variables
 •	Tool parameters: Adjust variant calling parameters
 •	Resource allocation: Modify SLURM parameters
 •	Output locations: Change directory paths as needed

Citation
If you use this pipeline in your research, please cite:
 •	The original tools used (minimap2, vg, freebayes, bcftools, etc.)
 •	Your LSDV research publication
 •	This pipeline (if published)

Contributing
 1.	Fork the repository
 2.	Create a feature branch
 3.	Make your changes
 4.	Add tests if applicable
 5.	Submit a pull request

Support
For questions or issues:
 •	Open an issue on GitHub
 •	Check the troubleshooting section
 •	Review error logs in ERROR_FILES/

Changelog

 Version 1.0
 •	Initial release 

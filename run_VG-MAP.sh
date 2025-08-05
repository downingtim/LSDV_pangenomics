#!/bin/bash

#################################################################################
# LSDV (Lumpy Skin Disease Virus) Variant Calling Pipeline
#################################################################################
#
# DESCRIPTION:
# This script performs comprehensive variant calling analysis on LSDV sequencing data
# using different reference genome sets (1, 3, 6, or 121 genomes). It supports both
# paired-end and single-end reads and uses multiple variant callers for comparison.
#
# WORKFLOW:
# [0] Index setup from GFA files (commented - run once per reference)
# [1] Quality control checks
# [2] Read mapping using vg map
# [3] BAM file processing and quality assessment
# [4] Variant calling with vg (graph-based)
# [5] Variant calling with freebayes
# [6] Variant calling with BCFtools
# [7] SNP counting and consensus generation
#
# USAGE:
#   ./run.sh <sample_name> <reference_type>
#
# PARAMETERS:
#   sample_name    : Sample identifier (e.g., SRR11470182)
#   reference_type : Reference genome set to use (1, 3, 6, or 121)
#
# EXAMPLE SLURM SUBMISSION:
#   more acc_list.txt | perl -e 'while(<>){ chomp; print "sbatch -c 40 --mem=250G -J $_.job run.sh $_ 3\n" ; }'
#
# REQUIREMENTS:
# - vg toolkit
# - samtools (v1.9+)
# - bcftools
# - freebayes
# - bgzip/tabix
# - qualimap (optional)
#
# OUTPUT DIRECTORIES:
# - GAM_FILES/     : vg alignment files
# - BAM_FILES/     : BAM alignment files
# - VG_VCF_FILES/  : vg variant calls
# - FB_VCF_FILES/  : freebayes variant calls
# - BCF_VCF_FILES/ : bcftools variant calls
# - CONSENSUS/     : consensus sequences
# - COVERAGE/      : coverage statistics
# - ERROR_FILES/   : error logs
#
#################################################################################

# Check if correct number of arguments provided
if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 <sample_name> [reference_type]"
    echo "Reference types: 1 (default), 3, 6, 121"
    exit 1
fi

# Set input parameters
SAMPLE=$1
REF_TYPE=${2:-1}  # Default to single genome if not specified

# Validate reference type
if [[ ! "$REF_TYPE" =~ ^(1|3|6|121)$ ]]; then
    echo "Error: Reference type must be 1, 3, 6, or 121"
    exit 1
fi

echo "Processing sample: $SAMPLE with reference type: $REF_TYPE genomes"

#################################################################################
# CONFIGURATION BASED ON REFERENCE TYPE
#################################################################################

case $REF_TYPE in
    1)
        FASTA_REF="1genomes.fasta"
        VG_PREFIX="LSDV1"
        echo "Using single genome reference"
        ;;
    3)
        FASTA_REF="3genomes.fasta"
        VG_PREFIX="LSDV3"
        echo "Using 3-genome reference"
        ;;
    6)
        FASTA_REF="6genomes.fasta"
        VG_PREFIX="LSDV6"
        echo "Using 6-genome reference"
        ;;
    121)
        FASTA_REF="LSDV_genomes.fasta"
        VG_PREFIX="LSDV2"  # Note: Script 4 uses LSDV2 for 121 genomes
        echo "Using 121-genome reference"
        ;;
esac

#################################################################################
# FILE NAME SETUP
#################################################################################

# Input files
file1="${SAMPLE}_1.fastq"
file2="${SAMPLE}_2.fastq"
trim1="${SAMPLE}_1_trim.fastq"
trim2="${SAMPLE}_2_trim.fastq"
fastp1="${SAMPLE}_1_trim_fastp.fastq"
fastp2="${SAMPLE}_2_trim_fastp.fastq"
fastp3="${SAMPLE}_trim_fastp.fastq"
fastphtml1="${SAMPLE}_1_trim_fastp.html"
fastphtml2="${SAMPLE}_2_trim_fastp.html"
fastphtml3="${SAMPLE}_trim_fastp.html"

# Output files
gam="${SAMPLE}.gam"
gams="${SAMPLE}.stat"
gaf="${SAMPLE}.gaf"
pack="${SAMPLE}.pack"
depth="${SAMPLE}.depth.csv"
bam="${SAMPLE}.bam"
sorted="${SAMPLE}.sorted.bam"
sorted2="${sorted}.2"
sorted3="${sorted}.3"
rmdup="${SAMPLE}.rmdup.sorted.bam"
rmdup2="${SAMPLE}.flagstat"
cov="${SAMPLE}.coverage.txt"
qual="${SAMPLE}.html"
vcf="${SAMPLE}.vg.vcf"
vcfgz="${vcf}.gz"
aug="${SAMPLE}.aug"
snarls="${SAMPLE}.snarls"
fb="${SAMPLE}.fb.vcf"
fbgz="${fb}.gz"
bcf="${SAMPLE}.bcf.vcf"
bcfgz="${bcf}.gz"
normvg="${SAMPLE}.norm.vg.vcf.gz"
normfb="${SAMPLE}.norm.fb.vcf.gz"
normbcf="${SAMPLE}.norm.bcf.vcf.gz"
fasta1="${SAMPLE}.LA.fasta"
fasta2="${SAMPLE}.LR.fasta"

# Display file configuration
echo "File configuration:"
echo "Input files: $file1 $file2 $trim1 $trim2 $fastp1 $fastp2 $fastphtml1 $fastphtml2"

# Check input files
echo "Checking input files..."
ls -lt /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1 2>/dev/null || echo "Warning: $fastp1 not found"
ls -lt /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2 2>/dev/null || echo "Warning: $fastp2 not found"
ls -lt /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3 2>/dev/null || echo "Warning: $fastp3 not found"

#################################################################################
# DIRECTORY SETUP
#################################################################################

echo "Setting up directories..."
mkdir -p AUG GAM_FILES2 ERROR_FILES BAM_FILES QUALIMAP COVERAGE
mkdir -p VG_VCF_FILES FB_VCF_FILES BCF_VCF_FILES CONSENSUS
mkdir -p GAM_FILES GAF_FILES PACK_FILES DEPTH_FILES

#################################################################################
# [0] GRAPH INDEX SETUP (COMMENTED - RUN ONCE PER REFERENCE)
#################################################################################

# Uncomment and run these commands once per reference type:
#
# For 1 genome:
# # Already indexed
#
# For 3 genomes:
# # pggb -i 3genomes.fasta -m -S -o LSDV3 -t 36 -p 90 -s 1k -n 3
# # vg convert -t 12 -g LSDV3.gfa -v > LSDV3.vg
# # vg convert -x LSDV3.vg > LSDV3.xg
# # vg convert LSDV3.xg -p > LSDV3.pg
# # vg autoindex -p LSDV3 -g LSDV3.gfa -t 22
# # samtools faidx 3genomes.fasta
#
# For 6 genomes:
# # pggb -i 6genomes.fasta -m -S -o LSDV6 -t 36 -p 90 -s 1k -n 121
# # vg convert -t 12 -g LSDV6.gfa -v > LSDV6.vg
# # vg convert -x LSDV6.vg > LSDV6.xg
# # vg autoindex -p LSDV6 -g LSDV6.gfa -t 22
# # samtools faidx 6genomes.fasta
# # vg convert LSDV6.xg -p > LSDV6.pg
#
# For 121 genomes:
# # pggb -i LSDV_genomes.fasta -m -S -o LSDV -t 36 -p 90 -s 1k -n 121
# # vg convert -t 12 -g LSDV.gfa -v > LSDV.vg
# # vg convert -x LSDV.vg > LSDV.xg
# # vg autoindex -p LSDV -g LSDV.gfa -t 22
# # samtools faidx LSDV_genomes.fasta
# # vg convert LSDV.xg -p > LSDV.pg

#################################################################################
# [1] QUALITY CONTROL
#################################################################################

echo "=== Quality Control Phase ==="
# QC steps would go here if needed

#################################################################################
# [2] READ MAPPING
#################################################################################

echo "=== Read Mapping Phase ==="
echo "Mapping reads to $VG_PREFIX graph..."

# Set thread count based on reference type
if [ "$REF_TYPE" == "121" ]; then
    THREADS=29
else
    THREADS=22
fi

# Determine input file type and map accordingly
if test -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2; then
    echo "Found paired-end files, performing paired-end mapping..."
    echo "Command: vg map -x ${VG_PREFIX}.xg -g ${VG_PREFIX}.gcsa -t $THREADS -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1 -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2 > GAM_FILES/$gam"
    vg map -x ${VG_PREFIX}.xg -g ${VG_PREFIX}.gcsa -t $THREADS \
        -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1 \
        -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2 \
        > GAM_FILES/$gam 2> ERROR_FILES/${SAMPLE}.gam2.errors.txt
    echo "Done paired-end mapping"
    
elif test -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3; then
    echo "Found single-end file, performing single-end mapping..."
    echo "Command: vg map -x ${VG_PREFIX}.xg -g ${VG_PREFIX}.gcsa -t $THREADS -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3 > GAM_FILES/$gam"
    vg map -x ${VG_PREFIX}.xg -g ${VG_PREFIX}.gcsa -t $THREADS \
        -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3 \
        > GAM_FILES/$gam 2> ERROR_FILES/${SAMPLE}.gam.errors.txt
    echo "Done single-end mapping"
else
    echo "Error: No valid input files found for sample $SAMPLE"
    exit 1
fi

# Generate mapping statistics
echo "Generating mapping statistics..."
vg stats -a GAM_FILES/$gam > GAM_FILES/$gams
echo "Done generating stats"

# Optional: Convert to GAF format (commented out)
# vg convert -G GAM_FILES/$gam ${VG_PREFIX}.xg > GAF_FILES/$gaf 2> ERROR_FILES/${SAMPLE}.gaf.errors.txt
# echo "Converted to GAF format"

#################################################################################
# [3] BAM FILE PROCESSING AND QUALITY ASSESSMENT
#################################################################################

echo "=== BAM Processing Phase ==="
echo "Converting to BAM format..."

# Convert GAM to BAM
if test -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2; then
    echo "Converting paired-end alignment to BAM..."
    vg map -x ${VG_PREFIX}.xg -g ${VG_PREFIX}.gcsa -t 22 \
        -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1 \
        -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2 \
        --surject-to bam > BAM_FILES/$bam 2> ERROR_FILES/${SAMPLE}.surject.errors.txt
elif test -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3; then
    echo "Converting single-end alignment to BAM..."
    vg map -x ${VG_PREFIX}.xg -g ${VG_PREFIX}.gcsa -t 22 \
        -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3 \
        --surject-to bam > BAM_FILES/$bam 2> ERROR_FILES/${SAMPLE}.map.errors.txt
fi

# Validate BAM file
echo "Validating BAM file with samtools quickcheck..."
~/samtools-1.9/samtools quickcheck -v BAM_FILES/$bam

# BAM processing pipeline
echo "Processing BAM file..."
echo "Step 1: Name sorting..."
~/samtools-1.9/samtools sort -n BAM_FILES/$bam -o BAM_FILES/$sorted

echo "Step 2: Fix mate information..."
~/samtools-1.9/samtools fixmate -m BAM_FILES/$sorted BAM_FILES/$sorted2

echo "Step 3: Position sorting..."
~/samtools-1.9/samtools sort BAM_FILES/$sorted2 -o BAM_FILES/$sorted

echo "Step 4: Remove duplicates..."
~/samtools-1.9/samtools markdup -r -s BAM_FILES/$sorted BAM_FILES/$rmdup

echo "Step 5: Indexing final BAM..."
~/samtools-1.9/samtools index BAM_FILES/$rmdup

echo "Step 6: Generating statistics..."
~/samtools-1.9/samtools flagstat BAM_FILES/$rmdup > BAM_FILES/$rmdup2

echo "Step 7: Calculating coverage..."
/usr/bin/samtools coverage BAM_FILES/$rmdup > COVERAGE/$cov 2> ERROR_FILES/${SAMPLE}.cov.errors.txt

# Cleanup intermediate files
echo "Cleaning up intermediate BAM files..."
rm -rf BAM_FILES/$bam BAM_FILES/$sorted BAM_FILES/$sorted2 BAM_FILES/$sorted3

# Optional: Quality assessment with Qualimap
if [ "$REF_TYPE" != "1" ]; then
    echo "Running Qualimap quality assessment..."
    cd QUALIMAP/
    ~/qualimap_v2.2.1/qualimap bamqc -bam ../BAM_FILES/$rmdup -outfile $qual -outdir $SAMPLE -outformat html
    cd ..
fi

#################################################################################
# [4] VARIANT CALLING WITH VG
#################################################################################

echo "=== VG Variant Calling Phase ==="
echo "Performing graph-based variant calling..."

# Determine correct .pg file based on reference type
if [ "$REF_TYPE" == "121" ]; then
    PG_FILE="LSDV.pg"
    VG_FILE="LSDV.vg"
else
    PG_FILE="${VG_PREFIX}.pg"
    VG_FILE="${VG_PREFIX}.vg"
fi

# Augment the graph
echo "Augmenting graph with sample alignments..."
vg augment $PG_FILE GAM_FILES/$gam -m 3 -q 10 -Q 10 -s -A GAM_FILES2/$gam > AUG/$aug

# Compute snarls (variation sites)
echo "Computing snarls..."
vg snarls AUG/$aug > AUG/$snarls

# Create pileup
echo "Creating pileup..."
vg pack -x AUG/$aug -g GAM_FILES2/$gam -Q 10 -o PACK_FILES/$pack

# Call variants
echo "Calling variants with vg..."
echo "Command: vg call -d 1 AUG/$aug -a -r AUG/$snarls -k PACK_FILES/$pack -s $SAMPLE > VG_VCF_FILES/$vcf"
vg call -d 1 AUG/$aug -a -r AUG/$snarls -k PACK_FILES/$pack -s $SAMPLE > VG_VCF_FILES/$vcf 2> ERROR_FILES/${SAMPLE}.vg.errors.txt

# Calculate depth
echo "Calculating depth across graph nodes..."
vg depth $VG_FILE -m 6 -k PACK_FILES/$pack > DEPTH_FILES/$depth 2> ERROR_FILES/${SAMPLE}.depth.errors.txt &

# Compress and index VCF
echo "Compressing and indexing VG VCF..."
rm -rf VG_VCF_FILES/$vcfgz

# Use appropriate bgzip/tabix based on reference type
if [ "$REF_TYPE" == "121" ]; then
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/bgzip -f VG_VCF_FILES/$vcf
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/tabix -f -p vcf VG_VCF_FILES/$vcfgz
else
    bgzip VG_VCF_FILES/$vcf
    tabix -f -p vcf VG_VCF_FILES/$vcfgz
fi

# Normalize variants
echo "Normalizing VG variants..."
bcftools norm -c w -f $FASTA_REF -m-both -Oz -o VG_VCF_FILES/$normvg VG_VCF_FILES/$vcfgz

if [ "$REF_TYPE" == "121" ]; then
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/tabix -f -p vcf VG_VCF_FILES/$normvg
else
    tabix -f -p vcf VG_VCF_FILES/$normvg
fi

#################################################################################
# [5] VARIANT CALLING WITH FREEBAYES
#################################################################################

echo "=== FreeBayes Variant Calling Phase ==="
echo "Calling variants with FreeBayes..."

rm -rf FB_VCF_FILES/$fbgz
~/freebayes/build/freebayes -f $FASTA_REF -F 0.01 -g 2000 -p 1 \
    --min-alternate-count 1 --min-alternate-fraction 0.001 \
    BAM_FILES/$rmdup > FB_VCF_FILES/$fb 2> ERROR_FILES/${SAMPLE}.fb.errors.txt

# Compress and index
if [ "$REF_TYPE" == "121" ]; then
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/bgzip -f FB_VCF_FILES/$fb
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/tabix -f -p vcf FB_VCF_FILES/$fbgz
else
    bgzip FB_VCF_FILES/$fb
    tabix -f -p vcf FB_VCF_FILES/$fbgz
fi

# Normalize FreeBayes variants
echo "Normalizing FreeBayes variants..."
bcftools norm -c w -f $FASTA_REF -m-both -Oz -o FB_VCF_FILES/$normfb FB_VCF_FILES/$fbgz

if [ "$REF_TYPE" == "121" ]; then
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/tabix -f -p vcf FB_VCF_FILES/$normfb
else
    tabix -f -p vcf FB_VCF_FILES/$normfb
fi

#################################################################################
# [6] VARIANT CALLING WITH BCFTOOLS
#################################################################################

echo "=== BCFtools Variant Calling Phase ==="
echo "Calling variants with BCFtools..."

rm -rf BCF_VCF_FILES/$bcfgz
echo "Command: bcftools mpileup -A -Ob -f $FASTA_REF BAM_FILES/$rmdup | bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES/$bcf"
bcftools mpileup -A -Ob -f $FASTA_REF BAM_FILES/$rmdup | bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES/$bcf

# Compress and index
if [ "$REF_TYPE" == "121" ]; then
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/bgzip -f BCF_VCF_FILES/$bcf
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/tabix -f -p vcf BCF_VCF_FILES/$bcfgz
else
    bgzip BCF_VCF_FILES/$bcf
    tabix -f -p vcf BCF_VCF_FILES/$bcfgz
fi

# Normalize BCFtools variants
echo "Normalizing BCFtools variants..."
bcftools norm -c w -f $FASTA_REF -m-both -Oz -o BCF_VCF_FILES/$normbcf BCF_VCF_FILES/$bcfgz

if [ "$REF_TYPE" == "121" ]; then
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/tabix -f -p vcf BCF_VCF_FILES/$normbcf
else
    tabix -f -p vcf BCF_VCF_FILES/$normbcf
fi

#################################################################################
# [7] SNP COUNTING AND CONSENSUS GENERATION
#################################################################################

echo "=== SNP Counting and Consensus Generation ==="

# Count SNPs from each caller
echo "Counting SNPs:"
echo -n "VG SNPs: "
gzip -dc VG_VCF_FILES/$normvg | grep -v -c "#"

echo -n "FreeBayes SNPs: "
gzip -dc FB_VCF_FILES/$normfb | grep -v -c "#"

echo -n "BCFtools SNPs: "
gzip -dc BCF_VCF_FILES/$normbcf | grep -v -c "#"

# Generate consensus sequences
echo "Generating consensus sequences..."
echo "Creating reference allele consensus..."
bcftools consensus -H LA -f $FASTA_REF BCF_VCF_FILES/$normbcf -o CONSENSUS/$fasta1 &

echo "Creating alternate allele consensus..."
bcftools consensus -H LR -f $FASTA_REF BCF_VCF_FILES/$normbcf -o CONSENSUS/$fasta2 &

wait  # Wait for consensus generation to complete

echo "=== Pipeline Complete ==="
echo "Sample $SAMPLE processed successfully with $REF_TYPE genome reference."
echo "Check output directories for results:"
echo "  - VG variants: VG_VCF_FILES/$normvg"
echo "  - FreeBayes variants: FB_VCF_FILES/$normfb"  
echo "  - BCFtools variants: BCF_VCF_FILES/$normbcf"
echo "  - Consensus sequences: CONSENSUS/$fasta1, CONSENSUS/$fasta2"
echo "  - Coverage statistics: COVERAGE/$cov"
echo "  - Error logs: ERROR_FILES/"

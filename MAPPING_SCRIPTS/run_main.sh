#!/bin/bash

#################################################################################
# LSDV (Lumpy Skin Disease Virus) Integrated Variant Calling Pipeline
#################################################################################
#
# DESCRIPTION:
# This unified pipeline performs comprehensive variant calling analysis on LSDV 
# sequencing data using three different mapping approaches:
# 1. Minimap2 - Traditional linear alignment
# 2. VG Giraffe - Fast graph-based alignment 
# 3. VG Map - Comprehensive graph-based alignment
#
# Each approach can use different reference genome sets (1, 3, 6, or 121 genomes)
# and supports both paired-end and single-end reads.
#
# WORKFLOW:
# [0] Parameter validation and setup
# [1] Reference genome configuration
# [2] Directory structure creation
# [3] Read mapping (user-selected method)
# [4] BAM processing and quality assessment
# [5] Variant calling (VG + FreeBayes + BCFtools)
# [6] Variant normalization and indexing
# [7] Statistics and consensus generation
#
# USAGE:
#   ./lsdv_integrated_pipeline.sh <sample_name> <reference_type> <mapping_method>
#
# PARAMETERS:
#   sample_name    : Sample identifier (e.g., SRR11470182)
#   reference_type : Reference genome set (1, 3, 6, or 121)
#   mapping_method : Alignment method (minimap2, giraffe, vgmap)
#
# EXAMPLES:
#   ./lsdv_integrated_pipeline.sh SRR11470182 3 giraffe
#   ./lsdv_integrated_pipeline.sh SRR10394925 121 minimap2
#   ./lsdv_integrated_pipeline.sh SRR12345678 6 vgmap
#
# SLURM BATCH SUBMISSION:
#   # For all samples with giraffe mapping to 3 genomes:
#   more acc_list.txt | perl -e 'while(<>){ chomp; print "sbatch -c 20 --mem=50G -J $_.giraffe.job lsdv_integrated_pipeline.sh $_ 3 giraffe\n"; }'
#   
#   # For all samples with minimap2 mapping to single genome:
#   more acc_list.txt | perl -e 'while(<>){ chomp; print "sbatch -c 20 --mem=50G -J $_.minimap2.job lsdv_integrated_pipeline.sh $_ 1 minimap2\n"; }'
#
# REQUIREMENTS:
# - vg toolkit (for giraffe and vgmap methods)
# - minimap2 (for minimap2 method)
# - samtools (v1.9+)
# - bcftools
# - freebayes
# - bgzip/tabix
# - qualimap (optional)
#
# OUTPUT STRUCTURE:
# MINIMAP2_<REF_TYPE>GENOME/  - Minimap2 results
# GIRAFFE_<REF_TYPE>GENOME/   - VG Giraffe results  
# VGMAP_<REF_TYPE>GENOME/     - VG Map results
# Each containing: SAM_FILES/, BAM_FILES/, *_VCF_FILES/, CONSENSUS/, etc.
#
#################################################################################

# Color codes for output formatting
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')] $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}[WARNING] $1${NC}"
}

print_error() {
    echo -e "${RED}[ERROR] $1${NC}"
}

print_info() {
    echo -e "${BLUE}[INFO] $1${NC}"
}

#################################################################################
# [0] PARAMETER VALIDATION AND SETUP
#################################################################################

# Check if correct number of arguments provided
if [ $# -ne 3 ]; then
    print_error "Incorrect number of arguments"
    echo "Usage: $0 <sample_name> <reference_type> <mapping_method>"
    echo ""
    echo "Parameters:"
    echo "  sample_name    : Sample identifier (e.g., SRR11470182)"
    echo "  reference_type : Reference genome set (1, 3, 6, or 121)"
    echo "  mapping_method : Alignment method (minimap2, giraffe, vgmap)"
    echo ""
    echo "Examples:"
    echo "  $0 SRR11470182 3 giraffe"
    echo "  $0 SRR10394925 121 minimap2"
    echo "  $0 SRR12345678 6 vgmap"
    exit 1
fi

# Set input parameters
SAMPLE=$1
REF_TYPE=$2
MAPPING_METHOD=$3

# Validate reference type
if [[ ! "$REF_TYPE" =~ ^(1|3|6|121)$ ]]; then
    print_error "Reference type must be 1, 3, 6, or 121"
    exit 1
fi

# Validate mapping method
if [[ ! "$MAPPING_METHOD" =~ ^(minimap2|giraffe|vgmap)$ ]]; then
    print_error "Mapping method must be minimap2, giraffe, or vgmap"
    exit 1
fi

print_status "Starting LSDV pipeline for sample: $SAMPLE"
print_info "Reference type: $REF_TYPE genomes"
print_info "Mapping method: $MAPPING_METHOD"

#################################################################################
# [1] REFERENCE GENOME CONFIGURATION
#################################################################################

case $REF_TYPE in
    1)
        FASTA_REF="1genomes.fasta"
        FASTA_ALT="1genome.fasta"  # Alternative naming from script 1
        VG_PREFIX="LSDV1"
        GENOME_COUNT="1"
        print_info "Using single genome reference"
        ;;
    3)
        FASTA_REF="3genomes.fasta"
        FASTA_ALT="3genomes.fasta"
        VG_PREFIX="LSDV3"
        GENOME_COUNT="3"
        print_info "Using 3-genome reference"
        ;;
    6)
        FASTA_REF="6genomes.fasta"
        FASTA_ALT="6genomes.fasta"
        VG_PREFIX="LSDV6"
        GENOME_COUNT="6"
        print_info "Using 6-genome reference"
        ;;
    121)
        FASTA_REF="LSDV_genomes.fasta"
        FASTA_ALT="LSDV_genomes.fasta"
        VG_PREFIX="LSDV"  # For 121 genomes, some scripts use LSDV, others LSDV2
        VG_PREFIX_ALT="LSDV2"  # Alternative prefix used in some scripts
        GENOME_COUNT="121"
        print_info "Using 121-genome reference"
        ;;
esac

#################################################################################
# [2] FILE NAMING SETUP
#################################################################################

# Input file patterns
file1="${SAMPLE}_1.fastq"
file2="${SAMPLE}_2.fastq"
file3="${SAMPLE}.fastq"
trim1="${SAMPLE}_1_trim.fastq"
trim2="${SAMPLE}_2_trim.fastq"
trim3="${SAMPLE}_trim.fastq"
fastp1="${SAMPLE}_1_trim_fastp.fastq"
fastp2="${SAMPLE}_2_trim_fastp.fastq"
fastp3="${SAMPLE}_trim_fastp.fastq"

# Create method-specific output directory
OUTPUT_DIR="${MAPPING_METHOD^^}_${REF_TYPE}GENOME"
print_info "Output directory: $OUTPUT_DIR"

# Output file naming (method-agnostic)
sam="${SAMPLE}.sam"
gam="${SAMPLE}.gam"
gaf="${SAMPLE}.gaf"
pack="${SAMPLE}.pack"
depth="${SAMPLE}.depth.csv"
bam="${SAMPLE}.bam"
sorted="${SAMPLE}.sorted.bam"
rmdup="${SAMPLE}.rmdup.sorted.bam"
cov="${SAMPLE}.coverage.txt"
qual="${SAMPLE}.html"

# VCF file naming
vg_vcf="${SAMPLE}.vg.vcf"
vg_vcfgz="${vg_vcf}.gz"
fb_vcf="${SAMPLE}.fb.vcf"
fb_vcfgz="${fb_vcf}.gz"
bcf_vcf="${SAMPLE}.bcf.vcf"
bcf_vcfgz="${bcf_vcf}.gz"

# Normalized VCF files
norm_vg="${SAMPLE}.norm.vg.vcf.gz"
norm_fb="${SAMPLE}.norm.fb.vcf.gz"
norm_bcf="${SAMPLE}.norm.bcf.vcf.gz"

# Consensus files
consensus_ref="${SAMPLE}.LA.fasta"
consensus_alt="${SAMPLE}.LR.fasta"

# Graph-specific files
aug="${SAMPLE}.aug"
snarls="${SAMPLE}.snarls"
gams="${SAMPLE}.stat"
flg="${SAMPLE}.flagstat"

#################################################################################
# [3] DIRECTORY STRUCTURE CREATION
#################################################################################

print_status "Setting up directory structure..."

# Create method-specific base directory
mkdir -p $OUTPUT_DIR

# Create subdirectories within method-specific directory
cd $OUTPUT_DIR
mkdir -p SAM_FILES BAM_FILES GAM_FILES GAM_FILES2 GAF_FILES
mkdir -p VG_VCF_FILES FB_VCF_FILES BCF_VCF_FILES
mkdir -p AUG PACK_FILES DEPTH_FILES SNARLS
mkdir -p QUALIMAP COVERAGE CONSENSUS
mkdir -p ERROR_FILES STATS

print_info "Created directory structure in $OUTPUT_DIR"

#################################################################################
# [4] INPUT FILE VALIDATION
#################################################################################

print_status "Validating input files..."

INPUT_PATH="/mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES"

# Check for input files and determine read type
PAIRED_END=false
SINGLE_END=false

if [[ -f "$INPUT_PATH/$fastp1" && -f "$INPUT_PATH/$fastp2" ]]; then
    PAIRED_END=true
    print_info "Found paired-end files: $fastp1, $fastp2"
    ls -lt $INPUT_PATH/$fastp1
    ls -lt $INPUT_PATH/$fastp2
elif [[ -f "$INPUT_PATH/$fastp3" ]]; then
    SINGLE_END=true
    print_info "Found single-end file: $fastp3"
    ls -lt $INPUT_PATH/$fastp3
else
    print_error "No valid input files found for sample $SAMPLE"
    print_error "Expected files:"
    print_error "  Paired-end: $INPUT_PATH/$fastp1 and $INPUT_PATH/$fastp2"
    print_error "  Single-end: $INPUT_PATH/$fastp3"
    exit 1
fi

#################################################################################
# GRAPH GENOME CONSTRUCTION COMMANDS (COMMENTED FOR REFERENCE)
#################################################################################

# Uncomment and run these commands once per reference type to build graph indices:

### For 1 genome reference:
# print_status "Building 1-genome graph indices..."
# samtools faidx 1genomes.fasta
# For minimap2 only, no graph construction needed

### For 3 genomes reference:
# print_status "Building 3-genome graph indices..."
# pggb -i 3genomes.fasta -m -S -o LSDV3 -t 36 -p 90 -s 1k -n 3
# vg convert -t 12 -g LSDV3.gfa -v > LSDV3.vg
# vg autoindex -p LSDV3 -g LSDV3.gfa -t 22  # creates xg and gcsa files
# samtools faidx 3genomes.fasta
# vg convert LSDV3.xg -p > LSDV3.pg
# vg snarls LSDV3.pg > LSDV3.snarls
# vg index -j LSDV3.dist -s LSDV3.snarls LSDV3.xg
# 
# # For giraffe method, additional indices:
# vg gbwt --gbz-format -g LSDV3.gbz -G LSDV3.gfa
# vg gbwt -o LSDV3.gbwt -Z LSDV3.gbz
# vg convert -x --drop-haplotypes LSDV3.gbz > LSDV3.xg
# vg minimizer LSDV3.gbz -d LSDV3.dist -o LSDV3.min

### For 6 genomes reference:
# print_status "Building 6-genome graph indices..."
# pggb -i 6genomes.fasta -m -S -o LSDV6 -t 36 -p 90 -s 1k -n 6
# vg convert -t 12 -g LSDV6.gfa -v > LSDV6.vg
# vg autoindex -p LSDV6 -g LSDV6.gfa -t 22
# samtools faidx 6genomes.fasta
# vg convert LSDV6.xg -p > LSDV6.pg
# vg snarls LSDV6.pg > LSDV6.snarls
# vg index -j LSDV6.dist -s LSDV6.snarls LSDV6.xg
#
# # For giraffe method:
# vg gbwt --gbz-format -g LSDV6.gbz -G LSDV6.gfa
# vg gbwt -o LSDV6.gbwt -Z LSDV6.gbz
# vg convert -x --drop-haplotypes LSDV6.gbz > LSDV6.xg
# vg minimizer LSDV6.gbz -d LSDV6.dist -o LSDV6.min

### For 121 genomes reference:
# print_status "Building 121-genome graph indices..."
# pggb -i LSDV_genomes.fasta -m -S -o LSDV -t 36 -p 90 -s 1k -n 121
# vg convert -t 12 -g LSDV.gfa -v > LSDV.vg
# vg autoindex -p LSDV -g LSDV.gfa -t 22
# samtools faidx LSDV_genomes.fasta
# vg convert LSDV.xg -p > LSDV.pg
# vg snarls LSDV.pg > LSDV.snarls
# vg index -j LSDV.dist -s LSDV.snarls LSDV.xg
#
# # For giraffe method:
# vg gbwt --gbz-format -g LSDV.gbz -G LSDV.gfa
# vg gbwt -o LSDV.gbwt -Z LSDV.gbz
# vg convert -x --drop-haplotypes LSDV.gbz > LSDV.xg
# vg minimizer LSDV.gbz -d LSDV.dist -o LSDV.min

#################################################################################
# [5] READ MAPPING - METHOD SELECTION
#################################################################################

print_status "Starting read mapping with $MAPPING_METHOD..."

case $MAPPING_METHOD in
    "minimap2")
        print_info "Using Minimap2 for linear alignment"
        
        if $PAIRED_END; then
            print_info "Performing paired-end mapping with minimap2..."
            minimap2 -ax sr ../$FASTA_ALT \
                $INPUT_PATH/$fastp1 \
                $INPUT_PATH/$fastp2 \
                > SAM_FILES/$sam 2> ERROR_FILES/${SAMPLE}.sam.errors.txt
        elif $SINGLE_END; then
            print_info "Performing single-end mapping with minimap2..."
            minimap2 -ax sr ../$FASTA_ALT \
                $INPUT_PATH/$fastp3 \
                > SAM_FILES/$sam 2> ERROR_FILES/${SAMPLE}.sam.errors.txt
        fi
        
        # Convert SAM to BAM for minimap2
        print_info "Converting SAM to BAM..."
        samtools view -bS SAM_FILES/$sam > BAM_FILES/$bam
        ;;
        
    "giraffe")
        print_info "Using VG Giraffe for fast graph-based alignment"
        
        # Set thread count based on reference type
        if [ "$REF_TYPE" == "121" ]; then
            THREADS=4
        else
            THREADS=22
        fi
        
        if $PAIRED_END; then
            print_info "Performing paired-end mapping with vg giraffe..."
            vg giraffe -Z ../${VG_PREFIX}.gbz -H ../${VG_PREFIX}.gbwt -p \
                -m ../${VG_PREFIX}.min -d ../${VG_PREFIX}.dist -t $THREADS \
                -f $INPUT_PATH/$fastp1 -f $INPUT_PATH/$fastp2 \
                > GAM_FILES/$gam 2> ERROR_FILES/${SAMPLE}.giraffe.errors.txt
                
            # Generate BAM directly
            print_info "Generating BAM output with vg giraffe..."
            if [ "$REF_TYPE" == "121" ]; then
                vg giraffe -x ../${VG_PREFIX}.xg -Z ../${VG_PREFIX}.gbz \
                    -H ../${VG_PREFIX}.gbwt -p -t $THREADS \
                    -m ../${VG_PREFIX}.min -d ../${VG_PREFIX}.dist \
                    -f $INPUT_PATH/$fastp1 -f $INPUT_PATH/$fastp2 \
                    -o BAM > BAM_FILES/$bam
            else
                vg giraffe -Z ../${VG_PREFIX}.gbz -H ../${VG_PREFIX}.gbwt -p \
                    -m ../${VG_PREFIX}.min -d ../${VG_PREFIX}.dist \
                    -f $INPUT_PATH/$fastp1 -f $INPUT_PATH/$fastp2 \
                    -o BAM > BAM_FILES/$bam
            fi
            
        elif $SINGLE_END; then
            print_info "Performing single-end mapping with vg giraffe..."
            vg giraffe -Z ../${VG_PREFIX}.gbz -H ../${VG_PREFIX}.gbwt -p \
                -m ../${VG_PREFIX}.min -d ../${VG_PREFIX}.dist -t $THREADS \
                -f $INPUT_PATH/$fastp3 \
                > GAM_FILES/$gam 2> ERROR_FILES/${SAMPLE}.giraffe.errors.txt
                
            # Generate BAM directly
            if [ "$REF_TYPE" == "121" ]; then
                vg giraffe -x ../${VG_PREFIX}.xg -Z ../${VG_PREFIX}.gbz \
                    -p -t $THREADS -m ../${VG_PREFIX}.min -d ../${VG_PREFIX}.dist \
                    -f $INPUT_PATH/$fastp3 -o BAM > BAM_FILES/$bam
            else
                vg giraffe -Z ../${VG_PREFIX}.gbz -H ../${VG_PREFIX}.gbwt -p \
                    -m ../${VG_PREFIX}.min -d ../${VG_PREFIX}.dist \
                    -f $INPUT_PATH/$fastp3 -o BAM > BAM_FILES/$bam
            fi
        fi
        
        # Generate mapping statistics
        vg stats -a GAM_FILES/$gam > GAM_FILES/$gams
        ;;
        
    "vgmap")
        print_info "Using VG Map for comprehensive graph-based alignment"
        
        # Set thread count based on reference type
        if [ "$REF_TYPE" == "121" ]; then
            THREADS=29
        else
            THREADS=22
        fi
        
        if $PAIRED_END; then
            print_info "Performing paired-end mapping with vg map..."
            vg map -x ../${VG_PREFIX}.xg -g ../${VG_PREFIX}.gcsa -t $THREADS \
                -f $INPUT_PATH/$fastp1 -f $INPUT_PATH/$fastp2 \
                > GAM_FILES/$gam 2> ERROR_FILES/${SAMPLE}.vgmap.errors.txt
                
            # Generate BAM with surjection
            print_info "Converting to BAM with surjection..."
            vg map -x ../${VG_PREFIX}.xg -g ../${VG_PREFIX}.gcsa -t 22 \
                -f $INPUT_PATH/$fastp1 -f $INPUT_PATH/$fastp2 \
                --surject-to bam > BAM_FILES/$bam 2> ERROR_FILES/${SAMPLE}.surject.errors.txt
                
        elif $SINGLE_END; then
            print_info "Performing single-end mapping with vg map..."
            vg map -x ../${VG_PREFIX}.xg -g ../${VG_PREFIX}.gcsa -t $THREADS \
                -f $INPUT_PATH/$fastp3 \
                > GAM_FILES/$gam 2> ERROR_FILES/${SAMPLE}.vgmap.errors.txt
                
            # Generate BAM with surjection
            vg map -x ../${VG_PREFIX}.xg -g ../${VG_PREFIX}.gcsa -t 22 \
                -f $INPUT_PATH/$fastp3 \
                --surject-to bam > BAM_FILES/$bam 2> ERROR_FILES/${SAMPLE}.surject.errors.txt
        fi
        
        # Generate mapping statistics
        vg stats -a GAM_FILES/$gam > GAM_FILES/$gams
        ;;
esac

print_status "Read mapping completed"

#################################################################################
# [6] BAM PROCESSING AND QUALITY ASSESSMENT
#################################################################################

print_status "Processing BAM files..."

# Validate BAM file
print_info "Validating BAM file..."
~/samtools-1.9/samtools quickcheck -v BAM_FILES/$bam

# BAM processing pipeline
print_info "Step 1: Name sorting..."
sorted2="${sorted}.2"
~/samtools-1.9/samtools sort -n BAM_FILES/$bam -o BAM_FILES/$sorted

print_info "Step 2: Fix mate information..."
~/samtools-1.9/samtools fixmate -m BAM_FILES/$sorted BAM_FILES/$sorted2

print_info "Step 3: Position sorting..."
~/samtools-1.9/samtools sort BAM_FILES/$sorted2 -o BAM_FILES/$sorted

print_info "Step 4: Remove duplicates..."
~/samtools-1.9/samtools markdup -r -s BAM_FILES/$sorted BAM_FILES/$rmdup

print_info "Step 5: Indexing final BAM..."
~/samtools-1.9/samtools index BAM_FILES/$rmdup

print_info "Step 6: Generating statistics..."
~/samtools-1.9/samtools flagstat BAM_FILES/$rmdup > BAM_FILES/$flg

print_info "Step 7: Calculating coverage..."
/usr/bin/samtools coverage BAM_FILES/$rmdup > COVERAGE/$cov 2> ERROR_FILES/${SAMPLE}.cov.errors.txt

# Cleanup intermediate BAM files
print_info "Cleaning up intermediate files..."
rm -rf SAM_FILES/$sam BAM_FILES/$bam BAM_FILES/$sorted BAM_FILES/$sorted2

# Quality assessment with Qualimap
print_info "Running Qualimap quality assessment..."
cd QUALIMAP/
~/qualimap_v2.2.1/qualimap bamqc -bam ../BAM_FILES/$rmdup -outfile $qual -outdir $SAMPLE -outformat html
cd ..

#################################################################################
# [7] VARIANT CALLING
#################################################################################

print_status "Starting variant calling phase..."

# Determine which bgzip/tabix to use based on reference type
if [ "$REF_TYPE" == "121" ]; then
    BGZIP_CMD="/mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/bgzip"
    TABIX_CMD="/mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/tabix"
else
    BGZIP_CMD="bgzip"
    TABIX_CMD="tabix"
fi

#################################################################################
# [7a] VG VARIANT CALLING (for graph-based methods only)
#################################################################################

if [[ "$MAPPING_METHOD" == "giraffe" || "$MAPPING_METHOD" == "vgmap" ]]; then
    print_info "Performing graph-based variant calling with VG..."
    
    # Determine correct .pg file based on reference type
    if [ "$REF_TYPE" == "121" ]; then
        PG_FILE="../LSDV.pg"
        VG_FILE="../LSDV.vg"
    else
        PG_FILE="../${VG_PREFIX}.pg"
        VG_FILE="../${VG_PREFIX}.vg"
    fi
    
    # Augment the graph
    print_info "Augmenting graph with sample alignments..."
    vg augment $PG_FILE GAM_FILES/$gam -m 3 -q 10 -Q 10 -A GAM_FILES2/$gam > AUG/$aug
    
    # Compute snarls (variation sites)
    print_info "Computing snarls..."
    vg snarls AUG/$aug > AUG/$snarls
    
    # Create pileup
    print_info "Creating pileup..."
    vg pack -x AUG/$aug -g GAM_FILES2/$gam -Q 10 -o PACK_FILES/$pack
    
    # Call variants
    print_info "Calling variants with vg..."
    vg call -d 1 AUG/$aug -a -r AUG/$snarls -k PACK_FILES/$pack -s $SAMPLE \
        > VG_VCF_FILES/$vg_vcf 2> ERROR_FILES/${SAMPLE}.vg.errors.txt
    
    # Calculate depth
    print_info "Calculating depth across graph nodes..."
    vg depth $VG_FILE -m 6 -k PACK_FILES/$pack > DEPTH_FILES/$depth 2> ERROR_FILES/${SAMPLE}.depth.errors.txt &
    
    # Compress and index VCF
    print_info "Compressing and indexing VG VCF..."
    rm -rf VG_VCF_FILES/$vg_vcfgz
    $BGZIP_CMD -f VG_VCF_FILES/$vg_vcf
    $TABIX_CMD -f -p vcf VG_VCF_FILES/$vg_vcfgz
    
    # Normalize variants
    print_info "Normalizing VG variants..."
    bcftools norm -c w -f ../$FASTA_REF -m-both -Oz -o VG_VCF_FILES/$norm_vg VG_VCF_FILES/$vg_vcfgz
    $TABIX_CMD -f -p vcf VG_VCF_FILES/$norm_vg
else
    print_info "Skipping VG variant calling (not applicable for minimap2)"
fi

#################################################################################
# [7b] FREEBAYES VARIANT CALLING
#################################################################################

print_info "Calling variants with FreeBayes..."
rm -rf FB_VCF_FILES/$fb_vcfgz
~/freebayes/build/freebayes -f ../$FASTA_REF -F 0.01 -g 2000 -p 1 \
    --min-alternate-count 1 --min-alternate-fraction 0.001 \
    BAM_FILES/$rmdup > FB_VCF_FILES/$fb_vcf 2> ERROR_FILES/${SAMPLE}.fb.errors.txt

# Compress and index
$BGZIP_CMD -f FB_VCF_FILES/$fb_vcf
$TABIX_CMD -f -p vcf FB_VCF_FILES/$fb_vcfgz

# Normalize FreeBayes variants
print_info "Normalizing FreeBayes variants..."
bcftools norm -c w -f ../$FASTA_REF -m-both -Oz -o FB_VCF_FILES/$norm_fb FB_VCF_FILES/$fb_vcfgz
$TABIX_CMD -f -p vcf FB_VCF_FILES/$norm_fb

#################################################################################
# [7c] BCFTOOLS VARIANT CALLING
#################################################################################

print_info "Calling variants with BCFtools..."
rm -rf BCF_VCF_FILES/$bcf_vcfgz
bcftools mpileup -A -Ob -f ../$FASTA_REF BAM_FILES/$rmdup | \
    bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES/$bcf_vcf

# Compress and index
$BGZIP_CMD -f BCF_VCF_FILES/$bcf_vcf
$TABIX_CMD -f -p vcf BCF_VCF_FILES/$bcf_vcfgz

# Normalize BCFtools variants
print_info "Normalizing BCFtools variants..."
bcftools norm -c w -f ../$FASTA_REF -m-both -Oz -o BCF_VCF_FILES/$norm_bcf BCF_VCF_FILES/$bcf_vcfgz
$TABIX_CMD -f -p vcf BCF_VCF_FILES/$norm_bcf

#################################################################################
# [8] STATISTICS AND CONSENSUS GENERATION
#################################################################################

print_status "Generating statistics and consensus sequences..."

# Count SNPs from each caller
print_info "Counting variants:"
if [[ "$MAPPING_METHOD" == "giraffe" || "$MAPPING_METHOD" == "vgmap" ]]; then
    VG_COUNT=$(gzip -dc VG_VCF_FILES/$norm_vg | grep -v -c "^#")
    echo "VG SNPs: $VG_COUNT"
else
    echo "VG SNPs: N/A (minimap2 method)"
fi

FB_COUNT=$(gzip -dc FB_VCF_FILES/$norm_fb | grep -v -c "^#")
echo "FreeBayes SNPs: $FB_COUNT"

BCF_COUNT=$(gzip -dc BCF_VCF_FILES/$norm_bcf | grep -v -c "^#")
echo "BCFtools SNPs: $BCF_COUNT"

# Generate consensus sequences
print_info "Generating consensus sequences..."
bcftools consensus -H

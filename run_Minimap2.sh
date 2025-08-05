#!/bin/bash

#============================================================================
# LSDV Read Mapping and Variant Calling Pipeline
#============================================================================
# 
# DESCRIPTION:
# This script performs read mapping and variant calling for LSDV (Lumpy Skin Disease Virus)
# samples against different reference genome sets. It supports mapping to 1, 3, 6, or 121 
# genome references and includes both graph-based (vg) and traditional alignment approaches.
#
# USAGE:
# ./merged_lsdv_pipeline.sh <sample_name> [reference_type]
#
# PARAMETERS:
# - sample_name: Base name of the sample (from acc_list.txt)
# - reference_type: Optional. Options are "1genome", "3genome", "6genome", or "121genome" 
#                   Default is "1genome"
#
# SLURM EXECUTION:
# Run all samples from acc_list.txt:
# more acc_list.txt | perl -e 'while(<>){ chomp; print "sbatch -c 20 --mem=50G -J $_.job merged_lsdv_pipeline.sh $_ [reference_type]\n"; }'
#
# PIPELINE STEPS:
# [0] Setup - Configure reference genome based on type
# [1] Read mapping with minimap2 (and vg for 1genome)
# [2] BAM processing (sort, fixmate, markdup, index)
# [3] Quality control with qualimap
# [4] Variant calling with vg (1genome only), freebayes, and bcftools
# [5] VCF normalization and indexing
# [6] Consensus sequence generation
#
# OUTPUT DIRECTORIES:
# - SAM_FILES/, BAM_FILES/, BAM_FILES2/: Alignment files
# - VG_VCF_FILES/, FB_VCF_FILES/, BCF_VCF_FILES/: Variant call files
# - QUALIMAP/, COVERAGE/: Quality control reports
# - CONSENSUS/: Consensus sequences
# - ERROR_FILES/: Error logs
#
#============================================================================

# Get input parameters
SAMPLE=$1
REF_TYPE=${2:-"1genome"}  # Default to 1genome if not specified

# Validate reference type
case $REF_TYPE in
    "1genome"|"3genome"|"6genome"|"121genome")
        echo "Processing sample $SAMPLE with reference type: $REF_TYPE"
        ;;
    *)
        echo "Error: Invalid reference type '$REF_TYPE'. Valid options: 1genome, 3genome, 6genome, 121genome"
        exit 1
        ;;
esac

#============================================================================
# [0] SETUP AND CONFIGURATION
#============================================================================

# Configure file naming based on sample
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

# Configure reference genome based on type
case $REF_TYPE in
    "1genome")
        fasta="1genomes.fasta"
        fasta2="1genome.fasta"
        ;;
    "3genome")
        fasta="3genomes.fasta"
        fasta2="3genomes.fasta"
        ;;
    "6genome")
        fasta="6genomes.fasta"
        fasta2="6genomes.fasta"
        ;;
    "121genome")
        fasta="LSDV_genomes.fasta"
        fasta2="LSDV_genomes.fasta"
        ;;
esac

# Output file naming
sam="${SAMPLE}.sam"
gaf="${SAMPLE}.gaf"
pack="${SAMPLE}.pack"
depth="${SAMPLE}.depth.csv"
bam="${SAMPLE}.bam"
sorted="${SAMPLE}.sorted.bam"
rmdup="${SAMPLE}.rmdup.sorted.bam"
cov="${SAMPLE}.coverage.txt"
qual="${SAMPLE}.html"
echo "Sample: $SAMPLE"
echo "Reference files: $fasta / $fasta2"
echo "Checking input files..."
ls -lt /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1
ls -lt /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2
ls -lt /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3

#============================================================================
# COMMENTED SETUP FOR VG GRAPH CONSTRUCTION (for reference)
#============================================================================
# [0] Set up index from GFA - simpler than it seemed!!
# For 1genome reference, these steps would have been run once:
#pggb -i 1genomes.fasta -m -S -o LSDV3 -t 36 -p 90 -s 1k -n 121
#vg convert -t 12 -g LSDV1.gfa -v > LSDV1.vg
#vg autoindex -p LSDV1 -g LSDV1.gfa -t 22 # creates xg and gcsa files
#samtools faidx 1genomes.fasta
#vg convert LSDV1.xg -p > LSDV1.pg

#============================================================================
# [1] READ MAPPING
#============================================================================

echo "=== READ MAPPING ==="

# Check for paired-end or single-end data
if test -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2; then
    echo "Processing paired-end data"
    
    # Minimap2 mapping
    echo "minimap2 -ax sr $fasta2 /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1 /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2 > SAM_FILES/$sam 2> ERROR_FILES/${SAMPLE}.sam.errors.txt"
    minimap2 -ax sr $fasta2 /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1 /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2 > SAM_FILES/$sam 2> ERROR_FILES/${SAMPLE}.sam.errors.txt
    echo "Done paired end mapping"
    
elif test -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3; then
    echo "Processing single-end data"
    
    minimap2 -ax sr $fasta2 /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3 > SAM_FILES/$sam 2> ERROR_FILES/${SAMPLE}.sam.errors.txt
    fi
    echo "Done single end mapping"
else
    echo "Error: No valid input files found"
    exit 1
fi

#============================================================================
# [2] BAM PROCESSING
#============================================================================

# Convert SAM to BAM
echo "Converting SAM -> BAM"
samtools view -bS SAM_FILES/$sam > BAM_FILES/$bam
# Sort BAM files
echo "Sorting BAM"
sorted2="${sorted}.2"
echo "~/samtools-1.9/samtools sort -n BAM_FILES/$bam -o BAM_FILES/$sorted"
~/samtools-1.9/samtools sort -n BAM_FILES/$bam -o BAM_FILES/$sorted
echo "~/samtools-1.9/samtools fixmate -m BAM_FILES/$sorted BAM_FILES/$sorted2"
~/samtools-1.9/samtools fixmate -m BAM_FILES/$sorted BAM_FILES/$sorted2
echo "~/samtools-1.9/samtools sort BAM_FILES/$sorted2 -o BAM_FILES/$sorted"
~/samtools-1.9/samtools sort BAM_FILES/$sorted2 -o BAM_FILES/$sorted

# Remove duplicates
echo "RMDUP BAM"
echo "samtools markdup -s BAM_FILES/$sorted BAM_FILES/$rmdup"
~/samtools-1.9/samtools markdup -r -s BAM_FILES/$sorted BAM_FILES/$rmdup

# Index BAM files
echo "Indexing BAM"
echo "samtools index BAM_FILES/$rmdup"
~/samtools-1.9/samtools index BAM_FILES/$rmdup

# Generate BAM statistics
echo "Getting BAM stats"
rmdup2="${SAMPLE}.flagstat"
~/samtools-1.9/samtools flagstat BAM_FILES/$rmdup > BAM_FILES/$rmdup2
/usr/bin/samtools coverage BAM_FILES/$rmdup > COVERAGE/$cov 2> ERROR_FILES/${SAMPLE}.cov.errors.txt

# Clean up intermediate files
rm -rf SAM_FILES/$sam BAM_FILES/$bam BAM_FILES/$sorted BAM_FILES/$sorted2 BAM_FILES/$sorted3

#============================================================================
# [3] QUALITY CONTROL
#============================================================================

echo "=== QUALITY CONTROL ==="

cd QUALIMAP/
echo "Doing Qualimap"
echo "~/qualimap_v2.2.1/qualimap bamqc -bam ../BAM_FILES/$rmdup -outfile $qual.html -outdir $SAMPLE -outformat html"
~/qualimap_v2.2.1/qualimap bamqc -bam ../BAM_FILES/$rmdup -outfile $qual.html -outdir $SAMPLE -outformat html
cd ..


#============================================================================
# [4] VARIANT CALLING WITH FREEBAYES
#============================================================================

echo "=== FREEBAYES VARIANT CALLING ==="

fb="${SAMPLE}.fb.vcf"
fbgz="${fb}.gz"
echo "rm -rf FB_VCF_FILES/$fbgz"
rm -rf FB_VCF_FILES/$fbgz
~/freebayes/build/freebayes -f $fasta -F 0.01 -g 2000 -p 1 --min-alternate-count 1 --min-alternate-fraction 0.001 BAM_FILES/$rmdup > FB_VCF_FILES/$fb 2> ERROR_FILES/${SAMPLE}.fb.errors.txt
bgzip FB_VCF_FILES/$fb
tabix -f -p vcf FB_VCF_FILES/$fbgz

# Normalize FB SNPs
echo "Normalise FB SNPs"
normfb="${SAMPLE}.norm.fb.vcf.gz"
echo "bcftools norm -c w -f $fasta -m-both -Oz -o FB_VCF_FILES/$normfb FB_VCF_FILES/$fbgz"
bcftools norm -c w -f $fasta -m-both -Oz -o FB_VCF_FILES/$normfb FB_VCF_FILES/$fbgz
echo "tabix -f -p vcf FB_VCF_FILES/$normfb"
tabix -f -p vcf FB_VCF_FILES/$normfb

# Additional freebayes run with fasta2 (if different from fasta)
if [[ "$fasta" != "$fasta2" ]]; then
    fb2="${SAMPLE}.fb.vcf"
    fbgz2="${fb2}.gz"
    
    echo "rm -rf FB_VCF_FILES/$fbgz2"
    rm -rf FB_VCF_FILES/$fbgz2
    ~/freebayes/build/freebayes -f $fasta2 -F 0.01 -g 2000 -p 1 --min-alternate-count 1 --min-alternate-fraction 0.001 BAM_FILES/$rmdup > FB_VCF_FILES/$fb2 2> ERROR_FILES/${SAMPLE}.fb2.errors.txt
    bgzip FB_VCF_FILES/$fb2
    tabix -f -p vcf FB_VCF_FILES/$fbgz2
        
    echo "Normalise FB SNPs 2"
    normfb2="${SAMPLE}.norm.fb.vcf.gz"
    bcftools norm -c w -f $fasta2 -m-both -Oz -o FB_VCF_FILES/$normfb2 FB_VCF_FILES/$fbgz2
    tabix -f -p vcf FB_VCF_FILES/$normfb2
fi

#============================================================================
# [6] VARIANT CALLING WITH BCFTOOLS
#============================================================================

echo "=== BCFTOOLS VARIANT CALLING ==="

bcf="${SAMPLE}.bcf.vcf"
bcfgz="${bcf}.gz"
echo "rm -rf BCF_VCF_FILES/$bcfgz"
rm -rf BCF_VCF_FILES/$bcfgz
echo "bcftools mpileup -A -Ob -f $fasta BAM_FILES/$rmdup | bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES/$bcf"
bcftools mpileup -A -Ob -f $fasta BAM_FILES/$rmdup | bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES/$bcf
bgzip BCF_VCF_FILES/$bcf
tabix -f -p vcf BCF_VCF_FILES/$bcfgz

# Normalize BCF SNPs
echo "Normalise BCF SNPs"
normbcf="${SAMPLE}.norm.bcf.vcf.gz"
echo "bcftools norm -c w -f $fasta -m-both -Oz -o BCF_VCF_FILES/$normbcf BCF_VCF_FILES/$bcfgz"
bcftools norm -c w -f $fasta -m-both -Oz -o BCF_VCF_FILES/$normbcf BCF_VCF_FILES/$bcfgz
echo "tabix -f -p vcf BCF_VCF_FILES/$normbcf"
tabix -f -p vcf BCF_VCF_FILES/$normbcf

# Additional bcftools run with fasta2 (if different or if using VG)
if [[ "$fasta" != "$fasta2" ]]; then
    bcf2="${SAMPLE}.bcf.vcf"
    bcfgz2="${bcf2}.gz"
    
    echo "rm -rf BCF_VCF_FILES2/$bcfgz2"
    rm -rf BCF_VCF_FILES2/$bcfgz2
    echo "bcftools mpileup -A -Ob -f $fasta2 BAM_FILES2/$rmdup | bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES2/$bcf2"
    bcftools mpileup -A -Ob -f $fasta2 BAM_FILES2/$rmdup | bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES2/$bcf2
    bgzip BCF_VCF_FILES2/$bcf2
    tabix -f -p vcf BCF_VCF_FILES2/$bcfgz2

    echo "rm -rf BCF_VCF_FILES/$bcfgz2"
    rm -rf BCF_VCF_FILES/$bcfgz2
    bcftools mpileup -A -Ob -f $fasta2 BAM_FILES/$rmdup | bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES/$bcf2
    bgzip BCF_VCF_FILES/$bcf2
    tabix -f -p vcf BCF_VCF_FILES/$bcfgz2
fi

#============================================================================
# [7] VARIANT COUNTS SUMMARY
#============================================================================

echo "=== VARIANT COUNTS SUMMARY ==="

echo "FB SNPs"
gzip -dc FB_VCF_FILES/$normfb | grep -v -c "#"

echo "BCF SNPs"
gzip -dc BCF_VCF_FILES/$normbcf | grep -v -c "#"

#============================================================================
# [8] CONSENSUS SEQUENCE GENERATION
#============================================================================

echo "=== CONSENSUS SEQUENCE GENERATION ==="

# Generate consensus sequences from BCF VCF
echo "Making consensus"
fasta1="${SAMPLE}.LA.fasta"
fasta2_consensus="${SAMPLE}.LR.fasta"

echo "bcftools consensus -H LA -f $fasta BCF_VCF_FILES/$normbcf -o CONSENSUS/$fasta1"
bcftools consensus -H LA -f $fasta BCF_VCF_FILES/$normbcf -o CONSENSUS/$fasta1  # ref allele

echo "bcftools consensus -H LR -f $fasta BCF_VCF_FILES/$normbcf -o CONSENSUS/$fasta2_consensus"
bcftools consensus -H LR -f $fasta BCF_VCF_FILES/$normbcf -o CONSENSUS/$fasta2_consensus  # alt allele

#============================================================================
# PIPELINE COMPLETION
#============================================================================

echo "=== PIPELINE COMPLETED ==="
echo "Sample: $SAMPLE"
echo "Reference type: $REF_TYPE"
echo "Output files generated:"
echo "- Alignments: BAM_FILES/$rmdup"
echo "- FreeBayes Variants: FB_VCF_FILES/$normfb"
echo "- BCFtools Variants: BCF_VCF_FILES/$normbcf"
echo "- Quality Reports: QUALIMAP/$SAMPLE/"
echo "- Coverage: COVERAGE/$cov"
echo "- Consensus: CONSENSUS/$fasta1, CONSENSUS/$fasta2_consensus"
echo ""
echo "Check ERROR_FILES/ directory for any error logs."

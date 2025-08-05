#!/bin/bash
#
# LSDV Genomics Pipeline - Unified Variant Calling Script
# 
# USAGE: ./run.sh <sample_id> <reference_type>
# 
# PARAMETERS:
#   sample_id      : Sample identifier (e.g., SRR10394925)
#   reference_type : Reference genome type (1|3|6|121)
#                   1   = 1genomes.fasta   (single genome reference)
#                   3   = 3genomes.fasta   (3-genome reference)
#                   6   = 6genomes.fasta   (6-genome reference)
#                   121 = LSDV_genomes.fasta (121-genome reference)
#
# DESCRIPTION:
#   This pipeline processes LSDV (Lumpy Skin Disease Virus) sequencing data through:
#   [0] Graph genome indexing (commented - run once per reference)
#   [1] Quality control and trimming
#   [2] Read mapping using vg giraffe to graph genomes
#   [3] BAM processing and quality assessment
#   [4] Variant calling with vg (graph-based)
#   [5] Variant calling with freebayes
#   [6] Variant calling with BCFtools
#   [7] SNP counting and statistics
#   [8] Consensus sequence generation
#
# REQUIREMENTS:
#   - vg toolkit for graph genomics
#   - samtools, bcftools, bgzip, tabix
#   - freebayes
#   - qualimap
#   - Input files in KRAKEN_VALID_FILES/ directory
#
# SLURM USAGE:
#   Generate job submissions from acc_list.txt:
#   more acc_list.txt | perl -e 'while(<>){ chomp; print "sbatch -c 20 --mem=50G -J $_.job run.sh $_ <ref_type>\n"; }'
#
# AUTHOR: Bioinformatics Pipeline for LSDV Analysis
# DATE: 2024

# Check arguments
if [ $# -lt 2 ]; then
    echo "Error: Missing arguments"
    echo "Usage: $0 <sample_id> <reference_type>"
    echo "Reference types: 1, 3, 6, or 121"
    exit 1
fi

# Parse command line arguments
SAMPLE_ID=$1
REF_TYPE=$2

# Validate reference type and set parameters
case $REF_TYPE in
    1)
        FASTA="1genomes.fasta"
        VG_PREFIX="LSDV1"
        AUG_SUFFIX="_LSDV6_aug.pg"  # Note: keeping original naming for compatibility
        ;;
    3)
        FASTA="3genomes.fasta"
        VG_PREFIX="LSDV3"
        AUG_SUFFIX="_LSDV3_aug.pg"
        ;;
    6)
        FASTA="6genomes.fasta"
        VG_PREFIX="LSDV6"
        AUG_SUFFIX="_LSDV6_aug.pg"
        ;;
    121)
        FASTA="LSDV_genomes.fasta"
        VG_PREFIX="LSDV"
        AUG_SUFFIX="_LSDV6_aug.pg"  # Note: keeping original naming for compatibility
        ;;
    *)
        echo "Error: Invalid reference type '$REF_TYPE'"
        echo "Valid options: 1, 3, 6, or 121"
        exit 1
        ;;
esac

# Set up file naming variables
file1="${SAMPLE_ID}_1.fastq"
file2="${SAMPLE_ID}_2.fastq"
file3="${SAMPLE_ID}.fastq"
trim1="${SAMPLE_ID}_1_trim.fastq"
trim2="${SAMPLE_ID}_2_trim.fastq"
trim3="${SAMPLE_ID}_trim.fastq"
fastp1="${SAMPLE_ID}_1_trim_fastp.fastq"
fastp2="${SAMPLE_ID}_2_trim_fastp.fastq"
fastp3="${SAMPLE_ID}_trim_fastp.fastq"
fastphtml1="${SAMPLE_ID}_1_trim_fastp.html"
fastphtml2="${SAMPLE_ID}_2_trim_fastp.html"
fastphtml3="${SAMPLE_ID}_trim_fastp.html"

echo "Processing sample: $SAMPLE_ID with reference type: $REF_TYPE ($FASTA)"
echo "Files: $file1 $file2 $trim1 $trim2 $fastp1 $fastp2 $fastphtml1 $fastphtml2"

# Create necessary directories (uncomment as needed)
# mkdir -p AUG GAM_FILES GAM_FILES2 ERROR_FILES PACK_FILES DEPTH_FILES
# mkdir -p BAM_FILES QUALIMAP COVERAGE VG_VCF_FILES FB_VCF_FILES BCF_VCF_FILES CONSENSUS

# [0] Set up index from GFA - Run once per reference type (COMMENTED - uncomment to rebuild indices)
################################################################################
# echo "Setting up graph genome index for $VG_PREFIX..."
# pggb -i $FASTA -m -S -o $VG_PREFIX -t 36 -p 90 -s 1k -n $(echo $REF_TYPE)
# vg convert -t 12 -g ${VG_PREFIX}.gfa -v > ${VG_PREFIX}.vg
# samtools faidx $FASTA
# 
# # Build GBWT of the walks in the GFA file
# vg gbwt --gbz-format -g ${VG_PREFIX}.gbz -G ${VG_PREFIX}.gfa
# vg gbwt -o ${VG_PREFIX}.gbwt -Z ${VG_PREFIX}.gbz
# vg convert -x --drop-haplotypes ${VG_PREFIX}.gbz > ${VG_PREFIX}.xg
# vg deconstruct ${VG_PREFIX}.xg -g ${VG_PREFIX}.gbwt -a -t 12 > ${VG_PREFIX}.2.vcf
# vg convert ${VG_PREFIX}.xg -p > ${VG_PREFIX}.pg
# vg snarls ${VG_PREFIX}.pg > ${VG_PREFIX}.snarls
# vg index -j ${VG_PREFIX}.dist -s ${VG_PREFIX}.snarls ${VG_PREFIX}.xg
# vg minimizer ${VG_PREFIX}.gbz -d ${VG_PREFIX}.dist -o ${VG_PREFIX}.min
################################################################################

# [1] Quality Control (COMMENTED - uncomment to enable trimming and QC)
################################################################################
# echo "Starting quality control for $SAMPLE_ID..."
# 
# # Download from SRA if needed
# # echo "~/sratoolkit.3.0.1-ubuntu64/bin/fasterq-dump -pe 24 $SAMPLE_ID"
# # cd /mnt/lustre/RDS-archive/downing/LSDV/FASTQ_FILES/
# # ~/sratoolkit.3.0.1-ubuntu64/bin/fasterq-dump -pe 24 $SAMPLE_ID
# # cd ..
# 
# if test -f /mnt/lustre/RDS-archive/downing/LSDV/FASTQ_FILES/$file1 ; then
#     # Paired-end processing
#     # fastqc /mnt/lustre/RDS-archive/downing/LSDV/FASTQ_FILES/$file1 &
#     # fastqc /mnt/lustre/RDS-archive/downing/LSDV/FASTQ_FILES/$file2 &
#     echo $file1 >> /mnt/lustre/RDS-archive/downing/LSDV/TRIM_FILES/trimming_stats.txt
#     echo $file2 >> /mnt/lustre/RDS-archive/downing/LSDV/TRIM_FILES/trimming_stats.txt
#     
#     # Quality trimming
#     # ~/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l 80 -i /mnt/lustre/RDS-archive/downing/LSDV/FASTQ_FILES/$file1 -o /mnt/lustre/RDS-archive/downing/LSDV/TRIM_FILES/$trim1 -v >> /mnt/lustre/RDS-archive/downing/LSDV/TRIM_FILES/trimming_stats.txt
#     # ~/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l 80 -i /mnt/lustre/RDS-archive/downing/LSDV/FASTQ_FILES/$file2 -o /mnt/lustre/RDS-archive/downing/LSDV/TRIM_FILES/$trim2 -v >> /mnt/lustre/RDS-archive/downing/LSDV/TRIM_FILES/trimming_stats.txt
#     
#     # Additional processing with fastp
#     # ~/FASTQC/fastp --overrepresentation_analysis --html /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastphtml1 -i /mnt/lustre/RDS-archive/downing/LSDV/TRIM_FILES/$trim1 -o /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1 &
#     # ~/FASTQC/fastp --overrepresentation_analysis --html /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastphtml2 -i /mnt/lustre/RDS-archive/downing/LSDV/TRIM_FILES/$trim2 -o /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2
#     # fastqc /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1 &
#     # fastqc /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2 &
# else
#     # Single-end processing
#     echo $file3 >> /mnt/lustre/RDS-archive/downing/LSDV/TRIM_FILES/trimming_stats.txt
#     # ~/FASTQC/fastq_quality_trimmer -Q 33 -t 30 -l 80 -i /mnt/lustre/RDS-archive/downing/LSDV/FASTQ_FILES/$file3 -o /mnt/lustre/RDS-archive/downing/LSDV/TRIM_FILES/$trim3 -v >> /mnt/lustre/RDS-archive/downing/LSDV/TRIM_FILES/trimming_stats.txt
#     # ~/FASTQC/fastp --overrepresentation_analysis --html /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastphtml3 -i /mnt/lustre/RDS-archive/downing/LSDV/TRIM_FILES/$trim3 -o /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3
#     # fastqc /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3 &
# fi
################################################################################

# [2] Map reads using vg giraffe
echo "Mapping reads for $SAMPLE_ID..."
gam="${SAMPLE_ID}.gam"
gaf="${SAMPLE_ID}.gaf"
pack="${SAMPLE_ID}.pack"
depth="${SAMPLE_ID}.depth.csv"
bam3="${SAMPLE_ID}.bam"

# Check for paired-end files first
if test -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1 ; then
    echo "  Paired-end mapping with vg giraffe to ${VG_PREFIX}..."
    echo "  vg giraffe -Z ${VG_PREFIX}.gbz -H ${VG_PREFIX}.gbwt -p -m ${VG_PREFIX}.min -d ${VG_PREFIX}.dist -t 22 -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1 -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2 > GAM_FILES/$gam"
    vg giraffe -Z ${VG_PREFIX}.gbz -H ${VG_PREFIX}.gbwt -p -m ${VG_PREFIX}.min -d ${VG_PREFIX}.dist -t 22 -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1 -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2 > GAM_FILES/$gam # 2> ERROR_FILES/${SAMPLE_ID}.gam2.errors.txt
    
    echo "  Generating BAM output..."
    if [ $REF_TYPE -eq 121 ]; then
        # For 121-genome reference, use additional -x parameter
        vg giraffe -x ${VG_PREFIX}.xg -Z ${VG_PREFIX}.gbz -H ${VG_PREFIX}.gbwt -p -t 4 -m ${VG_PREFIX}.min -d ${VG_PREFIX}.dist -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1 -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2 -o BAM > BAM_FILES/$bam3
    else
        vg giraffe -Z ${VG_PREFIX}.gbz -H ${VG_PREFIX}.gbwt -p -m ${VG_PREFIX}.min -d ${VG_PREFIX}.dist -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp1 -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp2 -o BAM > BAM_FILES/$bam3
    fi
    echo "  Done paired-end mapping"

elif test -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3 ; then
    echo "  Single-end mapping with vg giraffe to ${VG_PREFIX}..."
    echo "  vg giraffe -Z ${VG_PREFIX}.gbz -H ${VG_PREFIX}.gbwt -p -m ${VG_PREFIX}.min -d ${VG_PREFIX}.dist -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3 > GAM_FILES/$gam"
    vg giraffe -Z ${VG_PREFIX}.gbz -H ${VG_PREFIX}.gbwt -p -m ${VG_PREFIX}.min -d ${VG_PREFIX}.dist -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3 > GAM_FILES/$gam # 2> ERROR_FILES/${SAMPLE_ID}.gam2.errors.txt
    
    echo "  Generating BAM output..."
    if [ $REF_TYPE -eq 121 ]; then
        # For 121-genome reference, use additional -x parameter
        vg giraffe -x ${VG_PREFIX}.xg -Z ${VG_PREFIX}.gbz -p -t 22 -m ${VG_PREFIX}.min -d ${VG_PREFIX}.dist -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3 -o BAM > BAM_FILES/$bam3
    else
        vg giraffe -Z ${VG_PREFIX}.gbz -H ${VG_PREFIX}.gbwt -p -m ${VG_PREFIX}.min -d ${VG_PREFIX}.dist -f /mnt/lustre/RDS-archive/downing/LSDV/KRAKEN_VALID_FILES/$fastp3 -o BAM > BAM_FILES/$bam3
    fi
    echo "  Done single-end mapping"
else
    echo "  Error: No valid input files found for $SAMPLE_ID"
    exit 1
fi

# Generate mapping statistics
echo "Generating mapping statistics..."
gams="${SAMPLE_ID}.stat"
vg stats -a GAM_FILES/$gam > GAM_FILES/$gams

echo "Checking BAM file integrity..."
~/samtools-1.9/samtools quickcheck -v BAM_FILES/$bam3

# [3] BAM processing and quality assessment
echo "Processing BAM files..."
bam="${SAMPLE_ID}.bam"
sorted="${SAMPLE_ID}.sorted.bam"
rmdup="${SAMPLE_ID}.rmdup.sorted.bam"
cov="${SAMPLE_ID}.coverage.txt"
qual="${SAMPLE_ID}.html"

echo "Sorting BAM..."
~/samtools-1.9/samtools sort -n BAM_FILES/$bam3 -o BAM_FILES/$sorted
sorted2="${sorted}.2"
~/samtools-1.9/samtools fixmate -m BAM_FILES/$sorted BAM_FILES/$sorted2
~/samtools-1.9/samtools sort BAM_FILES/$sorted2 -o BAM_FILES/$sorted

echo "Removing duplicates..."
~/samtools-1.9/samtools markdup -r -s BAM_FILES/$sorted BAM_FILES/$rmdup

echo "Indexing BAM..."
~/samtools-1.9/samtools index BAM_FILES/$rmdup

echo "Generating BAM statistics..."
flg="${SAMPLE_ID}.flagstat"
~/samtools-1.9/samtools flagstat BAM_FILES/$rmdup > BAM_FILES/$flg

# Coverage analysis
/usr/bin/samtools coverage BAM_FILES/$rmdup > COVERAGE/$cov 2> ERROR_FILES/${SAMPLE_ID}.cov.errors.txt

# Quality assessment with Qualimap
cd QUALIMAP/
echo "Running Qualimap analysis..."
~/qualimap_v2.2.1/qualimap bamqc -bam ../BAM_FILES/$rmdup -outfile $qual -outdir $SAMPLE_ID -outformat html
cd ..

# [4] Variant calling with vg
echo "Calling variants with vg..."
vcf="${SAMPLE_ID}.vg.vcf"
vcfgz="${vcf}.gz"
aug="${SAMPLE_ID}${AUG_SUFFIX}"
snarls="${SAMPLE_ID}.snarls"

echo "Augmenting graph..."
vg augment ${VG_PREFIX}.pg GAM_FILES/$gam -m 3 -q 10 -Q 10 -A GAM_FILES2/$gam > AUG/$aug
ls -lt GAM_FILES/$gam GAM_FILES2/$gam AUG/$aug

echo "Computing snarls..."
vg snarls AUG/$aug > AUG/$snarls

echo "Packing reads..."
vg pack -x AUG/$aug -g GAM_FILES2/$gam -Q 10 -o PACK_FILES/$pack

echo "Calling variants..."
vg call -d 1 AUG/$aug -a -r AUG/$snarls -k PACK_FILES/$pack -s $SAMPLE_ID > VG_VCF_FILES/$vcf 2> ERROR_FILES/${SAMPLE_ID}.vg.errors.txt

echo "Computing depth..."
vg depth ${VG_PREFIX}.vg -m 6 -k PACK_FILES/$pack > DEPTH_FILES/$depth 2> ERROR_FILES/${SAMPLE_ID}.depth.errors.txt &

# Compress and index VCF
echo "Compressing and indexing VG VCF..."
rm -rf VG_VCF_FILES/$vcfgz
if [ $REF_TYPE -eq 121 ]; then
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/bgzip VG_VCF_FILES/$vcf
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/tabix -f -p vcf VG_VCF_FILES/$vcfgz
else
    bgzip VG_VCF_FILES/$vcf
    tabix -f -p vcf VG_VCF_FILES/$vcfgz
fi

# Normalize VG variants
echo "Normalizing VG variants..."
normvg="${SAMPLE_ID}.norm.vg.vcf.gz"
bcftools norm -c w -f $FASTA -m-both -Oz -o VG_VCF_FILES/$normvg VG_VCF_FILES/$vcfgz
if [ $REF_TYPE -eq 121 ]; then
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/tabix -f -p vcf VG_VCF_FILES/$normvg
else
    tabix -f -p vcf VG_VCF_FILES/$normvg
fi

# [5] Variant calling with FreeBayes
echo "Calling variants with FreeBayes..."
fb="${SAMPLE_ID}.fb.vcf"
fbgz="${fb}.gz"
rm -rf FB_VCF_FILES/$fbgz
~/freebayes/build/freebayes -f $FASTA -F 0.01 -g 2000 -p 1 --min-alternate-count 1 --min-alternate-fraction 0.001 BAM_FILES/$rmdup > FB_VCF_FILES/$fb 2> ERROR_FILES/${SAMPLE_ID}.fb.errors.txt

if [ $REF_TYPE -eq 121 ]; then
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/bgzip FB_VCF_FILES/$fb
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/tabix -f -p vcf FB_VCF_FILES/$fbgz
else
    bgzip FB_VCF_FILES/$fb
    tabix -f -p vcf FB_VCF_FILES/$fbgz
fi

# Normalize FreeBayes variants
echo "Normalizing FreeBayes variants..."
normfb="${SAMPLE_ID}.norm.fb.vcf.gz"
bcftools norm -c w -f $FASTA -m-both -Oz -o FB_VCF_FILES/$normfb FB_VCF_FILES/$fbgz
if [ $REF_TYPE -eq 121 ]; then
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/tabix -f -p vcf FB_VCF_FILES/$normfb
else
    tabix -f -p vcf FB_VCF_FILES/$normfb
fi

# [6] Variant calling with BCFtools
echo "Calling variants with BCFtools..."
bcf="${SAMPLE_ID}.bcf.vcf"
bcfgz="${bcf}.gz"
rm -rf BCF_VCF_FILES/$bcfgz
bcftools mpileup -A -Ob -f $FASTA BAM_FILES/$rmdup | bcftools call -cvO v --ploidy 1 -o BCF_VCF_FILES/$bcf

if [ $REF_TYPE -eq 121 ]; then
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/bgzip BCF_VCF_FILES/$bcf
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/tabix -f -p vcf BCF_VCF_FILES/$bcfgz
else
    bgzip BCF_VCF_FILES/$bcf
    tabix -f -p vcf BCF_VCF_FILES/$bcfgz
fi

# Normalize BCFtools variants
echo "Normalizing BCFtools variants..."
normbcf="${SAMPLE_ID}.norm.bcf.vcf.gz"
bcftools norm -c w -f $FASTA -m-both -Oz -o BCF_VCF_FILES/$normbcf BCF_VCF_FILES/$bcfgz
if [ $REF_TYPE -eq 121 ]; then
    /mnt/lustre/RDS-live/downing/miniconda3/envs/snippy/bin/tabix -f -p vcf BCF_VCF_FILES/$normbcf
else
    tabix -f -p vcf BCF_VCF_FILES/$normbcf
fi

# [7] Count variants and generate statistics
echo "Counting variants..."
echo "VG SNPs:"
gzip -dc VG_VCF_FILES/$normvg | grep -v -c "^#"
echo "FreeBayes SNPs:"
gzip -dc FB_VCF_FILES/$normfb | grep -v -c "^#"
echo "BCFtools SNPs:"
gzip -dc BCF_VCF_FILES/$normbcf | grep -v -c "^#"

# Generate transition/transversion statistics
echo $SAMPLE_ID >> TSTV.txt
bcftools stats VG_VCF_FILES/$normvg | grep TSTV >> TSTV.txt

# Log flagstat information for reference type 3
if [ $REF_TYPE -eq 3 ]; then
    echo $SAMPLE_ID >> flagstat.txt
fi

# [8] Generate consensus sequences
echo "Generating consensus sequences..."
fasta1="${SAMPLE_ID}.LA.fasta"
fasta2="${SAMPLE_ID}.LR.fasta"
bcftools consensus -H LA -f $FASTA BCF_VCF_FILES/$normbcf -o CONSENSUS/$fasta1 & # reference allele
bcftools consensus -H LR -f $FASTA BCF_VCF_FILES/$normbcf -o CONSENSUS/$fasta2 & # alternate allele

echo "Pipeline completed for sample $SAMPLE_ID with reference type $REF_TYPE"
echo "Results available in respective output directories:"
echo "  - Alignments: GAM_FILES/, BAM_FILES/"
echo "  - Variants: VG_VCF_FILES/, FB_VCF_FILES/, BCF_VCF_FILES/"
echo "  - Quality: QUALIMAP/, COVERAGE/"
echo "  - Consensus: CONSENSUS/"

#!/usr/bin/env perl
#
# PURPOSE: Parses VCF, BAM, and GAM files from a specified directory structure
#          to extract SNP counts and alignment metrics.
#
# USAGE:
# perl 02_BAM_GAM_metrics.pl --datadir /path/to/LSDV_data \
#                           --outdir ./metrics_output
#
# DESCRIPTION:
# This script is a rewritten version to be portable. It replaces hard-coded paths
# with command-line arguments. It expects the --datadir to contain subfolders
# for each mapping strategy (e.g., LSDV1, LSDVG, LSDVVG_ALL) and that these
# subfolders contain the VCF_FILES, BAM_FILES, and GAM_FILES directories.

use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use File::Path qw(make_path);

# --- Command Line Argument Parsing ---
my $data_dir;
my $out_dir = '.';
GetOptions(
    'datadir=s' => \$data_dir,
    'outdir=s'  => \$out_dir
) or die "Error in command line arguments\n";

die "Usage: $0 --datadir /path/to/data --outdir /path/to/output\n" unless $data_dir;

# Create output directory if it doesn't exist
make_path($out_dir) unless -d $out_dir;
print "Data directory: $data_dir\n";
print "Output directory: $out_dir\n";

# --- Helper Functions ---
# Safely read a list of samples from a file
sub read_sample_list {
    my ($file) = @_;
    my %samples;
    open(my $fh, '<', $file) or die "Cannot open sample list $file: $!";
    while (my $line = <$fh>) {
        chomp $line;
        $samples{$line} = 1;
    }
    close $fh;
    return \%samples;
}

# --- Load Sample Lists ---
# Note: These files must exist in the directory where the script is run.
my $wgs_samples        = read_sample_list('acc_list.wgs.txt');
my $metagenomic_samples = read_sample_list('acc_list.metagenomic.txt');
my $amplicon_samples   = read_sample_list('acc_list.amplicon.txt');
my $ont_samples        = read_sample_list('acc_list.ont.txt');

# --- Open Output Files ---
my %out_files;
my @categories = ('WGS', 'METAGENOMIC', 'AMPLICON', 'ONT');
my @stat_types = ('SNPS', 'BAMSTATS', 'GAMSTATS');

foreach my $cat (@categories) {
    foreach my $stat (@stat_types) {
        my $filename = File::Spec->catfile($out_dir, "$stat.$cat.txt");
        open(my $fh, '>', $filename) or die "Cannot open output file $filename: $!";
        $out_files{"${stat}_${cat}"} = $fh;
    }
}

# --- Process VCF Files for SNP Counts ---
print "Processing VCF files for SNP counts...\n";
my @vcf_files = glob("$data_dir/LSDV*/*_VCF_FILES/*norm*vcf*gz");

foreach my $file (@vcf_files) {
    # Skip excluded files
    next if $file =~ /SRR10394925t|255939|23347779|FINAL/;

    my @parts = File::Spec->splitdir($file);
    my $mapper_dir = $parts[-3]; # e.g., LSDVG_3
    my ($sample_name) = $parts[-1] =~ /^(.*?)\./;

    # Run command to count SNPs
    my $snp_count = `gzip -dc "$file" | grep -vc "^#"`;
    chomp $snp_count;

    # Determine sample type and write to correct file
    my $fh;
    if ($wgs_samples->{$sample_name}) {
        $fh = $out_files{'SNPS_WGS'};
    } elsif ($metagenomic_samples->{$sample_name}) {
        $fh = $out_files{'SNPS_METAGENOMIC'};
    } elsif ($amplicon_samples->{$sample_name}) {
        $fh = $out_files{'SNPS_AMPLICON'};
    } elsif ($ont_samples->{$sample_name}) {
        $fh = $out_files{'SNPS_ONT'};
    }

    if ($fh) {
        # Format: Sample  Mapper  Rate
        print $fh "$sample_name\t$mapper_dir\t$snp_count\n";
    }
}

# --- Process BAM Files for Mapping Stats ---
print "Processing BAM flagstat files...\n";
my @bam_files = glob("$data_dir/LSDV*/*_FILES*/*flagstat");

foreach my $file (@bam_files) {
    next if $file =~ /SRR10394925t|255939|23347779|FINAL/;

    my @parts = File::Spec->splitdir($file);
    my $mapper_dir = $parts[-3];
    my ($sample_name) = $parts[-1] =~ /^(.*?)\.flagstat/;

    my $mapped_line = `grep "mapped (" "$file"`;
    next unless $mapped_line;
    my ($rate) = $mapped_line =~ /\((\d+\.\d+)\% :/;
    
    my $fh;
    if ($wgs_samples->{$sample_name}) {
        $fh = $out_files{'BAMSTATS_WGS'};
    } elsif ($metagenomic_samples->{$sample_name}) {
        $fh = $out_files{'BAMSTATS_METAGENOMIC'};
    } elsif ($amplicon_samples->{$sample_name}) {
        $fh = $out_files{'BAMSTATS_AMPLICON'};
    } elsif ($ont_samples->{$sample_name}) {
        $fh = $out_files{'BAMSTATS_ONT'};
    }
    
    if ($fh) {
        print $fh "$sample_name\t$mapper_dir\t$rate\n";
    }
}

# --- Process GAM Files for Alignment Stats ---
print "Processing GAM stat files...\n";
my @gam_files = glob("$data_dir/LSDV*/*_FILES*/*stat");

foreach my $file (@gam_files) {
    next if $file =~ /\.flagstat$/; # Skip flagstat files from BAM processing
    next if $file =~ /SRR10394925t|255939|23347779/;
    
    my @parts = File::Spec->splitdir($file);
    my $mapper_dir = $parts[-3];
    my ($sample_name) = $parts[-1] =~ /^(.*?)\.stat/;

    my $aligned_line = `grep "aligned" "$file"`;
    my $primary_line = `grep "primary" "$file"`;
    next unless $aligned_line && $primary_line;

    my ($aligned_count) = $aligned_line =~ /(\d+)\s+aligned/;
    my ($primary_count) = $primary_line =~ /(\d+)\s+primary/;

    my $rate = 'NA';
    if ($primary_count > 0) {
        $rate = sprintf("%.2f", (100 * $aligned_count) / $primary_count);
    }

    my $fh;
    if ($wgs_samples->{$sample_name}) {
        $fh = $out_files{'GAMSTATS_WGS'};
    } elsif ($metagenomic_samples->{$sample_name}) {
        $fh = $out_files{'GAMSTATS_METAGENOMIC'};
    } elsif ($amplicon_samples->{$sample_name}) {
        $fh = $out_files{'GAMSTATS_AMPLICON'};
    } elsif ($ont_samples->{$sample_name}) {
        $fh = $out_files{'GAMSTATS_ONT'};
    }

    if ($fh) {
        print $fh "$sample_name\t$mapper_dir\t$rate\n";
    }
}

# --- Close All File Handles ---
foreach my $key (keys %out_files) {
    close($out_files{$key});
}

print "Script finished successfully. Metrics written to '$out_dir'.\n";

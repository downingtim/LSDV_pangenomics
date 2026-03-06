# Analysis Scripts

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Pipeline Status](https://img.shields.io/badge/Pipeline-Production-green.svg)](https://github.com/your-org/lsdv-pipeline)
[![Version](https://img.shields.io/badge/Version-1.0.0-blue.svg)](https://github.com/your-org/lsdv-pipeline/releases)

This directory contains the core scripts for the LSDV pangenome analysis. The LSDV pangenome variation graph (PVG) analysis is a set of workflows designed for analyzing Lumpy Skin Disease Virus (LSDV) sequencing libraries using PVG-based approaches. This repository contains the scripts to perform mapping, variant calling, and generate figures for publication.

The scripts are numbered to suggest a logical order of execution and have been revised from their original versions to be more portable, readable, and user-friendly by using command-line arguments instead of hard-coded paths and filenames.

**General Dependencies:** R (with packages optparse, dplyr, ggplot2, etc.), Perl (with Getopt::Long), Python 3, and command-line tools like bcftools and samtools.

## Installation

This pipeline requires Conda for managing dependencies.

1. Clone the repository (replace `your-username` with your GitHub username or organization):
   ```bash
   git clone https://github.com/your-username/lsdv-pangenomics.git
   cd lsdv-pangenomics
   ```

2. Create the Conda environment from the `environment.yml` file. This file lists all required software.
   ```bash
   conda env create -f environment.yml
   conda activate lsdv-pangenomics
   ```

## Usage

To run the pipeline, first obtain the necessary input data. An example is provided below.

1. **Download test data:**
   ```bash
   # Add your commands here to download example FASTQ files into a 'data/' directory.
   # For example:
   # wget -P data/ http://example.com/path/to/sample1.fastq.gz
   ```

2. **Run the full pipeline:**
   ```bash
   # This example assumes a master script named 'run_pipeline.sh'
   bash run_pipeline.sh
   ```
   Results will be generated in the `results/` directory, which is created automatically if it doesn't exist.

## Workflow Description

The pipeline consists of several numbered scripts executed in sequence:

1. **`01_qc.R`**: Performs quality control analysis on coverage data and generates QC plots
2. **`02_BAM_GAM_metrics.pl`**: Extracts SNP counts and alignment metrics from VCF, BAM, and GAM files
3. **`03_generate_summary_plots.R`**: Generates summary boxplots comparing variant callers and mapping strategies
4. **`04_*` **: <removed>
5. **`05_get_liftover.py`**: Extracts liftover failure statistics from log files
6. **`06_get_liftover.R`**: Visualizes liftover failure rates across samples
7. **`07_get_stats.pl`**: Calculates Transition/Transversion (Ti/Tv) ratios from VCF files
8. **`08_*`**: <removed>
9. **`09_plot_sites.screened.R`**: Visualizes genomic locations of filtered variants
10. **`10_finalise_snps.pl`**: Performs final SNP merging across callers and mappers
11. **`11_analyze_heterozygosity.R`**: Analyzes and visualizes allele frequencies from heterozygosity data
12. **`12_finalise_snps.R`**: Creates scatter plots comparing SNP counts between mapping methods
13. **`13_get.titv.ratios.R`**: Calculates and analyzes Ti/Tv ratios
14. **`14_get_snps.R`**: Aggregates statistics from multiple sources into comprehensive summary
15. **`15_plot_coverage.R`**: Investigates correlation between coverage and SNP detection
16. **`16_plot_snp_overlap.all.R`**: Visualizes SNP overlap between Minimap2 and PVG approaches
17. **`17_get.ratios.R`**: Calculates Ti/Tv ratios from screened and unscreened results
18. **`18_plot_snp_frequency_Figure2.R`**: Generates Figure 2 - mutation density across genome with CDS annotation
19. **`19_snp_variation_Figure5.R`**: Generates Figure 5 - multi-panel SNP distribution and method comparison
20. **`20_sim_reads_viz_Figure6.R`**: Create Figure 6 of the precision, recall and F1, along with other metrics
21. **`21_sim_libraries_viz_Figure_S16.R`**: Visualise read mapping SNP comparisons for Figure S16.

## Data Preparation & QC (Scripts 01-02)

### SCRIPT: 01_qc.R

**PURPOSE:** Reads raw coverage data from multiple mapping experiments, merges it with sample metadata, and generates quality control plots comparing depth, base quality (BQ), and mapping quality (MQ).

**USAGE:**
```bash
Rscript 01_qc.R --samples_list <path/to/samples.txt> --metadata <path/to/metadata.csv> --coverage_dir <path/to/parent_dir> --out_csv <output.csv> --out_prefix <plots_prefix>
```

**ARGUMENTS:**
- `--samples_list`: Path to a plain text file with one sample name per line.
- `--metadata`: Path to the tab-separated metadata file (cols: Sample, Library_type).
- `--coverage_dir`: Parent directory that contains the experiment subfolders (e.g., LSDV1, LSDVG).
- `--out_csv`: Path for the output summary CSV file.
- `--out_prefix`: Prefix for the output PDF plot files (e.g., qc_plots_Depth_BQ.pdf).

### SCRIPT: 02_BAM_GAM_metrics.pl

**PURPOSE:** Parses VCF, BAM, and GAM files from a specified directory structure to extract SNP counts and alignment metrics. This is a portable version of the original script.

**USAGE:**
```bash
perl 02_BAM_GAM_metrics.pl --datadir <path/to/LSDV_data> --outdir <path/to/output_dir>
```

**ARGUMENTS:**
- `--datadir`: The top-level directory containing mapper subfolders (e.g., LSDV1, LSDVG) which in turn contain VCF_FILES, BAM_FILES, etc.
- `--outdir`: Directory where the output metrics files (SNPS.*.txt, BAMSTATS.*.txt) will be saved.
- **Note:** Expects sample lists (e.g., acc_list.wgs.txt) to be in the current working directory.

## Core Analysis & Plotting (Scripts 03-18)

### SCRIPT: 03_generate_summary_plots.R (Consolidates 03, 04 and 08)

**PURPOSE:** Reads various statistics files (SNPS, GAMSTATS, BAMSTATS, TSTV) and generates numerous summary boxplots comparing performance across different variant callers and mapping strategies.

**USAGE:**
```bash
Rscript 03_generate_summary_plots.R --input_dir <path/to/metrics> --output_dir <path/to/plots>
```

**ARGUMENTS:**
- `--input_dir`: Directory containing the metrics files generated by script 02 (e.g., SNPS.ILLUMINA.txt, TSTV.WGS.txt).
- `--output_dir`: Directory where all output PDF plots will be saved.

### SCRIPT: 05_get_liftover.py

**PURPOSE:** Parses log files from a "liftover" process to extract the percentage of variants that failed to be transferred. Compiles this data into a single summary CSV.

**USAGE:**
```bash
python3 05_get_liftover.py --metadata <path/to/metadata.csv> --input_dir <path/to/liftover_parent_dir> --output_csv <output.csv>
```

**ARGUMENTS:**
- `--metadata`: Path to the metadata file.
- `--input_dir`: Parent directory containing the "LIFTOVER_*" subfolders.
- `--output_csv`: Path for the output summary CSV file.

### SCRIPT: 06_get_liftover.R

**PURPOSE:** Reads the summary CSV from script 05 and generates a multi-faceted plot visualising the fraction of variants that were not successfully lifted over.

**USAGE:**
```bash
Rscript 06_get_liftover.R --input_csv <path/to/get_liftover.csv> --output_pdf <output.pdf>
```

**ARGUMENTS:**
- `--input_csv`: Path to the summary CSV file generated by script 05.
- `--output_pdf`: Path for the output PDF plot.

### SCRIPT: 07_get_stats.pl

**PURPOSE:** Calculates Transition/Transversion (Ti/Tv) ratios from VCF files using 'bcftools stats'.

**USAGE:**
```bash
perl 07_get_stats.pl --input_dir <path/to/screened_vcfs> --outdir <path/to/output_dir>
```

**ARGUMENTS:**
- `--input_dir`: Directory containing the "SCREENED_MERGED_*" subfolders.
- `--outdir`: Directory to save the output TSTV.*.txt files.
- **Note:** Expects sample lists (e.g., acc_list.wgs.txt) to be in the current working directory.

### SCRIPT: 09_plot_sites.screened.R

**PURPOSE:** Visualises the genomic locations and frequency of variants that were screened out during filtering by plotting "hotspots" of filtered mutations.

**USAGE:**
```bash
Rscript 09_plot_sites.screened.R --metadata <path/to/metadata.csv> --input_dir <path/to/screened_data> --output_dir <path/to/plots>
```

**ARGUMENTS:**
- `--metadata`: Path to the metadata file.
- `--input_dir`: Parent directory containing the "SCREENED_MERGED_*" subfolders.
- `--output_dir`: Directory where output plots and summary CSV will be saved.

### SCRIPT: 10_finalise_snps.pl

**PURPOSE:** A core workflow script that performs the final merging of SNPs from different callers and mappers using 'bcftools'.

**USAGE:**
```bash
perl 10_finalise_snps.pl --metadata <path/to/metadata.csv> --input_dir <path/to/screened_data> --output_dir <path/to/final_snps>
```

**ARGUMENTS:**
- `--metadata`: Path to the metadata file.
- `--input_dir`: Parent directory containing the "SCREENED_MERGED_*" subfolders.
- `--output_dir`: Top-level directory where output subfolders (FINAL_SNPS, CONCAT, etc.) will be created.

### SCRIPT: 11_analyze_heterozygosity.R (Consolidates 11, 18)

**PURPOSE:** Analyses and visualises allele frequencies from heterozygosity data by processing .vcf.het files from multiple experiments.

**USAGE:**
```bash
Rscript 11_analyze_heterozygosity.R --metadata <path/to/metadata.csv> --input_dir <path/to/screened_data> --output_dir <path/to/het_plots>
```

**ARGUMENTS:**
- `--metadata`: Path to the metadata file.
- `--input_dir`: Parent directory containing the "SCREENED_MERGED_*" subfolders.
- `--output_dir`: Directory to save the output PDF plots.

### SCRIPT: 12_finalise_snps.R

**PURPOSE:** Reads SNP count summary files from script 10 and creates scatter plots comparing SNP counts between Minimap2 and various PVG-based methods.

**USAGE:**
```bash
Rscript 12_finalise_snps.R --input_dir <path/to/rates_files> --output_dir <path/to/plots>
```

**ARGUMENTS:**
- `--input_dir`: Directory containing the `rates.*.txt` files generated by script 10.
- `--output_dir`: Directory to save the output PDF plots.

### SCRIPT: 14_get_snps.R

**PURPOSE:** A large data aggregation script that combines statistics from many sources (logs, .sites files, VCFs) into a single comprehensive summary table.

**USAGE:**
```bash
Rscript 14_get_snps.R --datadir <path/to/project> --metadata <path/to/metadata.csv> --outfile <path/to/All_SNPs.csv>
```

**ARGUMENTS:**
- `--datadir`: The top-level project directory containing the UNSCREENED, SCREENED, and LIFTOVER folders.
- `--metadata`: Path to the metadata file.
- `--outfile`: Path for the final, comprehensive output CSV file.

### SCRIPT: 15_plot_coverage.R

**PURPOSE:** Investigates the correlation between sequencing coverage on specific reference contigs and the number of unique SNPs detected.

**USAGE:**
```bash
Rscript 15_plot_coverage.R --coverage_dir <path/to/coverage_folder> --rates_file <path/to/rates.txt> --metadata <path/to/metadata.csv> --output_pdf <output.pdf>
```

**ARGUMENTS:**
- `--coverage_dir`: The specific directory containing the .coverage.txt files (e.g., LSDVG_3/COVERAGE).
- `--rates_file`: Path to the SNP rates summary file (e.g., rates.new3.txt).
- `--metadata`: Path to the metadata file.
- `--output_pdf`: Path for the output multi-panel PDF plot.

### SCRIPT: 16_plot_snp_overlap.all.R

**PURPOSE:** Visualises the overlap of SNPs detected between Minimap2 and various pangenome (PVG) approaches using stacked bar plots and scatter plots.

**USAGE:**
```bash
Rscript 16_plot_snp_overlap.all.R --input_dir <path/to/rates_files> --output_dir <path/to/plots>
```

**ARGUMENTS:**
- `--input_dir`: Directory containing the `rates.*.txt` files from script 10.
- `--output_dir`: Directory to save the various output PDF plots.

### SCRIPT: 17_get.ratios.R

**PURPOSE:** Calculates Ti/Tv ratios by parsing both raw `.sites` files and caller-specific VCF files from screened and unscreened results folders.

**USAGE:**
```bash
Rscript 17_get.ratios.R --datadir <path/to/project> --metadata <path/to/metadata.csv> --output_dir <path/to/ratios_output>
```

**ARGUMENTS:**
- `--datadir`: The top-level project directory containing the UNSCREENED_MERGED_* and SCREENED_MERGED_* folders.
- `--metadata`: Path to the metadata file.
- `--output_dir`: Directory to save the output CSV and PDF files.

## Publication Figure Scripts

### SCRIPT: 18_plot_snp_frequency_Figure2.R

**PURPOSE:** Generates Figure 2 for publication, showing mutation density across the genome aligned with an annotation of coding sequences (CDS).

**USAGE:**
```bash
Rscript 18_plot_snp_frequency_Figure2.R --vcf <path/to/vcf> --genbank <path/to/genbank> --output_pdf <output.pdf>
```

**ARGUMENTS:**
- `--vcf`: Path to the input VCF file (e.g., vcf/gfavariants.vcf).
- `--genbank`: Path to the input GenBank annotation file (e.g., KX894508.gb).
- `--output_pdf`: Path for the final PDF figure.

### SCRIPT: 19_snp_variation_Figure5.R

**PURPOSE:** Generates the complex, multi-panel Figure 5 for publication, showing SNP distribution, method comparison, and zoomed-in views of specific genomic regions.

**USAGE:**
```bash
Rscript 19_snp_variation_Figure5.R --vcf_dir <path/to/final_vcfs> --unique_snp_dir <path/to/unique_snps> --output_pdf <output.pdf>
```

**ARGUMENTS:**
- `--vcf_dir`: Path to the directory containing final VCF files (e.g., FINAL_SNPS_2/).
- `--unique_snp_dir`: Path to the directory containing the "Unique.*.txt" files.
- `--output_pdf`: Path for the final composite PDF figure.

### SCRIPT: 20_sim_reads_viz_Figure6.R

**PURPOSE:** Reads SNP calls from simulated read datasets generated across a range of sequencing depths and compares variant detection performance between linear reference mapping and pangenome graph mapping approaches. Calculates true positives, false positives, false negatives, precision, recall, and F1 score relative to known ground-truth SNP sets. Generates plots showing how SNP detection accuracy changes with read depth and produces the visualisations used in Figure 6.

**USAGE:**
```bash
Rscript 20_sim_reads_viz_Figure6.R
```

**INPUTS:**
- Simulated VCF files produced by different mapping strategies (e.g. LSDV1, LSDVG, LSDVG_3, LSDVG_6) located within the expected directory structure.
- Ground-truth SNP VCF files corresponding to the simulated genomes.

**OUTPUTS:**
- sim_comparison_results.csv – summary table of SNP detection metrics across depths and mapping methods.
- Multiple PDF plots showing detection rates, false positive and false negative fractions, precision, recall, and F1 score as a function of read depth.
- Heatmaps showing which SNPs are detected by linear versus pangenome mapping approaches across genomic positions.


### SCRIPT: 21_sim_libraries_viz_Figure_S16.R

**PURPOSE:** Evaluates SNP detection performance across simulated sequencing libraries generated from different reference genomes and coverage depths. Compares linear mapping with several pangenome graph configurations by measuring overlap with known truth SNP sets and calculating precision, recall, and F1 metrics. Produces visualisations summarising detection rates and SNP recovery patterns used in Supplementary Figure S16.

**USAGE:**
```bash
Rscript 21_sim_libraries_viz_Figure_S16.R
```

**INPUTS:**
- Variant call files (VCF) generated by different mapping strategies (e.g. LSDV1, LSDVG, LSDVG_3, LSDVG_6).
- Ground-truth SNP VCF files corresponding to the simulated genomes.

**OUTPUTS:**
- sim_comparison_results.csv – summary table of SNP detection metrics across depths, references, and mapping strategies.
- PDF plots illustrating SNP detection rates, precision–recall behaviour, and heatmaps showing SNP detection patterns across genomic positions.


## License

This project is licensed under the MIT License - see the `LICENSE` file for details.
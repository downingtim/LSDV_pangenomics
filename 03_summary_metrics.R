
The median rate of reads mapping was higher for the three-sample PVG for Giraffe (98.6%) compared to Minimap2 (92.2%) with VG-MAP much lower, and the same for the six- (98.6% vs 92.5%) and all-sample (98.6% vs 91.6%) PVGs.

The rate of read mapping using Giraffe was much higher than VG-MAP for Illumina data. Giraffe mapping for amplicon libraries was effectively the same across PVG options (medians 99.7-99.8%) and was much higher than the VG-MAP rates (53.2-54.8%). This was also true for WGS datasets (Giraffe 76.0-78.% vs VG-MAP 49.5-49.1%). For ONT datasets, performances were about the same with VG-MAP slightly better (medians 99.8-99.9%). 

_____

amplicon TS/TV
The average TS/TV rates for amplicon and WGS libraries were higher for mutations detected with BCFtools (2.18 and 2.42) ahead of VG (0.69 and 1.71) and Freebayes (0.32 and 0.29). Similarly, TS/TV rates were also higher for mutations detected from mapping with Minimap2 (1.49 and 1.71) more so than Freebayes (1.02 and 1.77) and VG (0.72 and 0.99).
For ONT data, Freebayes had the highest TS/TV rate (3.69) above VG (2.73) and VCFtools (0.50). Across mapping tools, Giraffe (2.62) and Minimap2 (2.61) had rates above VG-MAP (1.94).
This provided a rationale 
- can concat Giraffe and VG-MAP within BCFtools, within VG, within BCF
- bcftools concat -a -D vcfs -O v | bacftools norm -d all > out.vcf
- 

Amplicon:   Giraffe1_BCFtools, Minimap2_BCFtools, VG-MAP[1|3|6]_BCFtools 
[1] bcftools concat Giraffe1_BCFtools, VG-MAP[1|3|6]_BCFtools 
[2] bcftools merge [1] with Minimap2_BCFtools

# Resolve duplicate sample names, allowing multiallelic records

WGS/metagenomic: use BCFtools + Giraffe1, Minimap2 any, VG1/3/6; OR VG + VG-MAP any
as above with:
Any 2 of:  Giraffe1_BCFtools, Minimap2_BCFtools, VG-MAP[1|3|6]_BCFtools, VG-MAP_VG
SCREENED_MERGED_1/VG_BOKKAS_S10.vcf.gz
SCREENED_MERGED_3/VG_BOKKAS_S10.vcf.gz
SCREENED_MERGED_6/VG_BOKKAS_S10.vcf.gz
SCREENED_MERGED_ALL/VG_BOKKAS_S10.vcf.gz

ONT: use FB + any; VG + VG-MAP any (small sample size) 
Any 2 of:  Giraffe[any]_FB, Minimap2_FB, VG-MAP[any]_FB, VG-MAP_VG

__

library(dplyr) 
library(tidyr)
library(ggplot2)
library(ggpubr)

# Files to examine:

# SNPS: SNPS.AMPLICON.SE.txt SNPS.AMPLICON.PE.txt SNPS.METAGENOMIC.txt
#         SNPS.WGS.txt SNPS.ONT.txt
#        Format: ID, SNP_calller, Source, Rate

# BAMSTATS: BAMSTATS.AMPLICON.SE.txt BAMSTATS.AMPLICON.PE.txt
#           BAMSTATS.METAGENOMIC.txt BAMSTATS.WGS.txt BAMSTATS.ONT.txt
#        Format: ID, Source, Rate

# GAMSTATS: GAMSTATS.AMPLICON.SE.txt GAMSTATS.AMPLICON.PE.txt
#           GAMSTATS.METAGENOMIC.txt GAMSTATS.WGS.txt GAMSTATS.ONT.txt
#        Format: ID, Source, Rate

# Groups 1: amplicon_SE, amplicon_PE, metagenomic, wgs, ont
# Groups 2: illumina, ont
#

#######################################
# TSTV Rates

make_plots <- function(dataset, name1, label1, scale1){

dataset$Source[dataset$Source == 'LSDVVG'] <- 'VG-MAP_1'
dataset$Source[dataset$Source == 'LSDVVG_3'] <- 'VG-MAP_3'
dataset$Source[dataset$Source == 'LSDVVG_6'] <- 'VG-MAP_6'
dataset$Source[dataset$Source == 'LSDVVG_ALL'] <- 'VG-MAP_All'
dataset$Source[dataset$Source == 'LSDVG'] <- 'Giraffe_1'
dataset$Source[dataset$Source == 'LSDVG_3'] <- 'Giraffe_3'
dataset$Source[dataset$Source == 'LSDVG_6'] <- 'Giraffe_6'
dataset$Source[dataset$Source == 'LSDVG_ALL'] <- 'Giraffe_All'
dataset$Source[dataset$Source == 'LSDV1'] <- 'Minimap2_1'
  print(dataset %>% group_by(Source) %>% summarise(Median = median(na.omit(Rate))))

  if(scale1==1){ dataset$Rate = log10(dataset$Rate) } # for SNPs

  nameSummary <- paste("summary_", name1, ".pdf", sep="")
  nameSamples <- paste("samples_", name1, ".pdf", sep="")

  j <- ggplot(data = dataset, aes(x=Caller, y=Rate, color=Caller)) +
    geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1,alpha=0.4)+
    theme(axis.text.x=element_text(size=10, angle=90)) +ylab(label1)+ylim(0,8)
  pdf(nameSummary, width=4, height=5)
  print( j )
  dev.off()

  p <- ggplot(data = dataset, aes(x=Caller, y=Rate, color=Caller)) +
    geom_point( ) + theme(axis.text.x=element_text(size=10, angle=90))+
    ylab(label1) +ylim(0,8)
  pdf(nameSamples, width=5 + dim(dataset)[1]/99, height=5 + dim(dataset)[1]/111)
  print( p + facet_wrap(~Sample, ncol = 8) )
  dev.off()

  nameSummary2 <- paste("summary_source", name1, ".pdf", sep="")
  nameSamples2 <- paste("samples_source", name1, ".pdf", sep="")

  pdf(nameSummary2, width=5, height=5)
  print( ggplot(data = dataset, aes(x=Source, y=Rate, color=Source)) +
    geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1,alpha=0.4)+
    theme(axis.text.x=element_text(size=10, angle=90)) +ylab(label1)+ylim(0,8) )
  dev.off()

  p <- ggplot(data = dataset, aes(x=Source, y=Rate, color=Source)) +
    geom_point( ) + theme(axis.text.x=element_text(size=10, angle=90))+
    ylab(label1) +ylim(0,8)
#  pdf(nameSamples2, width=5 + dim(dataset)[1]/99, height=5 + dim(dataset)[1]/111)
#  print( p + facet_wrap(~Sample, ncol = 8) )
#  dev.off()

  nameSummary3 <- paste("summary_BY1_source", name1, ".pdf", sep="")
  nameSamples3 <- paste("samples_BY1_source", name1, ".pdf", sep="")

  pdf(nameSummary3, width=9, height=4)
  x2 <- ggplot(data = dataset, aes(x=Source, y=Rate, color=Source)) +
    geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1,alpha=0.4)+
    theme(axis.text.x=element_text(size=10, angle=90)) +ylab(label1)+ylim(0,8)
  print( x2 + facet_wrap(~Caller, ncol = 3) )
  dev.off()

  p2 <- ggplot(data = dataset, aes(x=Source, y=Rate, color=Source)) +
    geom_point( ) + theme(axis.text.x=element_text(size=10, angle=90))+
    ylab(label1) +ylim(0,8)
  p3 <- p2 + facet_wrap(~Caller, ncol = 3)
#  pdf(nameSamples3, width=15 + dim(dataset)[1]/99, height=15 + dim(dataset)[1]/111)
 # print( p3 + facet_wrap(~Sample, ncol = 8) )
 # dev.off()

  nameSummary4 <- paste("summary_BY2_source", name1, ".pdf", sep="")
  nameSamples4 <- paste("samples_BY2_source", name1, ".pdf", sep="")

  pdf(nameSummary4, width=5.5, height=4)
  x3 <- ggplot(data = dataset, aes(x=Caller, y=Rate, color=Caller)) +
    geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1,alpha=0.4)+
    theme(axis.text.x=element_text(size=10, angle=90)) +ylab(label1)
  print( x3 + facet_wrap(~Source, ncol = 4) )
  dev.off()

  q2 <- ggplot(data = dataset, aes(x=Caller, y=Rate, color=Caller)) +
    geom_point( ) + theme(axis.text.x=element_text(size=10, angle=90))+
    ylab(label1)
  q3 <- q2 + facet_wrap(~Source, ncol = 3)
 # pdf(nameSamples4, width=14 + dim(dataset)[1]/111, height=14 + dim(dataset)[1]/122)
#  print( q3 + facet_wrap(~Sample, ncol = 8) )
#  dev.off()
}

#######################################
# SNPs
#######################################

# Just do summary comparison of Illumina vs ONT

snps_illumina <- data.frame(read.csv("SNPS.ILLUMINA.txt", sep="\t"))
colnames(snps_illumina) <- c("Sample", "Caller", "Source", "Rate")
str(snps_illumina) # Sample, Caller, Source, Rate
snps_illumina <- subset(snps_illumina, Rate>0)
# snps_illumina <- subset(snps_illumina, Caller !="FB2" & Caller !="BCF2")
median(na.omit(subset(snps_illumina, Caller=="BCF" )$Rate)) # 60
median(na.omit(subset(snps_illumina, Caller=="FB"  )$Rate)) # 11827
median(na.omit(subset(snps_illumina, Caller=="VG")$Rate)) # 2843.5

name1 <- "snps_illumina"
label1 <- "Log10 of numbers of SNPs found - Illumina"
dataset <- snps_illumina
pdf(paste("summary_", name1, ".pdf", sep=""), width=5, height=5)
print( ggplot(data = dataset, aes(x=Caller, y=log10(Rate), color=Caller)) +
    geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1,alpha=0.4)+
    theme(axis.text.x=element_text(size=10, angle=90))+ylab(label1)+ylim(1,6))
dev.off()

dataset$Source[dataset$Source == 'LSDV1'] <- 'LSDV_1_map'
dataset$Source[dataset$Source == 'LSDV6'] <- 'LSDV_6_map'
dataset$Source[dataset$Source == 'LSDV6_GBWT'] <- 'LSDV_6_giraffe'
dataset$Source[dataset$Source == 'LSDV_GBWT'] <- 'LSDV_1_giraffe'
dataset$Source[dataset$Source == 'LSDV_FINAL_GBWT'] <- 'LSDV_all_giraffe'
dataset$Source[dataset$Source == 'LSDVVG_ALL'] <- 'LSDV_all_map'
name1 <- "snps_illumina_source"
pdf(paste("summary_", name1, ".pdf", sep=""), width=5, height=5)
print( ggplot(data = dataset, aes(x=Source, y=log10(Rate), color=Source)) +
    geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1,alpha=0.4)+
    theme(axis.text.x=element_text(size=10, angle=90))+ylab(label1)+ylim(1,5))
dev.off()

print(dataset %>% group_by(Source) %>% summarise(Median = median(na.omit(Rate))))

make_plots(snps_illumina, "snps_illumina", "Number of SNPs - Illumina", 1)

# ONT

snps_ont <- data.frame(read.csv("SNPS.ONT.txt", sep="\t"))
colnames(snps_ont) <- c("Sample", "Caller", "Source", "Rate")
str(snps_ont) # Sample, Caller, Source, Rate
snps_ont <- subset(snps_ont, Rate>0)
snps_ont <- subset(snps_ont, Caller !="FB2" & Caller !="BCF2")
median(na.omit(subset(snps_ont, Caller=="BCF"  )$Rate)) # 20758
median(na.omit(subset(snps_ont, Caller=="FB"  )$Rate)) # 76379
median(na.omit(subset(snps_ont, Caller=="VG")$Rate)) # 3,622

name1 <- "snps_ont"
label1 <- "Log10 of numbers of SNPs found - ONT"
dataset <- snps_ont
dataset$Source[dataset$Source == 'LSDV1'] <- 'LSDV_1_map'
dataset$Source[dataset$Source == 'LSDV6'] <- 'LSDV_6_map'
dataset$Source[dataset$Source == 'LSDV6_GBWT'] <- 'LSDV_6_giraffe'
dataset$Source[dataset$Source == 'LSDV_GBWT'] <- 'LSDV_1_giraffe'
dataset$Source[dataset$Source == 'LSDV_FINAL_GBWT'] <- 'LSDV_all_giraffe'
dataset$Source[dataset$Source == 'LSDV_FINAL'] <- 'LSDV_all_map'
pdf(paste("summary_", name1, ".pdf", sep=""), width=4, height=5)
print( ggplot(data = dataset, aes(x=Caller, y=log10(Rate), color=Caller)) +
    geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1,alpha=0.4)+
    theme(axis.text.x=element_text(size=10, angle=90))+ylab(label1)+ylim(1,6))
dev.off()

name1 <- "snps_ont_source"
pdf(paste("summary_", name1, ".pdf", sep=""), width=5, height=5)
print( ggplot(data = dataset, aes(x=Source, y=log10(Rate), color=Source)) +
    geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1,alpha=0.4)+
    theme(axis.text.x=element_text(size=10, angle=90))+ylab(label1)+ylim(1,6))
dev.off()

print(dataset %>% group_by(Source) %>% summarise(Median = median(na.omit(Rate))))

make_plots(snps_ont, "snps_ont", "Number of SNPs - ONT", 1)


#######################################
# GAMSTATS
#######################################

make_gam <- function(dataset, name1, label1){

  nameSummary <- paste("summary_", name1, ".pdf", sep="")
 nameSamples <- paste("samples_", name1, ".pdf", sep="")
dataset$Rate[dataset$Rate == ''] <- 0 # reassign empty as zero

dataset$Source[dataset$Source == 'LSDVVG'] <- 'VG-MAP_1'
dataset$Source[dataset$Source == 'LSDVVG_3'] <- 'VG-MAP_3'
dataset$Source[dataset$Source == 'LSDVVG_6'] <- 'VG-MAP_6'
dataset$Source[dataset$Source == 'LSDVVG_ALL'] <- 'VG-MAP_All'
dataset$Source[dataset$Source == 'LSDVG'] <- 'Giraffe_1'
dataset$Source[dataset$Source == 'LSDVG_3'] <- 'Giraffe_3'
dataset$Source[dataset$Source == 'LSDVG_6'] <- 'Giraffe_6'
dataset$Source[dataset$Source == 'LSDVG_ALL'] <- 'Giraffe_All'
dataset$Source[dataset$Source == 'LSDV1'] <- 'Minimap2'
dataset <- subset(dataset, Source !="LSDV1")
dataset <- subset(dataset, Source !="Minimap2" )

print(dataset %>% group_by(Source) %>% summarise(Median = median(na.omit(Rate))))

pdf(nameSummary, width=4, height=4)
print( ggplot(dataset, aes(x = Source, y = Rate, color=Source)) +
  geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1, alpha=0.4)+
  theme(axis.text.x=element_text(size=10, angle=90))+
  ylab(label1) + ylim(0, 100) )
dev.off()

pdf(nameSamples, width=8 +dim(dataset)[1]/99, height=8 +dim(dataset)[1]/111)
p <- ggplot(dataset, aes(x = Source, y = Rate, color=Source))  +
  geom_point( ) + theme(axis.text.x=element_text(size=10, angle=90))+
  ylab(label1) + ylim(0, 100)
print( p + facet_wrap(~Sample, ncol = 6) )
dev.off()
}

gam_amplicon <- data.frame(read.csv("GAMSTATS.AMPLICON.txt", sep="\t"))
colnames(gam_amplicon) <- c("Sample", "Source", "Rate") 
make_gam(gam_amplicon, "gam_amplicon", "GAM - % of reads mapped - Amplicon")

gam_wgs <- data.frame(read.csv("GAMSTATS.WGS.txt", sep="\t"))
colnames(gam_wgs) <- c("Sample", "Source", "Rate") 
make_gam(gam_wgs, "gam_wgs", "GAM - % of reads mapped - WGS")

gam_metagenomic <- data.frame(read.csv("GAMSTATS.METAGENOMIC.txt", sep="\t"))
colnames(gam_metagenomic) <- c("Sample", "Source", "Rate") 
make_gam(gam_metagenomic, "gam_metagenomic", "GAM - % of reads mapped - Metagenomic")

gam_illumina <- rbind(gam_amplicon,gam_wgs,gam_metagenomic)

gam_illumina$Source[gam_illumina$Source == 'LSDVVG'] <- 'VG-MAP_1'
gam_illumina$Source[gam_illumina$Source == 'LSDVVG_3'] <- 'VG-MAP_3'
gam_illumina$Source[gam_illumina$Source == 'LSDVVG_6'] <- 'VG-MAP_6'
gam_illumina$Source[gam_illumina$Source == 'LSDVVG_ALL'] <- 'VG-MAP_All'
gam_illumina$Source[gam_illumina$Source == 'LSDVG'] <- 'Giraffe_1'
gam_illumina$Source[gam_illumina$Source == 'LSDVG_3'] <- 'Giraffe_3'
gam_illumina$Source[gam_illumina$Source == 'LSDVG_6'] <- 'Giraffe_6'
gam_illumina$Source[gam_illumina$Source == 'LSDVG_ALL'] <- 'Giraffe_All'

gam_illumina <- subset(gam_illumina, Source !="LSDV1" )
gam_illumina <- subset(gam_illumina, Source !="LSDV1_6" )
gam_illumina <- subset(gam_illumina, Source !="LSDV1_3" )
gam_illumina <- subset(gam_illumina, Source !="LSDV1_ALL" )

gam_ont <- data.frame(read.csv("GAMSTATS.ONT.txt", sep="\t"))
colnames(gam_ont) <- c("Sample", "Source", "Rate")
gam_ont2 <- subset(gam_ont, Rate>=40)
make_gam(gam_ont, "gam_ont", "GAM - % of reads mapped - ONT")
dataset <- subset(dataset, Source !="LSDV" )
gam_ont$Source[gam_ont$Source == 'LSDVVG'] <- 'VG-MAP_1'
gam_ont$Source[gam_ont$Source == 'LSDVVG_3'] <- 'VG-MAP_3'
gam_ont$Source[gam_ont$Source == 'LSDVVG_6'] <- 'VG-MAP_6'
gam_ont$Source[gam_ont$Source == 'LSDVVG_ALL'] <- 'VG-MAP_All'
gam_ont$Source[gam_ont$Source == 'LSDVG'] <- 'Giraffe_1'
gam_ont$Source[gam_ont$Source == 'LSDVG_3'] <- 'Giraffe_3'
gam_ont$Source[gam_ont$Source == 'LSDVG_6'] <- 'Giraffe_6'
gam_ont$Source[gam_ont$Source == 'LSDVG_ALL'] <- 'Giraffe_All'
gam_ont <- subset(gam_ont, Source !="LSDV1" )
gam_ont <- subset(gam_ont, Source !="LSDV1_6" )
gam_ont <- subset(gam_ont, Source !="LSDV1_3" )
gam_ont <- subset(gam_ont, Source !="LSDV1_ALL" )

pdf("GAM_ILLUMINA.pdf", width=5, height=6)
print( ggplot(gam_illumina, aes(x =Source, y = Rate, color=Source)) +
  geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1, alpha=0.4)+
  theme(axis.text.x=element_text(size=10, angle=90))+
  ylab("GAM - % of reads mapped - Illumina") + ylim(0, 100) )
dev.off()

gam_ont2 <- subset(gam_ont2, Source !="LSDV1" )
pdf("GAM_ONT.pdf", width=5, height=6)
print( ggplot(gam_ont2, aes(x =Source, y = Rate, color=Source)) +
  geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1, alpha=0.4)+
  theme(axis.text.x=element_text(size=10, angle=90))+
  ylab("GAM - % of reads mapped - ONT")   )
dev.off()

#######################################
# BAMSTATS
#######################################

make_bam <- function(dataset, name1, label1){

  nameSummary <- paste("summary_", name1, ".pdf", sep="")
 nameSamples <- paste("samples_", name1, ".pdf", sep="")
dataset$Rate[dataset$Rate == ''] <- 0 # reassign empty as zero
dataset$Source[dataset$Source == 'LSDVVG'] <- 'VG-MAP_1'
dataset$Source[dataset$Source == 'LSDVVG_3'] <- 'VG-MAP_3'
dataset$Source[dataset$Source == 'LSDVVG_6'] <- 'VG-MAP_6'
dataset$Source[dataset$Source == 'LSDVVG_ALL'] <- 'VG-MAP_All'
dataset$Source[dataset$Source == 'LSDVG'] <- 'Giraffe_1'
dataset$Source[dataset$Source == 'LSDVG_3'] <- 'Giraffe_3'
dataset$Source[dataset$Source == 'LSDVG_6'] <- 'Giraffe_6'
dataset$Source[dataset$Source == 'LSDVG_ALL'] <- 'Giraffe_All'
dataset$Source[dataset$Source == 'LSDV1'] <- 'Minimap2_1'
dataset <- subset(dataset, Source !="Minimap2" )
dataset <- subset(dataset, Source !="LSDV1_3" )
dataset <- subset(dataset, Source !="LSDV1_6" )
dataset <- subset(dataset, Source !="LSDV1_ALL" )

print(dataset %>% group_by(Source) %>% summarise(Median = median(na.omit(Rate))))

pdf(nameSummary, width=4, height=4)
print( ggplot(dataset, aes(x = Source, y = Rate, color=Source)) +
  geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1, alpha=0.4)+
  theme(axis.text.x=element_text(size=10, angle=90))+
  ylab(label1)  )
dev.off()

pdf(nameSamples, width=8 +dim(dataset)[1]/99, height=5 +dim(dataset)[1]/99)
p <- ggplot(dataset, aes(x = Source, y = Rate, color=Source))  +
  geom_point( ) + theme(axis.text.x=element_text(size=10, angle=90))+
  ylab(label1)
print( p + facet_wrap(~Sample, ncol = 6) )
dev.off()
}

bam_amplicon <- data.frame(read.csv("BAMSTATS.AMPLICON.txt", sep="\t"))
colnames(bam_amplicon) <- c("Sample", "Source", "Rate") 
make_bam(bam_amplicon, "bam_amplicon", "BAM - % of reads mapped - Amplicon")

bam_wgs <- data.frame(read.csv("BAMSTATS.WGS.txt", sep="\t"))
colnames(bam_wgs) <- c("Sample", "Source", "Rate") 
make_bam(bam_wgs, "bam_wgs", "BAM - % of reads mapped - WGS")

bam_metagenomic <- data.frame(read.csv("BAMSTATS.METAGENOMIC.txt", sep="\t"))
colnames(bam_metagenomic) <- c("Sample", "Source", "Rate") 
make_bam(bam_metagenomic, "bam_metagenomic", "BAM - % of reads mapped - Metagenomic")

bam_ont <- data.frame(read.csv("BAMSTATS.ONT.txt", sep="\t"))
colnames(bam_ont) <- c("Sample", "Source", "Rate")
make_bam(bam_ont, "bam_ont", "BAM - % of reads mapped - ONT")

bam_illumina <- rbind(bam_amplicon,bam_wgs,bam_metagenomic)

bam_illumina$Source[bam_illumina$Source == 'LSDVVG'] <- 'VG-MAP_1'
bam_illumina$Source[bam_illumina$Source == 'LSDVVG_3'] <- 'VG-MAP_3'
bam_illumina$Source[bam_illumina$Source == 'LSDVVG_6'] <- 'VG-MAP_6'
bam_illumina$Source[bam_illumina$Source == 'LSDVVG_ALL'] <- 'VG-MAP_All'
bam_illumina$Source[bam_illumina$Source == 'LSDVG'] <- 'Giraffe_1'
bam_illumina$Source[bam_illumina$Source == 'LSDVG_3'] <- 'Giraffe_3'
bam_illumina$Source[bam_illumina$Source == 'LSDVG_6'] <- 'Giraffe_6'
bam_illumina$Source[bam_illumina$Source == 'LSDVG_ALL'] <- 'Giraffe_All'
bam_illumina$Source[bam_illumina$Source == 'LSDV1'] <- 'Minimap2_1'
bam_illumina$Source[bam_illumina$Source == 'LSDV1_3'] <- 'Minimap2_3'
bam_illumina$Source[bam_illumina$Source == 'LSDV1_6'] <- 'Minimap2_6'
bam_illumina$Source[bam_illumina$Source == 'LSDV1_ALL'] <- 'Minimap2_All'
bam_illumina <- subset(bam_illumina, Source !="Minimap2" )
bam_illumina <- subset(bam_illumina, Source !="LSDV" )
bam_illumina <- subset(bam_illumina, Source !="Minimap2_3" )
bam_illumina <- subset(bam_illumina, Source !="Minimap2_6" )
bam_illumina <- subset(bam_illumina, Source !="Minimap2_All" )

bam_ont$Source[bam_ont$Source == 'LSDVVG'] <- 'VG-MAP_1'
bam_ont$Source[bam_ont$Source == 'LSDVVG_3'] <- 'VG-MAP_3'
bam_ont$Source[bam_ont$Source == 'LSDVVG_6'] <- 'VG-MAP_6'
bam_ont$Source[bam_ont$Source == 'LSDVVG_ALL'] <- 'VG-MAP_All'
bam_ont$Source[bam_ont$Source == 'LSDVG'] <- 'Giraffe_1'
bam_ont$Source[bam_ont$Source == 'LSDVG_3'] <- 'Giraffe_3'
bam_ont$Source[bam_ont$Source == 'LSDVG_6'] <- 'Giraffe_6'
bam_ont$Source[bam_ont$Source == 'LSDVG_ALL'] <- 'Giraffe_All'
bam_ont$Source[bam_ont$Source == 'LSDV1'] <- 'Minimap2_1'
bam_ont$Source[bam_ont$Source == 'LSDV1_3'] <- 'Minimap2_3'
bam_ont$Source[bam_ont$Source == 'LSDV1_6'] <- 'Minimap2_6'
bam_ont$Source[bam_ont$Source == 'LSDV1_ALL'] <- 'Minimap2_All'
bam_ont <- subset(bam_ont, Source !="Minimap2" )
bam_ont <- subset(bam_ont, Source !="LSDV" )
bam_ont <- subset(bam_ont, Source !="Minimap2_3" )
bam_ont <- subset(bam_ont, Source !="Minimap2_6" )
bam_ont <- subset(bam_ont, Source !="Minimap2_All" )

bam_illumina %>% group_by(Source) %>% summarise(Median = median(na.omit(Rate)))
bam_illumina  %>% group_by(Source) %>% summarise(mean(na.omit(Rate)))

pdf("BAM_ILLUMINA_giraffe.pdf", width=5, height=6)
print( ggplot(bam_illumina, aes(x =Source, y = Rate, color=Source)) +
  geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1, alpha=0.4)+
  theme(axis.text.x=element_text(size=10, angle=90))+
  ylab("BAM - % of reads mapped - Illumina")  )
dev.off()

pdf("BAM_ONT.pdf", width=5, height=6)
print( ggplot(bam_ont, aes(x =Source, y = Rate, color=Source)) +
  geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1, alpha=0.4)+
  theme(axis.text.x=element_text(size=10, angle=90))+
  ylab("BAM - % of reads mapped - ONT") )
dev.off()

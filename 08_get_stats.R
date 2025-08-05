#######################################################
#install.packages("dplyr") # Let's install  dplyr
library(dplyr)
#install.packages("tidyr") # Let's install tidyr to make our data tidy
library(tidyr)
library(ggplot2)
library(ggpubr)

####################################### TSTV Rates

make_plots <- function(dataset, name1, label1){ 
  nameSummary <- paste("summary_", name1, ".pdf", sep="")
  nameSamples <- paste("samples_", name1, ".pdf", sep="")

dataset$Source[dataset$Source == '1'] <- 'VG-MAP_1'
dataset$Source[dataset$Source == '3'] <- 'VG-MAP_3'
dataset$Source[dataset$Source == '6'] <- 'VG-MAP_6'
dataset$Source[dataset$Source == 'ALL'] <- 'VG-MAP_All'
dataset$Source[dataset$Source == '1_GBWT'] <- 'Giraffe_1'
dataset$Source[dataset$Source == '3_GBWT'] <- 'Giraffe_3'
dataset$Source[dataset$Source == '6_GBWT'] <- 'Giraffe_6'
dataset$Source[dataset$Source == 'ALL_GBWT'] <- 'Giraffe_All'
dataset$Source[dataset$Source == 'M'] <- 'Minimap2_1'
dataset$Source[dataset$Source == 'M3'] <- 'Minimap2_3'
dataset$Source[dataset$Source == 'M6'] <- 'Minimap2_6'
dataset$Source[dataset$Source == 'MALL'] <- 'Minimap2_All'

  j <- ggplot(data = dataset, aes(x=Caller, y=Rate, color=Caller)) +
    geom_boxplot(alpha=0.5,outlier.shape=NA) + geom_jitter(size=1,alpha=0.4)+ 
    theme(axis.text.x=element_text(size=10, angle=90)) +ylab(label1)+
    ylim(0,8)+ facet_wrap(~Source, ncol = 4) 
  pdf(nameSummary, width=5.5, height=4)
  print( j )
  dev.off()

  p <- ggplot(data = dataset, aes(x=Caller, y=Rate, color=Caller)) +
    geom_point( ) + theme(axis.text.x=element_text(size=10, angle=90))+
    ylab(label1) +ylim(0,8)
  pdf(nameSamples, width=5 + dim(dataset)[1]/99, height=5 + dim(dataset)[1]/111)
  print( p + facet_wrap(~Sample, ncol = 8) )
  dev.off() }

tstv_amplicon <- data.frame(read.csv("TSTV.AMPLICON.txt", sep="\t"))
colnames(tstv_amplicon) <- c("Sample", "Caller", "Source", "Rate")
tstv_amplicon <- subset(tstv_amplicon, Rate>0)
tstv_amplicon <- subset(tstv_amplicon, Caller !="FB2" & Caller !="BCF2")
print("tstv_amplicon")
median(na.omit(subset(tstv_amplicon, Caller=="BCF"  )$Rate)) # 2
median(na.omit(subset(tstv_amplicon, Caller=="FB"  )$Rate)) # 0.44
median(na.omit(subset(tstv_amplicon, Caller=="VG")$Rate)) # 0.43
make_plots(tstv_amplicon, "tstv_amplicon", "Ti/Tv rate Amplicon")

tstv_metagenomic <- data.frame(read.csv("TSTV.METAGENOMIC.txt", sep="\t"))
colnames(tstv_metagenomic) <- c("Sample", "Caller", "Source", "Rate")
tstv_metagenomic <- subset(tstv_metagenomic, Rate>0)
tstv_metagenomic <- subset(tstv_metagenomic, Caller !="FB2" & Caller !="BCF2")
print("tstv_metagenomic")
median(na.omit(subset(tstv_metagenomic, Caller=="BCF"  )$Rate)) # 2.07
median(na.omit(subset(tstv_metagenomic, Caller=="FB"  )$Rate)) # 0.16
median(na.omit(subset(tstv_metagenomic, Caller=="VG")$Rate)) # 0.93
make_plots(tstv_metagenomic, "tstv_metagenomic", "Ti/Tv rate Metgenomic")

tstv_wgs <- data.frame(read.csv("TSTV.WGS.txt", sep="\t"))
colnames(tstv_wgs) <- c("Sample", "Caller", "Source", "Rate")
tstv_wgs <- subset(tstv_wgs, Rate>0)
tstv_wgs <- subset(tstv_wgs, Caller !="FB2" & Caller !="BCF2")
print("tstv_wgs")
median(na.omit(subset(tstv_wgs, Caller=="BCF"  )$Rate)) # 2.83
median(na.omit(subset(tstv_wgs, Caller=="FB"  )$Rate)) # 0.82
median(na.omit(subset(tstv_wgs, Caller=="VG")$Rate)) # 0.99
make_plots(tstv_wgs, "tstv_wgs", "Ti/Tv rate WGS")

tstv_ont <- data.frame(read.csv("TSTV.ONT.txt", sep="\t"))
colnames(tstv_ont) <- c("Sample", "Caller", "Source", "Rate")
tstv_ont <- subset(tstv_ont, Rate>0)
tstv_ont <- subset(tstv_ont, Caller !="FB2" & Caller !="BCF2")
print("tstv_ont")
median(na.omit(subset(tstv_ont, Caller=="BCF" | Caller=="BCF2")$Rate)) # 0.59
median(na.omit(subset(tstv_ont, Caller=="FB" | Caller=="FB2")$Rate)) # 3.55
median(na.omit(subset(tstv_ont, Caller=="VG")$Rate)) # 2.78
make_plots(tstv_ont, "tstv_ont", "Ti/Tv rate ONT")

options(dplyr.print_max = 1e9)
print("amplicon")
tstv_amplicon %>% group_by(Caller, Source) %>% summarise(Median_Rate = median(Rate, na.rm =T)) %>% arrange(Source, Caller)
print("wgs")
tstv_wgs %>% group_by(Caller, Source) %>% summarise(Median_Rate = median(Rate, na.rm =T)) %>% arrange(Source, Caller)
print("ont")
tstv_ont %>% group_by(Caller, Source) %>% summarise(Median_Rate = median(Rate, na.rm =T)) %>% arrange(Source, Caller)

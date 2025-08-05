#!/usr/bin/perl
use strict;
use warnings;

# VG-MAP Array of datasets and folders
my @vgmap_datasets = qw( giraffe3 ); 
my @folders = qw(   SCREENED_MERGED_1_GBWT SCREENED_MERGED_3_GBWT SCREENED_MERGED_6_GBWT );

for my $i (0..$#vgmap_datasets) {
    my $vgmap = $vgmap_datasets[$i];
    my $folder = $folders[$i];
    open(OUT, ">rates.new4.$vgmap.txt") or die "Cannot open output file: $!";
    print OUT "Sample\tLibrary_type\tPVG\tMinimap2\tTotal_all\tPercent\tPercentM\n";
    open my $metadata, '<', 'metadata2.csv' or die "Could not open 'metadata2.csv': $!";
    
    # Skip the header line if necessary
    my $x = 0;
    while (my $line = <$metadata>) {
        chomp $line;
        print $line, " $folder\n";
            my @r = split(/\s+/, $line, 2);
            my $sample = $r[0];  # Define the sample name
            $r[1] =~ s/_SE//g;
            $r[1] =~ s/_PE//g;
            
	    # Minimap2         CONCAT/$sample.$vgmap.m.vcf.gz
      # PVG              FINAL_SNPS_2/$sample.$vgmap.vcf.gz FINAL_SNPS_2/$sample.giraffe1.vcf.gz
	    system("bcftools merge --force-samples -m all -O v FINAL_SNPS_2/$sample.$vgmap.vcf.gz FINAL_SNPS_2/$sample.giraffe1.vcf.gz | bcftools view -v snps > FINAL_SNPS_2/$sample.giraffe13.vcf");
            system("bgzip -f FINAL_SNPS_2/$sample.giraffe13.vcf");
            system("tabix -p vcf -f FINAL_SNPS_2/$sample.giraffe13.vcf.gz");
            system("ls -lt FINAL_SNPS_2/$sample.giraffe13.vcf.gz"); # SNPS in 1 and 3
	
            # PVG + Minimap2   
	    system("bcftools merge --force-samples -m all -O v FINAL_SNPS_2/$sample.giraffe13.vcf.gz CONCAT/$sample.$vgmap.m.vcf.gz  | bcftools view -v snps > FINAL_SNPS_BOTH/$sample.giraffe13m.vcf");
            system("bgzip -f FINAL_SNPS_BOTH/$sample.giraffe13m.vcf");
            system("tabix -p vcf -f FINAL_SNPS_BOTH/$sample.giraffe13m.vcf.gz");
            system("ls -lt FINAL_SNPS_BOTH/$sample.giraffe13m.vcf.gz");  # SNPs in 13 and m
            
            # get 13 SNPs not in m
            system("bcftools isec FINAL_SNPS_2/$sample.giraffe13.vcf.gz CONCAT/$sample.$vgmap.m.vcf.gz -O z -n~10 -o FINAL_SNPS_2/$sample.giraffe13.unique.vcf.gz ");
            # get 1 SNPs not in m
            system("bcftools isec FINAL_SNPS_2/$sample.giraffe1.vcf.gz CONCAT/$sample.$vgmap.m.vcf.gz -O z -n~10 -o FINAL_SNPS_2/$sample.giraffe1.unique.vcf.gz ");  
             # get m SNPs not in 1
            system("bcftools isec CONCAT/$sample.$vgmap.m.vcf.gz FINAL_SNPS_2/$sample.giraffe1.vcf.gz -O z -n~10 -o CONCAT/$sample.$vgmap.m.unique1.vcf.gz  ");                  
            # get m SNPs not in 13
            system("bcftools isec CONCAT/$sample.$vgmap.m.vcf.gz FINAL_SNPS_2/$sample.giraffe13.vcf.gz -O z -n~10 -o CONCAT/$sample.$vgmap.m.unique.vcf.gz ");
                     
            my $mini = `gzip -dc CONCAT/$sample.m.vcf.gz | grep -cv "#"`;
            chomp($mini);  # Minimap2
            
            my $mergedall2 = `gzip -dc FINAL_SNPS_2/$sample.giraffe13.vcf.gz | grep -cv "#"`;
            chomp($mergedall2); # PVG 
                       
              my $mergedall = `gzip -dc FINAL_SNPS_BOTH/$sample.giraffe13m.vcf.gz | grep -cv "#"`;
             chomp($mergedall); # PVG and M
            
            my $pc="NA"; # PVG % 
            my $pcm="NA"; # Minimap2 %
            if($mergedall>0){ $pc = sprintf("%.2f",(100*($mergedall-$mini))/($mergedall) );   # all - M = PVG
                       $pcm= sprintf("%.2f",(100*($mergedall-$mergedall2))/($mergedall) ) ; } # all - PVG = M
            print OUT "$sample\t$r[1]\t$mergedall2\t$mini\t$mergedall\t$pc\t$pcm\n";   }
    close $metadata;
    close OUT;
}

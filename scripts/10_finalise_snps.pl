#!/usr/bin/perl
use strict;
use warnings;

# VG-MAP Array of datasets and folders
my @vgmap_datasets = qw( giraffe3 giraffe6 );
 my @vgmap_datasets = qw(vgmap1 vgmap3 vgmap6 vgmapall giraffe1 giraffe3 giraffe6 giraffeall);
# my @vgmap_datasets = qw( giraffe1 );
my @folders = qw(SCREENED_MERGED_1 SCREENED_MERGED_3 SCREENED_MERGED_6 SCREENED_MERGED_ALL
                 SCREENED_MERGED_1_GBWT SCREENED_MERGED_3_GBWT SCREENED_MERGED_6_GBWT SCREENED_MERGED_ALL_GBWT);

for my $i (0..$#vgmap_datasets) {
    my $vgmap = $vgmap_datasets[$i];
    my $folder = $folders[$i];
    
    # Open output file for each vgmap dataset
    open(OUT, ">rates.new3.$vgmap.txt") or die "Cannot open output file: $!";
    print OUT "Sample\tLibrary_type\tPVG\tMinimap2\tTotal_all\tPercent\tPercentM\n";
    
    # Open metadata file
    open my $metadata, '<', 'metadata2.csv' or die "Could not open 'metadata2.csv': $!";
    
    # Skip the header line if necessary
    my $x = 0;
    while (my $line = <$metadata>) {
        chomp $line;
        print $line, " $folder\n";
       #  if ($x > 0) {  # Skip header
            my @r = split(/\s+/, $line, 2);
            my $sample = $r[0];  # Define the sample name
            $r[1] =~ s/_SE//g;
            $r[1] =~ s/_PE//g;
            
            # Get Minimap2 SNPs - $sample.m.vcf.gz
            system("bcftools view -v snps SCREENED_MERGED_M/BCF_$sample.vcf.gz > CONCAT/$sample.m.bcf.vcf");
            system("bgzip -f CONCAT/$sample.m.bcf.vcf");
            system("tabix -p vcf -f CONCAT/$sample.m.bcf.vcf.gz");
            system("ls -lt CONCAT/$sample.m.bcf.vcf.gz");

            # Get Minimap2 SNPs - $sample.m.vcf.gz -FB
            system("bcftools view -v snps SCREENED_MERGED_M/FB_$sample.vcf.gz > CONCAT/$sample.m.fb.vcf");
            system("bgzip -f CONCAT/$sample.m.fb.vcf");
            system("tabix -p vcf -f CONCAT/$sample.m.fb.vcf.gz");
            system("ls -lt CONCAT/$sample.m.fb.vcf.gz");

	    # Merge Minimap2 BCF and FB
            system("bcftools merge --force-samples -m all -O v CONCAT/$sample.m.bcf.vcf.gz CONCAT/$sample.m.fb.vcf.gz | bcftools view -v snps > CONCAT/$sample.$vgmap.m.vcf");
            system("bgzip -f CONCAT/$sample.$vgmap.m.vcf");
            system("tabix -p vcf -f CONCAT/$sample.$vgmap.m.vcf.gz");
            system("ls -lt CONCAT/$sample.$vgmap.m.vcf.gz");
                 
            # Create concat PVG BCF SNPs file from BCFtools
            system("bcftools concat -a -D $folder/BCF_$sample.vcf.gz -O v | bcftools norm -d all | bcftools view -v snps > CONCAT/$sample.$vgmap.temp");
            system("bgzip -f CONCAT/$sample.$vgmap.temp");
            system("tabix -p vcf -f CONCAT/$sample.$vgmap.temp.gz");
            system("ls -lt CONCAT/$sample.$vgmap.temp.gz");
                  
            # Merge PVG-BCF above with Minimap2
            system("bcftools merge --force-samples -m all -O v CONCAT/$sample.m.vcf.gz CONCAT/$sample.$vgmap.temp.gz | bcftools view -v snps > FINAL_SNPS/$sample.$vgmap.vcf");
            system("bgzip -f FINAL_SNPS/$sample.$vgmap.vcf");
            system("tabix -p vcf -f FINAL_SNPS/$sample.$vgmap.vcf.gz");
            system("ls -lt FINAL_SNPS/$sample.$vgmap.vcf.gz");
                 
            # Create concat PVG SNPs file from FB
            system("bcftools concat -a -D $folder/FB_$sample.vcf.gz -O v | bcftools norm -d all | bcftools view -v snps > CONCAT_FB/$sample.$vgmap.temp");
            system("bgzip -f CONCAT_FB/$sample.$vgmap.temp");
            system("tabix -p vcf -f CONCAT_FB/$sample.$vgmap.temp.gz");
            system("ls -lt CONCAT_FB/$sample.$vgmap.temp.gz");
                
            # Merge PVG-FB above with Minimap2
            system("bcftools merge --force-samples -m all -O v CONCAT/$sample.m.vcf.gz CONCAT_FB/$sample.$vgmap.temp.gz | bcftools view -v snps > FINAL_SNPS_FB/$sample.$vgmap.vcf");
            system("bgzip -f FINAL_SNPS_FB/$sample.$vgmap.vcf");
            system("tabix -p vcf -f FINAL_SNPS_FB/$sample.$vgmap.vcf.gz");
            system("ls -lt FINAL_SNPS_FB/$sample.$vgmap.vcf.gz");
            
            # Create concat PVG SNPs file from VG 
#        system("bcftools concat -a -D $folder/VG_$sample.vcf.gz -O v | bcftools norm -d all | bcftools view -v snps > CONCAT_VG/$sample.$vgmap.temp"); 
#        system("bgzip -f CONCAT_VG/$sample.$vgmap.temp");
#        system("tabix -p vcf -f CONCAT_VG/$sample.$vgmap.temp.gz");
#        system("ls -lt CONCAT_VG/$sample.$vgmap.temp.gz");
                       
            # VG Edit the VCF files to replace '##FORMAT=<ID=AD,Number=R,' with '##FORMAT=<ID=AD,Number=.,
#            system("zcat CONCAT_VG/$sample.$vgmap.temp.gz | sed 's/##FORMAT=<ID=AD,Number=R,/##FORMAT=<ID=AD,Number=.,/' | sed 's/KX894507,/KX894508,/' > CONCAT_VG/$sample.$vgmap.temp.vcf ");
#            system("bgzip -f CONCAT_VG/$sample.$vgmap.temp.vcf");
#            system("mv CONCAT_VG/$sample.$vgmap.temp.vcf.gz CONCAT_VG/$sample.$vgmap.temp.gz");
#            system("tabix -p vcf -f CONCAT_VG/$sample.$vgmap.temp.gz");
        
            # Merge PVG-VG above with Minimap2
#            system("bcftools merge --force-samples -m all -O v CONCAT/$sample.m.vcf.gz CONCAT_VG/$sample.$vgmap.temp.gz | bcftools view -v snps > FINAL_SNPS_VG/$sample.$vgmap.vcf"); 
#            system("bgzip -f FINAL_SNPS_VG/$sample.$vgmap.vcf");
#            system("tabix -p vcf -f FINAL_SNPS_VG/$sample.$vgmap.vcf.gz");
#            system("ls -lt FINAL_SNPS_VG/$sample.$vgmap.vcf.gz");
            
            # VG Edit the VCF files to replace '##FORMAT=<ID=AD,Number=R,' with '##FORMAT=<ID=AD,Number=.,
#            system("zcat CONCAT_FB/$sample.$vgmap.temp.gz | sed 's/##FORMAT=<ID=AD,Number=R,/##FORMAT=<ID=AD,Number=.,/' > CONCAT_FB/$sample.$vgmap.temp.vcf ");
#            system("bgzip -f CONCAT_FB/$sample.$vgmap.temp.vcf");
#            system("mv CONCAT_FB/$sample.$vgmap.temp.vcf.gz CONCAT_FB/$sample.$vgmap.temp.gz");
#            system("tabix -p vcf -f CONCAT_FB/$sample.$vgmap.temp.gz");
        
            # Merge PVG_BCF above with PVG_FB  
            system("bcftools merge --force-samples -m all -O v CONCAT/$sample.$vgmap.temp.gz CONCAT_FB/$sample.$vgmap.temp.gz  | bcftools view -v snps | sed 's/KX894507,/KX894508,/' > FINAL_SNPS_2/$sample.$vgmap.vcf"); 
            system("bgzip -f FINAL_SNPS_2/$sample.$vgmap.vcf ");
            system("tabix -p vcf -f FINAL_SNPS_2/$sample.$vgmap.vcf.gz");
            system("ls -lt FINAL_SNPS_2/$sample.$vgmap.vcf.gz");
                               
            # Merge PVG_BCF above with PVG_FB with Minimap2
            system("bcftools merge --force-samples -m all -O v CONCAT/$sample.m.vcf.gz CONCAT/$sample.$vgmap.temp.gz  CONCAT_FB/$sample.$vgmap.temp.gz | bcftools view -v snps > FINAL_SNPS_BOTH/$sample.$vgmap.vcf"); 
            system("bgzip -f FINAL_SNPS_BOTH/$sample.$vgmap.vcf");
            system("tabix -p vcf -f FINAL_SNPS_BOTH/$sample.$vgmap.vcf.gz");
            system("ls -lt FINAL_SNPS_BOTH/$sample.$vgmap.vcf.gz");
                      
            # Get the counts for BCFtools only
            my $pvg = `gzip -dc CONCAT/$sample.$vgmap.temp.gz | grep -cv "#"`;
            chomp($pvg);  # Count SNPs PVG-based method
            
            my $pvgfb = `gzip -dc CONCAT_FB/$sample.$vgmap.temp.gz | grep -cv "#"`;
            chomp($pvgfb);  # Count SNPs PVG-based method
            
          #    my $pvgvg = `gzip -dc CONCAT_VG/$sample.$vgmap.temp.gz | grep -cv "#"`;
          #    chomp($pvgvg);  # Count SNPs PVG-based method
            
            my $mini = `gzip -dc CONCAT/$sample.m.vcf.gz | grep -cv "#"`;
            chomp($mini);  # Minimap2
            
            my $merged =   `gzip -dc FINAL_SNPS/$sample.$vgmap.vcf.gz | grep -cv "#"`;
            chomp($merged);  # PVG BCF only
            
            my $mergedfb = `gzip -dc FINAL_SNPS_FB/$sample.$vgmap.vcf.gz | grep -cv "#"`;
            chomp($mergedfb); # PVG FB only

           #   my $mergedvg = `gzip -dc FINAL_SNPS_VG/$sample.$vgmap.vcf.gz | grep -cv "#"`;
           #   chomp($mergedvg); # PVG VG only        
            
            my $mergedall2 = `gzip -dc FINAL_SNPS_2/$sample.$vgmap.vcf.gz | grep -cv "#"`;
            chomp($mergedall2); # PVG all
                       
              my $mergedall = `gzip -dc FINAL_SNPS_BOTH/$sample.$vgmap.vcf.gz | grep -cv "#"`;
             chomp($mergedall); # PVG and M
            
            my $pc="NA"; # PVG % 
            my $pcm="NA"; # Minimap2 %
            if($mergedall>0){ $pc = sprintf("%.2f",(100*($mergedall-$mini))/($mergedall) );   # all - M = PVG
                       $pcm= sprintf("%.2f",(100*($mergedall-$mergedall2))/($mergedall) ) ; } # all - PVG = M
                   
            # Write results to the output file
            print OUT "$sample\t$r[1]\t$mergedall2\t$mini\t$mergedall\t$pc\t$pcm\n";
       #  }
        #  $x++;
    }
    close $metadata;
    close OUT;
}

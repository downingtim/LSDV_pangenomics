# TS/TV rates

open(ALL, "ls SCREENED_*/*vcf*z | grep -v \"255939\" | grep -v \"23347779 \" | ");

open(META,"acc_list.metagenomic.txt");
@a3=<META>;
%h3={};
for $e (0..$#a3){ chomp($a3[$e]);    $h3{$a3[$e]}=1; }
close(META);
open(OUT3, ">./TSTV.METAGENOMIC.txt");

open(W,"acc_list.wgs.txt");
@a4=<W>;
%h4={};
for $e (0..$#a4){ chomp($a4[$e]);    $h4{$a4[$e]}=1; }
close(W);
open(OUT4, ">./TSTV.WGS.txt");

open(AMP,"acc_list.amplicon.txt");
@a22=<AMP>;
%h22={};
for $e (0..$#a22){ chomp($a22[$e]);    $h22{$a22[$e]}=1; }
#print "$h22{$a22[$e]} -> $a22[$e] -> $e \n"; }
close(AMP);
open(OUT22, ">./TSTV.AMPLICON.txt"); 

while(<ALL>){
    chomp;
    @r=split/\//,$_;
    $r[0]=~ s/SCREENED_MERGED_//g;
    @z=split/\./,$r[1];
    @z2=split/_/,$z[0]; 
    open(HIT, "bcftools stats $_ | grep TSTV | grep -v \"#\" | ");
    $in =<HIT>;
    @in2=split/\s+/,$in;

    #print "bcftools stats $_ | grep TSTV | grep -v \"#\" | -> $in2[4] \n";
    
    if($h3{$z2[1]}) {  # if metagenomic 
      	print OUT3 "$z2[1]\t$z2[0]\t$r[0]\t$in2[4]\n"; } 
    elsif($h4{$z2[1]}) {  # if wgs 
      	print OUT4 "$z2[1]\t$z2[0]\t$r[0]\t$in2[4]\n";  } 
    if($h22{$z2[1]}) {  # if in amplicon   
      	print OUT22 "$z2[1]\t$z2[0]\t$r[0]\t$in2[4]\n";   } 
    close(HIT);}
close(ALL);
close(OUT1);
close(OUT3);
close(OUT4);

open(ONT, "ls SCREENED_*/*vcf*z | grep -E \"255939|23347779 \" | ");  
open(ONT1, ">./TSTV.ONT.txt");

while(<ONT>){
    chomp;
    open(HIT, "bcftools stats $_ | grep TSTV | grep -v \"#\"  |");
    $in =<HIT>;
    @in2=split/\s+/,$in;
    @r=split/\//,$_;
    $r[0]=~ s/SCREENED_MERGED_//g;
    @z=split/\./,$r[1];
    @z2=split/_/,$z[0];
    print ONT1 "$z2[1]\t$z2[0]\t$r[0]\t$in2[4]\n";
    close(HIT);}
close(ONT1); 

# 

open(ALL, "ls /mnt/lustre/RDS-ephemeral/downing/LSDV/LSDV*/*_VCF_FILES/*norm*vcf*z | grep -v \"SRR10394925t\" | grep -v \"255939\" | grep -v \"23347779 \" | grep -v FINAL | "); # all files

open(META,"acc_list.metagenomic.txt");
@a3=<META>;
%h3={};
for $e (0..$#a3){ chomp($a3[$e]);    $h3{$a3[$e]}=1; }
close(META);
open(SNPS3, ">./SNPS.METAGENOMIC.txt");
open(B3, ">./BAMSTATS.METAGENOMIC.txt");
open(G3, ">./GAMSTATS.METAGENOMIC.txt");

open(W,"acc_list.wgs.txt");
@a4=<W>;
%h4={};
for $e (0..$#a4){ chomp($a4[$e]);    $h4{$a4[$e]}=1; }
close(W);
open(SNPS4, ">./SNPS.WGS.txt");
open(B4, ">./BAMSTATS.WGS.txt");
open(G4, ">./GAMSTATS.WGS.txt");

open(AMP,"acc_list.amplicon.txt");
@a22=<AMP>;
%h22={};
for $e (0..$#a22){ chomp($a22[$e]);    $h22{$a22[$e]}=1; }
close(AMP);
open(SNPS22, ">./SNPS.AMPLICON.txt");
open(B22, ">./BAMSTATS.AMPLICON.txt");
open(G22, ">./GAMSTATS.AMPLICON.txt");

print B1 "Sample\tSource\tRate\n";
print B2 "Sample\tSource\tRate\n";
print B22 "Sample\tSource\tRate\n";
print B3 "Sample\tSource\tRate\n";
print B4 "Sample\tSource\tRate\n";
print G1 "Sample\tSource\tRate\n";
print G2 "Sample\tSource\tRate\n";
print G22 "Sample\tSource\tRate\n";
print G3 "Sample\tSource\tRate\n";
print G4 "Sample\tSource\tRate\n";

while(<ALL>){
    chomp;
    @r=split/\//,$_;
    $r[7]=~ s/_VCF_FILES//g;
    @z=split/\./,$r[8];
    
    $comm=`gzip -dc $_ | awk ' { print \$2 } ' | sort | uniq | grep -c -v "#" > temp.illumina `;
    open(HIT2, "temp.illumina");
    $in2 =<HIT2>;
    @r2=split/\//,$_;
    $r2[7]=~ s/_VCF_FILES//g;
    @z=split/\./,$r2[8];

    if($h{$z[0]}) {  # if in amplicon PE      # print "$z[0] -> $h{$z[0]}\n";
	print OS1 "$z[0]\t$r[7]\t$r[6]\t$in2[4]\t$in2";
        print SNPS1 "$z[0]\t$r2[7]\t$r2[6]\t$in2";  } 
    elsif($h2{$z[0]}) {  # if in amplicon SE
	print OS2 "$z[0]\t$r[7]\t$r[6]\t$in2[4]\t$in2";
        print SNPS2 "$z[0]\t$r2[7]\t$r2[6]\t$in2";  } 
    elsif($h3{$z[0]}) {  # if metagenomic
	print OS3 "$z[0]\t$r[7]\t$r[6]\t$in2[4]\t$in2";
        print SNPS3 "$z[0]\t$r2[7]\t$r2[6]\t$in2";  } 
    elsif($h4{$z[0]}) {  # if wgs
	print OS4 "$z[0]\t$r[7]\t$r[6]\t$in2[4]\t$in2";
        print SNPS4 "$z[0]\t$r2[7]\t$r2[6]\t$in2";  } 

    if($h22{$z[0]}) {  # if in amplicon SE
	print OS22 "$z[0]\t$r[7]\t$r[6]\t$in2[4]\t$in2";
        print SNPS22 "$z[0]\t$r2[7]\t$r2[6]\t$in2";  } 
    close(HIT);
    close(HIT2);
}
close(ALL);
close(OUT);
close(SNPS1);
close(SNPS2);
close(SNPS22);
close(SNPS3);
close(SNPS4);
close(OS1);
close(OS2);
close(OS22);
close(OS3);
close(OS4);

# also: cat SNPS.AMPLICON.SE.txt SNPS.AMPLICON.PE.txt SNPS.WGS.txt SNPS.METAGENOMIC.txt > SNPS.ILLUMINA.txt

open(BAM, "ls /mnt/lustre/RDS-ephemeral/downing/LSDV/LSDVV*/BAM_FILE*/*flagstat /mnt/lustre/RDS-ephemeral/downing/LSDV/LSDVG*/BAM_FILE*/*flagstat /mnt/lustre/RDS-ephemeral/downing/LSDV/LSDV1/BAM_FILES2/*flagstat | grep -v \"SRR10394925t\" | grep -v \"255939\" | grep -v \"23347779 \"  | grep -v FINAL | "); # all files

while(<BAM>){
    chomp;
    open(TT, "grep mapped $_ | grep \% | ");
    @x=split/\//,$_;
    $x[8]=~ s/.flagstat//g;
    @a=<TT>;
    chomp($a[0]);
    @r=split/\s+/,$a[0];
    $r[4]=~ s/\(//g;
    $r[4]=~ s/\%//g;
  #  if($x[7]=~ /S2/){ $x[6]=""; $x[8]=""; $r[4]="";}     else { 
    if($h{$x[8]}) {  # if in amplicon PE     
    # print "$z[0] -> $h{$z[0]}\n";
	print B1 "$x[8]\t$x[6]\t$r[4]\n"; }
    elsif($h2{$x[8]}) {  # if in amplicon SE
	print B2 "$x[8]\t$x[6]\t$r[4]\n"; }
    elsif($h3{$x[8]}) {  # if metagenomic
	print B3 "$x[8]\t$x[6]\t$r[4]\n"; }
    elsif($h4{$x[8]}) {  # if wgs
	print B4 "$x[8]\t$x[6]\t$r[4]\n"; }# }

    if($h22{$x[8]}) {  # if in amplicon 
	     print B22 "$x[8]\t$x[6]\t$r[4]\n"; }
    close(TT);
}

open(GAM, "ls /mnt/lustre/RDS-ephemeral/downing/LSDV/LSDVG*/GAM_FILES/*stat  /mnt/lustre/RDS-ephemeral/downing/LSDV/LSDVVG*/GAM_FILES/*stat| grep -v \"SRR10394925t\" | grep -v \"255939\" | grep -v \"23347779 \"  | "); # all files

while(<GAM>){
    chomp;
    @x=split/\//,$_;
    $x[8]=~ s/.stat//g;

    open(TT, "grep aligne $_ | ");
    @a=<TT>;
    chomp($a[0]);
    @r=split/\s+/,$a[0];
    close(TT);
    open(TT2, "grep primary $_ | ");
    @a2=<TT2>;
    chomp($a2[0]);
    @r2=split/\s+/,$a2[0];
    close(TT2);

    if($h{$x[8]}) {  # if in amplicon PE      # print "$z[0] -> $h{$z[0]}\n";
        print G1 "\n$x[8]\t$x[6]\t";
          if(($r[2]>0)&&($r2[2])){ print G1 sprintf("%.2f", (100*$r[2])/$r2[2]); }
	else { print G1 "NA"; } }
    elsif($h2{$x[8]}) {  # if in amplicon SE
        print G2 "\n$x[8]\t$x[6]\t";
          if(($r[2]>0)&&($r2[2])){ print G2 sprintf("%.2f", (100*$r[2])/$r2[2]); }
	else { print G2 "NA"; } }
    elsif($h3{$x[8]}) {  # if metagenomic
        print G3 "\n$x[8]\t$x[6]\t";
          if(($r[2]>0)&&($r2[2])){ print G3 sprintf("%.2f", (100*$r[2])/$r2[2]); }
	else { print G3 "NA"; } }
    elsif($h4{$x[8]}) {  # if wgs
        print G4 "\n$x[8]\t$x[6]\t";
          if(($r[2]>0)&&($r2[2])){ print G4 sprintf("%.2f", (100*$r[2])/$r2[2]); }
	else { print G4 "NA"; } }

    if($h22{$x[8]}) {  # if in amplicon 
        print G22 "\n$x[8]\t$x[6]\t";
          if(($r[2]>0)&&($r2[2])){ print G22 sprintf("%.2f", (100*$r[2])/$r2[2]); }
	else { print G22 "NA"; } }

    close(TT); }

# FB ONT only #ls  FB_VCF_FILES/*.norm*vcf.gz | grep -E "255939|23347779" | 

open(ONT, "ls /mnt/lustre/RDS-ephemeral/downing/LSDV/LSDVV*/*VCF_FILES/*norm*vcf*z /mnt/lustre/RDS-ephemeral/downing/LSDV/LSDVG*/*VCF_FILES/*norm*vcf*z | grep -v \"SRR10394925t\" | grep -E \"255939|23347779 \" | "); # all files
open(ONT2, ">./SNPS.ONT.txt");

while(<ONT>){
    chomp;
    $comm=`gzip -dc $_ | awk ' { print \$2 } ' | sort | uniq | grep -c -v "#" > temp.ont `;
    open(HIT, "temp.ont");
    $in =<HIT>;
    @r=split/\//,$_;
    $r[7]=~ s/_VCF_FILES//g;
    @z=split/\./,$r[8];
    print ONT2 "$z[0]\t$r[7]\t$r[6]\t$in"; 
}
close(ONT1);
close(ONT2);

open(INTO, "acc_list.ont.txt");
open(O, ">./BAMSTATS.ONT.txt");
open(O2, ">./GAMSTATS.ONT.txt");
print O  "Sample\tSource\tRate\n";
print O2 "Sample\tSource\tRate\n";

@xx=<INTO>;
@type=qw(LSDVVG LSDV1 LSDV1_3 LSDV1_6 LSDV1_ALL LSDVG LSDVVG_3 LSDVG_3 LSDVVG_6 LSDVG_6 LSDVVG_ALL LSDVG_ALL );

for $e (0..$#xx){
      chomp($xx[$e]);
      for $type1 (0..$#type){
          $name = "/mnt/lustre/RDS-e*/downing/LSDV/".$type[$type1]."/BAM_FILES/".$xx[$e].".flagstat";
          open(TT, "grep mapped $name | grep \% | ");
          @a=<TT>;
          chomp($a[0]);
          @r=split/\s+/,$a[0];
          $r[4]=~ s/\(//g;
          $r[4]=~ s/\%//g;
          print O "$xx[$e]\t$type[$type1]\t$r[4]\n";
          close(TT);      }
}

for $e (0..$#xx){
      chomp($xx[$e]);
      for $type1 (0..$#type){
          $name = "/mnt/lustre/RDS-e*/downing/LSDV/".$type[$type1]."/GAM_FILES/".$xx[$e].".stat";
          open(TT, "grep aligne $name | ");
          @a=<TT>;
          chomp($a[0]);
          @r=split/\s+/,$a[0];
          close(TT);
          open(TT2, "grep primary $name | ");
          @a2=<TT2>;
          chomp($a2[0]);
          @r2=split/\s+/,$a2[0];
          close(TT2);
          print O2 "$xx[$e]\t$type[$type1]\t";
          if(($r[2]>0)&&($r2[2])){ print O2 sprintf("%.2f", (100*$r[2])/$r2[2]); }
          else { print OUT2 "NA"; }
          print O2 "\n";      }
}

close(INTO);
close(O);
close(O2);

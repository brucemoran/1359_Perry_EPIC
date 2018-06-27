#! /bin/bash
WORKDIR="$1"

if [[ ! -d "$WORKDIR" ]]; then
  mkdir -p "$WORKDIR";
fi
cd "$WORKDIR"

#get files
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip
wget http://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/infinium-methylationepic-v1-0-missing-legacy-cpg-b3-vs-b2-annotations.zip

##(g)unzip those files
unzip -o infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip
unzip -o infinium-methylationepic-v1-0-missing-legacy-cpg-b3-vs-b2-annotations.zip

#########################################
##create annotation files in bed format##
#########################################
##Infinium EPIC bed:

##making MethylationEPIC_v-1-0_B4.csv into MethylationEPIC_v-1-0_B4.bed
##want 4 ';' delim fields, along with 1,2,3 of chr, start, end
##anno delim ';' -> probe;gene;transcript;feature_type
##NB if your chips start 2011...something... then you don't want to use the Legacy B2 Missing probes
##ours start 2007...something... so I do use it
cp MethylationEPIC_v-1-0_B4.csv MethylationEPIC_v-1-0_B4.Legacy_B2.csv
tail -n+8 "MethylationEPIC Missing Legacy CpG (v1.0_B3 vs. v1.0_B2) Annotations.csv" >> MethylationEPIC_v-1-0_B4.Legacy_B2.csv
tail -n+8 MethylationEPIC_v-1-0_B4.Legacy_B2.csv | perl -ane '@s=split(/\,/);
  $cn=$s[11]; $cn=~s/chr//;
  $st=$s[12]; $en=$st; $en++;
  $str=$s[14];
  if($str eq "R"){$stro="-"} if($str eq "F"){$stro="+"}
  if(($str ne "F") && ($str ne "R")){$stro="*"}
  $gn=$s[15];
  $tn=$s[16];
  $bd=$s[17];
  $in=$s[19];
  if($in eq ""){$in="Open_Sea";}
  $gn=~s/\;/\,/g; $tn=~s/\;/\,/g; $bd=~s/\;/\,/g; $in=~s/\;/\,/g;
  if($cn ne "CHR"){print "$cn\t$s[12]\t$en\t$stro\t$s[0];$gn;$tn;$bd;$in\n";}' | gsort -V | perl -ane 'if(scalar(@F)==5){if($F[0] ne "CHR"){print $_;}}' > MethylationEPIC_v-1-0_B4.Legacy_B2.bed

##zip CSVs and delete
zip EPIC_annotation.zip *csv *csv.* *zip
rm *csv *csv.* infin*

#!/bin/sh

data=$1

echo ${data}

#plink --bfile split_by_chr/${data}.chr6 --snp rs17671350 --window 50 --out region1 --make-bed
#plink --bfile split_by_chr/${data}.chr6 --snp rs11752328 --window 50 --out region2 --make-bed

#cat region1.bim region2.bim > region.bim
#awk {'print $2'} region.bim > snps.txt
#start=`head -1 snps.txt | gawk '{ print \$0 }'`
#end=`tail -1 snps.txt | gawk '{print \$0}'`

#extract region covered by the Mcaroll Lab reference panel +/- 50kb. Coordinates from build 7.
plink --bfile split_by_chr/${data}.chr6 --from-bp 24944177 --to_bp 33940574 --out region.${data}.chr6 --make-bed 

#write snps that have only one allele to excluded_snps.txt
#gawk '$5==0 { print $2 }' all.region.camh.chr6.bim > excluded_snps.txt

#exclude selected snps
#plink --bfile all.region.camh.chr6 --exclude excluded_snps.txt --out region.camh.chr6 --make-bed

#clean up files
#rm region1*
#rm region2*
#rm snps.txt
#rm region.bim

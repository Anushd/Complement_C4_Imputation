library(xlsx)
library(data.table)
library(dplyr)
library(plyr)

setwd("/Users/anush/Documents/Projects/MGH/imputation/ref_panel")

ref <- data.frame(fread("MHC_haplotypes_CEU_HapMap3_ref_panel.bgl", header=T ))
ordered_ref <- ref[order(ref$id),]

kgp1 <- data.frame(fread("snp147_chr6_bp_only.txt", header=F))
colnames(kgp1) <- c("id","bp")

kgp2 <- data.frame(fread("1kgp3_chr6.bim", header=F))
colnames(kgp2) <- c("chr","id","-","bp","a1","a2")
kgp2 <- kgp2[c("id","bp")]

kgp3 <- data.frame(fread("1kgp1_chr6.bim", header=F))
colnames(kgp3) <- c("chr","id","-","bp","a1","a2")
kgp3 <- kgp3[c("id","bp")]

missing <- data.frame(fread("missing_snps.txt", header=F))
colnames(missing) <- c("id","bp")

kgp <- rbind(kgp2,kgp1,kgp3,missing)
kgp <- kgp[!duplicated(kgp$id), ]
ordered_kgp <- kgp[order(kgp$id),]

merged <- merge(ordered_kgp, ordered_ref, by=c("id"))
merged[ , c(1,2,3)] <- merged[ , c(3,1,2)]
colnames(merged)[c(1,2,3)] <- colnames(merged)[c(3,1,2)]

#REORDER BY BASEPAIR COORDINATES
merged <- merged[order(merged$bp),]

#make reference panel file
ref_panel <- merged[ , !(names(merged) %in% "bp" )]
#write reference panel file
write.table(ref_panel, file = "ref_panel_phased.bgl", quote=FALSE, sep = "  ", row.names = FALSE)

#make and write markers file
alleles <- merged[ , !(names(merged) %in% c("I") )]
i=1
for (i in 1:nrow(alleles)){
  allele_list <- unique(c(alleles[i,]))
  allele_list <- unlist(allele_list)
  write(allele_list,file="markers_file.markers",append=TRUE,ncolumns = length(allele_list))
  i=i+1
}

library(xlsx)
library(data.table)
library(dplyr)
library(plyr)

ref <- data.frame(fread("Documents/Projects/MGH/imputation/ref_panel/MHC_haplotypes_CEU_HapMap3_ref_panel.bgl", header=T ))
ordered_ref <- ref[order(ref$id),]

kgp1 <- data.frame(fread("Documents/Projects/MGH/imputation/ref_panel/snp147_chr6_bp_only.txt", header=F))
colnames(kgp1) <- c("id","bp")

kgp2 <- data.frame(fread("Documents/Projects/MGH/imputation/ref_panel/1kgp3_chr6.bim", header=F))
colnames(kgp2) <- c("chr","id","-","bp","a1","a2")
kgp2 <- kgp2[c("id","bp")]

kgp3 <- data.frame(fread("Documents/Projects/MGH/imputation/ref_panel/1kgp1_chr6.bim", header=F))
colnames(kgp3) <- c("chr","id","-","bp","a1","a2")
kgp3 <- kgp3[c("id","bp")]

missing <- data.frame(fread("Documents/Projects/MGH/imputation/ref_panel/missing_snps.txt", header=F))
colnames(missing) <- c("id","bp")

kgp <- rbind(kgp2,kgp1,kgp3,missing)
kgp <- kgp[!duplicated(kgp$id), ]
ordered_kgp <- kgp[order(kgp$id),]

merged <- merge(ordered_kgp, ordered_ref, by=c("id"))
merged[ , c(1,2,3)] <- merged[ , c(3,1,2)]
colnames(merged)[c(1,2,3)] <- colnames(merged)[c(3,1,2)]

#make REF, ALT, and genotypes 
alleles <- subset(merged, select=-c(id, bp, I))
allele_list <- data.frame()
genotypes <- data.frame()
i=1
for (i in 1:nrow(alleles)){
  allele <- unique(c(alleles[i,]))
  allele_list[i,1] = allele[1]
  allele_list[i,2] = as.character(paste(unlist(allele), collapse = ","))
  
  gen_count=1
  for (j in seq(1,ncol(alleles),2)){
    allele_1 = match(alleles[i,j],allele)
    allele_2 = match(alleles[i,j+1],allele)
    genotypes[i,gen_count] <- as.character(paste(c(allele_1, allele_2), collapse = "|"))
    gen_count = gen_count+1
  }
  
}
colnames(allele_list) <- c('REF', 'ALT')
gen_names <- colnames(alleles)[seq(1,ncol(alleles),2)]

#make CHR 
chr = matrix(rep(6,each=nrow(merged)),nrow=nrow(merged))

#make POS 
pos <- merged['bp']

#make ID
id <- merged['id']

#make QUAL
qual <- matrix(rep('.',each=nrow(merged)),nrow=nrow(merged))

#make FILTER
filter <- matrix(rep('PASS',each=nrow(merged)),nrow=nrow(merged))

#make INFO
info <- matrix(rep('.',each=nrow(merged)),nrow=nrow(merged))

#make FORMAT
format <- matrix(rep('GT',each=nrow(merged)),nrow=nrow(merged))

#make vcf format
ref_panel = data.frame(chr,pos,id,allele_list,qual,filter,info,format,genotypes)
colnames(ref_panel) <- c('CHR','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',gen_names)
rownames(ref_panel) <- NULL
write(c('##fileformat=VCFv4.2'),file="ref_panel.vcf")
write.table(ref_panel,file="ref_panel.vcf", quote=FALSE, row.names = FALSE, sep='\t', append=TRUE)

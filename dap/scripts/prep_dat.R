#!/usr/bin/env Rscript
Sys.time()
library(optparse)
library(data.table)
library(stringr)
#supposed to be run in cis_dap folder
#mol, genotype and cov are all comming from the ../cis/data/folder
#make sure cov file has exactly the people want to include

#example usage prep_dat.R --mol_trait ENSG00000174595.4
#--pos ~/muscle_trans_snk/scripts/pos_strand.tab.gz

option_list <- list(
  make_option(c("-m", "--mol_trait"), type="character", default= NULL, 
              help="molecular trait for which dat file is"),
  
  make_option(c("-p", "--pos"),  type="character", default=NULL,
              help='position file of molecular traits must have start end')
 )
opt = parse_args(OptionParser(option_list=option_list))

pos = data.frame(fread(paste("zcat", opt$pos), stringsAsFactors=F,check.names=F))
pos = subset(pos,pos[,1]==opt$mol_trait)
gene_chr = pos[,'chr']
gene_start = pos[,'start']
gene_end = pos[,'end']
gene_coor = paste(gene_chr,gene_start,sep=':')
gene_coor = paste(gene_coor,gene_end,sep='-')
cmd = paste("tabix -h ../cis/data/moltraits.bed.gz", gene_coor)
exp = data.table(fread(cmd,sep="\t",header=T, stringsAsFactors=F,
                 check.names=FALSE))             
exp = data.frame(exp,check.names = F,stringsAsFactors=F)
exp = subset(exp,exp$pid == opt$mol_trait)

snp_chr = gene_chr
snp_start = gene_start-1000000
snp_end = gene_end+1000000
snp_coor = paste(snp_chr,snp_start,sep=':')
snp_coor = paste(snp_coor,snp_end,sep='-')

cmd = paste("tabix ../cis/data/genotypes.vcf.gz", snp_coor)
vcf = data.frame(fread(cmd,sep="\t",header=F, stringsAsFactors=F,
                              check.names=FALSE))

if (nrow(exp) >0 & nrow(vcf) >0) {
cmd = "tabix -H ../cis/data/genotypes.vcf.gz | tail -1"
vcf_header = data.table(fread(cmd,sep="\t",header=F, stringsAsFactors=F,
                              check.names=FALSE))             
vcf_header = as.character(vcf_header)

colnames(vcf) = vcf_header

cov = read.table('../cis/data/covariates.txt.gz',sep="\t",header=T, 
                 stringsAsFactors=F,check.names = F,row.names = 1)
sample_name = colnames(cov)

vcf_snp = vcf$ID
vcf = vcf[,colnames(cov)]
exp = exp[,colnames(cov)]

get_dose = function(x){
  x=t(x)
  x = data.frame(str_split(x, pattern=":", n = Inf, simplify = T),stringsAsFactors = F,check.names = F)
  return(x[,2])
}
vcf = as.data.frame(t(apply(vcf,1,get_dose)),stringsAsFactors = F,check.names = F)
colnames(vcf)=colnames(cov)
rownames(vcf)=vcf_snp

exp$type = 'pheno'
exp$trait = opt$mol_trait
exp$group = 'fusion'
vcf$type = 'geno'
vcf$trait = vcf_snp
vcf$group = 'fusion'
cov$type = 'controlled'
cov$trait = rownames(cov)
cov$group = 'fusion'

out = data.frame(rbind(exp,vcf,cov),stringsAsFactors = F, check.names = F)
outfile = paste0('sbams_data/',opt$mol_trait,'.dat')
write.table(out[,c('type','trait','group',sample_name)], 
            file = outfile, row.names=F, col.names=F, quote=F, sep = " ")
}else{
  cat('No Snps within 1Mb of ', opt$mol_trait, '\n')
}

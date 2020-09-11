#!/usr/bin/env Rscript
Sys.time()
library(optparse)
library(data.table)
library(gam)
library(ggplot2)
library(gridExtra)
library(grid)
library(stringr)

#opt$in_file='~/mus_exeRpt/exceRpt_outlier_eQTL/cis/qtltools_nominal.tsv.gz'
#opt$exp_level = '/net/snowwhite/home/FUSION/Tissue/miRNA/exceRpt/rpmmm/muscle_exceRpt_miRNA_RPMMM.txt.gz'
#opt$bin_file = '~/fusion_mus_adi_cmp/gam_muscle_mirna_same_denomitor_outlier_bin20_gam_pred_prop.tab.gz'
#opt$bin_number = '8,9,10,11,12,13,14,15'
#opt$begin_col=2
#opt$gene_name_col='miRNA'
#opt$threshold_list = '0.05,0.01'
#opt$outprx = "~/fusion_mus_adi_cmp/mus_mirna_outlier_bin20"

#get the lower boundary of the 1st bin to save
#remve them
#calculate 5% and 1% FDR
#see the number of sig genes

option_list <- list(
  make_option(c("-i", "--in_file"), type="character", default= NULL, 
              help="Input file: eQTL res"),
  
  make_option(c("-e", "--exp_level"), type="character", default= NULL, 
              help="Expression level file"),
    
  make_option(c("-p", "--bin_file"), type="character", default= NULL, 
              help="The txt output from the draw figure script"),
    
  make_option(c("-n", "--bin_number"), type="character", default= NULL, 
              help="a list of the first bins want to keep"),
  
  make_option(c("-b", "--begin_col"), type="numeric", default= NULL, 
              help="the column in exp_level file where the number begins"),
  
  make_option(c("-g", "--gene_name_col"), type="character", default= NULL, 
              help="the column in exp_level file where gene names are"),
  
  make_option(c("-t", "--threshold_list"), type="character", default= NULL, 
              help="the list of threshold that want to calculate res for"),

  make_option(c("-o", "--outprx"),  type="character", default=NULL,
              help='Prefix for the output file')
)

opt = parse_args(OptionParser(option_list=option_list))

sum_ori = data.frame(fread(paste("zcat", opt$bin_file), stringsAsFactors=F))

fdrs = as.numeric(str_split(opt$threshold_list, pattern=",", n = Inf, simplify = T)[1,])

report_ini = matrix(ncol=5)
report_ini = data.frame(report_ini,stringsAsFactors=F)
colnames(report_ini) = c('start_bin','fdr','pvalue','num_gene','num_pairs')

bins = as.numeric(str_split(opt$bin_number, pattern=",", n = Inf, simplify = T)[1,])

tryCatch({
  exp = data.frame(fread(paste("zcat", opt$exp_level), header=T, stringsAsFactors=F,check.names=F))
}, error=function(e){
  exp = data.frame(fread(opt$exp_level, header=T, stringsAsFactors=F,check.names=F))
})

mat = as.matrix(exp[,as.numeric(opt$begin_col):ncol(exp)])
rownames(mat) = exp[,as.character(opt$gene_name_col)]

mat = cbind('mean'=rowMeans(as.matrix(mat)),mat)
mat = data.frame(mat,stringsAsFactors=F)
mat = mat[order(mat$mean,decreasing = F),]

res = data.frame(fread(paste("zcat", opt$in_file), stringsAsFactors=F))

for (j in 1:length(bins)) {
 sum = subset(sum_ori, sum_ori$bin >= bins[j])
 mat_sub = subset(mat, mat$mean >= min(sum$mean))
 keep_genes = rownames(mat_sub)
 res_sub = subset(res, res$pheno_id %in% keep_genes)
 res_sub$BH_FDR = p.adjust(res_sub$p_nominal, method = 'BH')
 
 report = matrix(nrow=length(fdrs),ncol=5)
 report = data.frame(report,stringsAsFactors=F)
 colnames(report) = c('start_bin','fdr','pvalue','num_gene','num_pairs')
 
 report$start_bin = bins[j]
 report$fdr = fdrs
    for (i in 1:length(fdrs)) {
      res_fdr = subset(res_sub, res_sub$BH_FDR <= fdrs[i])
      report[i,'pvalue'] = max(res_fdr$p_nominal)
      report[i,'num_gene'] = length(unique(res_fdr$pheno_id))
      report[i,'num_pairs'] = nrow(res_fdr)  
}
     report_ini = rbind(report_ini, report) 
}

write.table(report_ini[-1,],file=paste0(opt$outprx,'_diff_fdr.tab'), row.names = F, col.names =T, quote=F, sep="\t")
Sys.time()

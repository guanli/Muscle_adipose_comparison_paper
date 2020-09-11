#!/usr/bin/env Rscript
Sys.time()
library(optparse)
library(data.table)
library(gam)
library(ggplot2)
library(gridExtra)
library(grid)

#opt$in_file='~/mus_exeRpt/exceRpt_outlier_eQTL/cis/qtltools_nominal.tsv.gz'
#opt$exp_level = '/net/snowwhite/home/FUSION/Tissue/miRNA/exceRpt/rpmmm/muscle_exceRpt_miRNA_RPMMM.txt.gz'
#opt$bin_number = 10
#opt$begin_col=2
#opt$gene_name_col='miRNA'
#opt$outprx = "~/fusion_mus_adi_cmp/gam_muscle_mirna_same_denomitor_outlier"

option_list <- list(
  make_option(c("-i", "--in_file"), type="character", default= NULL, 
              help="Input file: eQTL res"),
  
  make_option(c("-e", "--exp_level"), type="character", default= NULL, 
              help="Expression level file"),
  
  make_option(c("-n", "--bin_number"), type="numeric", default= NULL, 
              help="the number of bins wanted"),
  
  make_option(c("-b", "--begin_col"), type="numeric", default= NULL, 
              help="the column in exp_level file where the number begins"),
  
  make_option(c("-g", "--gene_name_col"), type="character", default= NULL, 
              help="the column in exp_level file where gene names are"),

  make_option(c("-o", "--outprx"),  type="character", default=NULL,
              help='Prefix for the output file')
)

opt = parse_args(OptionParser(option_list=option_list))

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
mat = subset(mat,mat$mean >0)

bin_size = ceiling(nrow(mat)/as.numeric(opt$bin_number))
mat$bin = head(rep(1:as.numeric(opt$bin_number),times=1,each=bin_size),nrow(mat))

res = data.frame(fread(paste("zcat", opt$in_file), stringsAsFactors=F))
res = cbind(res, 'bin'=mat[match(res$pheno_id,rownames(mat)),'bin'])
res = res[,c('pheno_id','var_id','p_nominal','bin')]
#p_nominal, phenotype_d,bin are needed
mat$bin = as.factor(mat$bin)
sig_list = list()
#for (i in 90:100){
for (i in 1:length(levels(mat$bin))){
  piece = subset(res, res$bin == levels(mat$bin)[i])
  piece$BH = p.adjust(piece$p_nominal, method = 'BH')
  if (min(piece$BH) <= 0.05) {
    sig_list[[i]] = unique(subset(piece, piece$BH <= 0.05)[,'pheno_id'])
    } else {
    sig_list[[i]] = NA
    } 
}

sig_list = unlist(sig_list)
sig_list = sig_list[!is.na(sig_list)]

mat$sig = ifelse(rownames(mat) %in% sig_list,1,0)
rm(res)

res = data.frame(cbind('gene'=rownames(mat),'sig'= mat$sig),stringsAsFactors=F)
rownames(res) = res$gene
res = cbind('mean'=mat[match(rownames(res),rownames(mat)),'mean'], res)
res = subset(res, res$mean > 0)
res$log10mean=log10(res$mean)
res = res[order(res$log10mean,decreasing = F),]

res$bin = head(rep(1:as.numeric(opt$bin_number),times=1,each=bin_size),nrow(res))
res = res[,-c(match('gene',colnames(res)))]

res$sig = as.numeric(res$sig)
res$log10mean = as.numeric(res$log10mean)

model = paste0("sig"," ~ s(log10mean)")
fit = gam(as.formula(model), data = res, family = "binomial")
res = cbind(res, fitted(fit))
colnames(res)[ncol(res)]=paste0('gam_pred_',"sig")

#draw a fig for each trait
dat = aggregate(res[,"sig"], by=list(bin=res$bin), FUN=mean)
dat[,1]=as.numeric(dat[,1])
line_y = paste0('gam_pred_',"sig")

png(paste0(opt$outprx,'_',"sig",'.png'))
p1=ggplot(dat, aes(x=bin, y=x))+
  geom_bar(stat="identity") +
  ylab("Proportion of significant molecular traits") 
p2=ggplot(data=res, aes(log10mean, get(line_y)))+
  geom_line(color='red')+
  ylab("Predicted probability of significant molecular traits")
grid.arrange(p1,p2,nrow=1)
dev.off()

gzfh = gzfile(paste0(opt$outprx,'_gam_pred_prop.tab.gz'), "w")
write.table(cbind('molecular_traits'=rownames(res),res), gzfh, row.names=F, col.names=T, quote=F, sep = "\t")
close(gzfh)



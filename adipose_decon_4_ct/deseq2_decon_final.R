#!/usr/bin/env Rscript
#.libPaths("/net/snowwhite/home/guanli/transFusion/packrat/lib/x86_64-pc-linux-gnu/3.4.4")
Sys.time()

#example useage:
#Rscript deseq2_decon.R -r ~/AdiposeProject/tissDecon/cibersort/kerrin_4_type_blood_mo/ref.txt
                      # -m /net/snowwhite/home/FUSION/Tissue/freeze5/adipose/gene/freeze5.adipose.analysis.tpm.all_genes.tsv.gz
                      # -c ~/AdiposeProject/tissDecon/cibersort/kerrin_4_type_blood_mo/class.txt
                      # -v ~/AdiposeProject/tissDecon/unmix/july62018/shift_value.txt
                      # -o /net/snowwhite/home/guanli/AdiposeProject/tissDecon/unmix/july62018/

library(optparse)
library(data.table)
library(DESeq2)
#library(pbapply)
library(readr)
library(data.table)
library(vsn)
library(Biobase)
library(dplyr)

option_list <- list(
    make_option(c("-r", "--reference"), type="character", action="store", default= NULL, metavar="character",
    help="Reference file"),
    
    make_option(c("-m", "--mixture"), type="character", action="store", default= NULL, metavar="character",
    help="Mixture file"),

    make_option(c("-c", "--class"), type="character", action="store", default= NULL, metavar="character",
    help="Classification file for reference"),
    
    make_option(c("-v", "--value"), type="character", action="store", default= NULL, metavar="character",
    help="Shift value file"),
    
    make_option(c("-o", "--output_dir"), type="character", default= NULL, 
    help="output directory")
)

opt = parse_args(OptionParser(option_list=option_list))

adi = data.frame(fread(paste("zcat",opt$mixture), sep="\t", header=T, stringsAsFactors=F))
rownames(adi) = adi[,1]
adi = adi[,-1]

ref = read.table(opt$reference, sep="\t", header=T, stringsAsFactors=F)
rownames(ref) = ref[,1]
ref = ref[,-1]
cat('Head of reference file is: \n')
print(head(ref))

class = read.table(opt$class, sep="\t", header=F, stringsAsFactors=F)
cat('Classification of reference file is: \n')
print(class)

rownames(class)= class[,1]
class=class[,-1]
class=as.matrix(t(class))
class[class==2] = 0
idx =factor((class %*% (1:ncol(class))),labels=colnames(class))

deseq_ref <- data.frame(cbind(rowMedians(as.matrix(ref[,idx==colnames(class)[1]]))))
for (i in(2:ncol(class))) {
  deseq_ref <- data.frame(cbind(deseq_ref, rowMedians(as.matrix(ref[,idx==colnames(class)[i]]))))
}
colnames(deseq_ref) = colnames(class)                      
rownames(deseq_ref) = rownames(ref)
cat('Head of deseq_ref file is: \n')
print(head(deseq_ref))

deseq_ref$max = rowMaxs(as.matrix(deseq_ref))
deseq_ref = subset(deseq_ref, deseq_ref$max > 0)
common <- intersect(rownames(deseq_ref),rownames(adi)) #24136 genes left
adi = adi[common,]
deseq_ref = deseq_ref[common,1:ncol(class)]

deseq_ref <- sweep(deseq_ref, 2, 1e6/colSums(deseq_ref), "*")  
adi <- sweep(adi, 2, 1e6/colSums(adi), "*")  

dir = opt$output_dir

#study_sample is a matrix with all numbers
test_shift<-function(shift, study_sample, ...){
  png(paste0(dir, "meanSdPlot_",shift,".png"))
  meanSdPlot(log(study_sample + shift))
  dev.off()
  return(meanSdPlot(log(study_sample + shift), plot=F))
  
}

vals = read.table(opt$value, sep="\t", header=F, stringsAsFactors=F)$V1

cv<-vector()
for (i in 1:length(vals)) {
  temp<-test_shift(vals[i], as.matrix(adi[,-1]))
  #Take the points on the red line of meanSdPlot, calculate SD and mean, summarize with CV = SD/mean
  tempsd<-temp$sd
  cv[i]<-sd(tempsd)/mean(tempsd)
}	
plotds<-cbind.data.frame(vals,cv)
write.table(plotds,file=paste0(dir,"test_shift_cvs.txt"),row.names = F, col.names = T, quote=F, sep="\t")


vals<-c(plotds[which.min(plotds[,2]),]$vals)
for (i in 1:length(vals)) {
  #  unmix takes about 11 minutes on snowwhite.
  mix <- unmix(as.matrix(adi[,-1]), pure=as.matrix(deseq_ref), shift=vals[i])
  rownames(mix) = colnames(adi)[-1]
  write.table(cbind('sample'=rownames(mix),mix),paste0(dir,"DESeq2_tissue_proportions_shift_",vals[i],".txt"), col.names=T, row.names=F, quote=F, sep="\t")
}

gzfh = gzfile(paste0(dir,"unmix_pure.tab.gz"),'w')
write.table(cbind('gene'=rownames(deseq_ref),deseq_ref), gzfh,row.names = F, col.names = T, quote=F, sep="\t")
close(gzfh)

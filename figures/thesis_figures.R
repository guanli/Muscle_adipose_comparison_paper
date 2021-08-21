library(ggplot2)
library(scales)
library(gridExtra)
library(grid)
library(ggrepel)
library(scales)
library(ggpubr) 
library(Hmisc)
library(data.table)
library(RColorBrewer)
library(stringr)
library(cowplot)
library(R.utils)
library(ggforestplot)
library(GGally)
#in ~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011
text_size =  20#if pdf, do text_size = 20 or 25, if tiff 7
text_font = 'bold'
col_palette = 'Dark2'
go_axis_size = 6

ggcust <- function(...){
  ggplot(...) +
    theme_bw()+
    theme(axis.title = element_text(size=text_size),
          axis.text=element_text(size=text_size),
          legend.text=element_text(size=text_size),
          strip.text=element_text(size=text_size),
          #legend.title = element_text(size=text_size),
          legend.title =element_blank(),
          plot.title = element_text(hjust = 0.5,size=text_size))
}

####Figure 2.3.1 mRNA and miRNA cis-QTL discovery
################################################################################
load('~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/keep_traits.RData')
f = '~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/muscle_mrna_rc_gam_pred_prop.tab.gz'
mus_mrna = read.table(f,header = T, sep = '\t', stringsAsFactors = F)
mus_mrna = subset(mus_mrna, mus_mrna$molecular_traits %in% keep_traits$mus_mrna_qtl)
mus_mrna$Tissue = 'Muscle'
mus_mrna$Type = 'mRNA'
f = '~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/adipose_mrna_rc_gam_pred_prop.tab.gz'
adi_mrna = read.table(f,header = T, sep = '\t', stringsAsFactors = F)
adi_mrna = subset(adi_mrna, adi_mrna$molecular_traits %in% keep_traits$adi_mrna_qtl)
adi_mrna$Tissue = 'Adipose'
adi_mrna$Type = 'mRNA'
f = '~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/muscle_mirna_rc_gam_pred_prop.tab.gz'
mus_mirna = read.table(f,header = T, sep = '\t', stringsAsFactors = F)
mus_mirna = subset(mus_mirna, mus_mirna$molecular_traits %in% keep_traits$mus_mirna_qtl)
mus_mirna$Tissue = 'Muscle'
mus_mirna$Type = 'miRNA'
f = '~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/adipose_mirna_rc_gam_pred_prop.tab.gz'
adi_mirna = read.table(f,header = T, sep = '\t', stringsAsFactors = F)
adi_mirna = subset(adi_mirna, adi_mirna$molecular_traits %in% keep_traits$adi_mirna_qtl)
adi_mirna$Tissue = 'Adipose'
adi_mirna$Type = 'miRNA'

df = data.frame(rbind(mus_mrna,adi_mrna,mus_mirna,adi_mirna),stringsAsFactors = F)
df$Type = factor(df$Type,levels = c('mRNA','miRNA'))
df$Tissue = factor(df$Tissue,levels = c('Muscle','Adipose'))

p1=ggcust(data=df, aes(x=log10mean,y=gam_pred_sig,color=Type))+
  geom_line(aes(linetype=Tissue, color=Type),size=1.5)+
  labs(x = "log10 mean read count",
       y = "Probability of having QTL")+
  scale_color_brewer(palette = col_palette,drop=FALSE)+
  theme(legend.position = c(0.1, 0.8),
        legend.key = element_blank(),
        legend.key.size = unit(0.5,'in'))+
  scale_y_continuous(breaks=seq(0.0, 1, 0.1), limits=c(0, 0.5))

f='~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/df_new.tab.gz'
df_new= data.frame(fread(f, sep="\t", header=T, stringsAsFactors=F)) 
brks<-c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5)
labs=c("0.02","0.04","0.06","0.08","0.1","0.2","0.3","0.4","0.5")
tmp = subset(df_new,Tissue == 'Muscle' & Type == 'mRNA')
p4=ggcust(tmp,aes(x=maf)) + 
  geom_histogram(data=subset(tmp, snp_type == 'All'),aes(y=..count../sum(..count..) * 100),
                 fill = "#666666", alpha = 0.8, binwidth=0.01, boundary = 0.02) +
  geom_histogram(data=subset(tmp, snp_type == 'Lead'),aes(y=..count../sum(..count..) * 100),
                 fill = "#1B9E77", alpha = 0.5, binwidth=0.01, boundary = 0.02) +
  scale_x_continuous(breaks=brks, labels=labs)+
  labs(x="Minor allele frequence", y="Proportions",title = "mRNA lead variants vs all tested varints") +
  scale_y_continuous(breaks = seq(0,20,2),limits = c(0,20))+
  theme(axis.text.x=element_text(angle=-90, hjust=0))

tmp = subset(df_new,Tissue == 'Muscle' & Type == 'miRNA')
p5=ggcust(tmp,aes(x=maf)) + 
  geom_histogram(data=subset(tmp, snp_type == 'All'),aes(y=..count../sum(..count..) * 100),
                 fill = "#666666", alpha = 0.8, binwidth=0.01, boundary = 0.02) +
  geom_histogram(data=subset(tmp, snp_type == 'Lead'),aes(y=..count../sum(..count..) * 100),
                 fill = "#D95F02", alpha = 0.5, binwidth=0.01, boundary = 0.02) +
  scale_x_continuous(breaks=brks, labels=labs)+
  labs(x="Minor allele frequence", y="Proportions",title = "miRNA lead variants vs all tested varints") +
  scale_y_continuous(breaks = seq(0,20,2),limits = c(0,20))+
  theme(axis.text.x=element_text(angle=-90, hjust=0))

tmp = subset(df_new,Tissue == 'Muscle' & snp_type == 'Lead')
p6=ggcust(tmp,aes(x=maf)) + 
  geom_histogram(data=subset(tmp, Type == 'mRNA'),aes(y=..count../sum(..count..) * 100),
                 fill = "#1B9E77", alpha = 0.8, binwidth=0.01, boundary = 0.02) +
  geom_histogram(data=subset(tmp, Type == 'miRNA'),aes(y=..count../sum(..count..) * 100),
                 fill = "#D95F02", alpha = 0.5, binwidth=0.01, boundary = 0.02) +
  scale_x_continuous(breaks=brks, labels=labs)+
  labs(x="Minor allele frequence", y="Proportions",title = "mRNA lead variants vs miRNA lead variants") +
  scale_y_continuous(breaks = seq(0,20,2),limits = c(0,20))+
  theme(axis.text.x=element_text(angle=-90, hjust=0))

pdf('figs/fig_1_mrna_mirna_qtl_4figversion.pdf')
grid.arrange(p1,p4,p5,p6, ncol = 2) 
dev.off()


####Figure 2.3.2 Probability of having eQTL with confidence interval
################################################################################
load('~/fusion_mus_adi_cmp/figs/csg_CI_curve.RData')

df$Type = factor(df$Type,levels = c('mRNA','miRNA'))
df$Tissue = factor(df$Tissue,levels = c('Muscle','Adipose'))

pdf('~/fusion_mus_adi_cmp/figs/csg_curve_ci_mus_2.pdf')
ggcust(data=subset(df,Tissue=='Muscle'), aes(x=log10mean,y=gam_pred_sig))+
  geom_line(aes(color=Type),size=1.5)+
  geom_ribbon(data=subset(df,Tissue=='Muscle' & Type=='mRNA'),
              aes(ymin=lower, ymax=upper), alpha=0.2, linetype = 0)+
  geom_ribbon(data=subset(df,Tissue=='Muscle' & Type=='miRNA'),
              aes(ymin=lower, ymax=upper), alpha=0.2, linetype = 0)+
  labs(x = "log10 mean read count",
       y = "Probability of having QTL")+
  scale_color_brewer(palette = col_palette,drop=FALSE)+
  theme(legend.position = c(0.15, 0.8),legend.key = element_blank(),legend.key.size = unit(0.3,'in'))+
  scale_y_continuous(breaks=seq(-0.1, 1, 0.1), limits=c(-0.1, 0.6))
dev.off()

pdf('~/fusion_mus_adi_cmp/figs/csg_curve_ci_adi.pdf')
ggcust(data=subset(df,Tissue=='Adipose'), aes(x=log10mean,y=gam_pred_sig,color=Type))+
  geom_line(aes(color=Type),size=1.5)+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, colour = NA)+
  labs(x = "log10 mean read count",
       y = "Probability of having QTL")+
  scale_color_brewer(palette = col_palette,drop=FALSE)+
  theme(legend.position = c(0.1, 0.8),legend.key = element_blank(),legend.key.size = unit(0.3,'in'))+
  scale_y_continuous(breaks=seq(-0.1, 1, 0.1), limits=c(-0.1, 0.6))+
  scale_x_continuous(breaks=seq(0, 6, 2), limits=c(min(subset(df,Tissue=='Adipose')[,'log10mean']), 6))
dev.off()


####Figure 2.3.3 Multiple independent QTL discovery (DAP results)
################################################################################
dap = data.frame(cbind('number'=c(1754,5860,2408,558,123,23,8,1,1,
                                  23,94,7,1,0,0,0,0,0,
                                  24910,98160,21004,3189,512,90,21,3,0,
                                  1901,6596,2673,688,149,46,10,4,1,
                                  36,101,18,3,1,0,0,0,0,
                                  18727,88350,15664,1995,313,54,14,5,0),
                       'category'=rep(seq(0,8),6),
                       'type'=c(rep('mRNA',9),rep('miRNA',9),rep('DNAme',9),
                                rep('mRNA',9),rep('miRNA',9),rep('DNAme',9)),
                       'tissue'=c(rep('Muscle',27),rep('Adipose',27))),stringsAsFactors = F)
dap$number = as.numeric(dap$number)
dap$category = factor(dap$category)
dap$tissue = factor(dap$tissue,levels=c('Muscle','Adipose'))
dap$type = factor(dap$type, levels = c('mRNA','miRNA','DNAme'))
pdf('~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/figs/DAP.pdf',height=cm(2),width = cm(4.5))
ggcust(dap, aes(x=category, y=number, fill=tissue)) +
  geom_bar(stat='identity', position=position_dodge(preserve = "single")) + 
  geom_text(aes(label=number), position=position_dodge(width=0.9), vjust=-0.25,angle=90)+
  facet_wrap(~ type,scales="free_y")+
  labs(x="Number of independent QTLs", 
       y="Counts of genes/DNAme sites") + 
  theme(legend.position="bottom") + 
  scale_fill_brewer(palette = col_palette,drop=FALSE) 
dev.off()


####Figure 2.3.11 luciferase assay plots
################################################################################
luci=data.frame('y'=c(2.5680, 2.9867, 3.5399, 4.4935, 4.1586,
                      1.9092, 1.7860, 1.6487, 2.2033, 2.1614,
                      2.3551, 2.6612, 2.9240, 3.3962, 3.0681,
                      2.0438, 2.4145, 1.9438, 1.7055, 1.8367,
                      20.4710, 8.2673, 4.1485, 8.8494, 9.3987,
                      1.9989, 7.7128, 4.1040, 3.3006, 4.5003,
                      10.4935,17.4137,NA, 8.34,9.4343,
                      5.7740,5.1348,6.9734,1.9299,1.7775,
                      1.1052,0.6968,1.0832,1.1148,0.7446,2.0147,1.1193,1.1361),
                'Orientation'=c(rep('Forward',10),rep('Reverse',10),rep('Forward',10),rep('Reverse',10),rep('EV',8)),
                'Cell_types'=c(rep('Preadipocyte',20),rep('Adipocyte',20),rep('Preadipocyte',4),rep('Adipocyte',4)),
                'Allele'=c(rep('G',5),rep('C',5),rep('G',5),rep('C',5),
                           rep('G',5),rep('C',5),rep('G',5),rep('C',5),rep('EV',8)))

luci$Allele = factor(luci$Allele,levels = c('EV','C','G'))
luci$Cell_types = factor(luci$Cell_types,levels = c('Preadipocyte','Adipocyte'))
luci$Orientation= factor(luci$Orientation,levels = c('EV','Forward','Reverse'))

pdf('figs/sub_Preadipocyte_luciferase.pdf',height=cm(2.5),width = cm(4.5))
ggcust(subset(luci,Cell_types=='Preadipocyte'), aes(x=Allele, y=y)) + 
  geom_dotplot(aes(fill=Allele),drop = TRUE, color='white',
               binaxis='y',stackdir='center')+
  facet_wrap(~ Orientation, scales = "free",nrow=1)+
  scale_y_continuous(limits = c(0,5),breaks = seq(0,5,1))+
  labs(x = "rs11688682",y = "Relative luciferase activity",title='Preadipocyte')+
  scale_fill_manual(values = c("#666666","#1B9E77","#D95F02"),labels = c("Empty vector", "C", "G"))+
  theme(strip.text.x = element_text(size = 20),
        plot.title = element_text(size = 20))
dev.off()
pdf('figs/sub_adipocyte_luciferase.pdf',height=cm(2.5),width = cm(4.5))
ggcust(subset(luci,Cell_types=='Adipocyte'), aes(x=Allele, y=y)) + 
  geom_dotplot(aes(fill=Allele),drop = TRUE, color='white',
               binaxis='y',stackdir='center')+
  facet_wrap(~ Orientation, scales = "free",nrow=1)+
  scale_y_continuous(limits = c(0,21),breaks = seq(0,21,5))+
  labs(x = "rs11688682",y = "Relative luciferase activity",title='Preadipocyte')+
  scale_fill_manual(values = c("#666666","#1B9E77","#D95F02"),labels = c("Empty vector", "C", "G"))+
  theme(strip.text.x = element_text(size = 20),
        plot.title = element_text(size = 20))
dev.off()


####Figure 2.3.12 Percent of mRNAs/miRNAs/DNAme sites associated with physiological trait: muscle-adipose comparison
################################################################################
load('~/fusion_mus_adi_cmp/trait_asso/res/res_fdr01_bon.RData')
res$pct = (res$number/res$totmolec)*100
res=subset(res,res$thresh=='FDR01')
res=subset(res, model %in% c("TissueFiber","Component_17"))

#48 physiological traits
int=c("bmi","DIo","fS_C_pept","fS_C_pept_30","fS_Kol_HDL","fS_Trigly",                    
      "GL120","GL30","GL60","glu_2h_biopsy","Glu_AUC_0to30","HIP" ,"HOMA",
      "Ins_AUC_0to30","InsGenIn","InsSec30","matsuda_3pt","matsuda_4pt",                  
      "p_insu","p_insu_120","p_insu_30","p_insu_60","RFM","S_ALAT",                       
      "S_GT","S_hs_CRP","S_Insu","S_Insu_30","S_LipoA1","S_Uraat",                    
      "T2D_NGT","WAIST","WEIGHT","whr","GL0","glu_fast_biopsy","S_LipoB",
      "ApoB_A1_ratio","whradjbmi","B_GHb_A1C","fS_Kol_LDL_c","B_HbA1c",                     
      "sbp","fS_Krea","HEIGHT","fS_Kol","dbp","CpepGenIn")   

res=subset(res,trait %in% int)

res$trait_long[res$trait=='HOMA'] = 'HOMA'
res$trait_long[res$trait=='T2D_NGT'] = 'T2D'
res$trait_long[res$trait=='WAIST'] = 'Waist'
res$trait_long[res$trait=='WEIGHT'] = 'Weight'
res$trait_long[res$trait=='bmi'] = 'BMI'
res$trait_long[res$trait=='DIo'] = 'Disposition index'
res$trait_long[res$trait=='fS_C_pept'] = 'Fasting serum C peptide'
res$trait_long[res$trait=='fS_C_pept_30'] = 'Fasting serum C peptide 30min'
res$trait_long[res$trait=='fS_Kol_HDL'] = 'HDL cholesterol'
res$trait_long[res$trait=='fS_Trigly'] = 'Triglycerides'
res$trait_long[res$trait=='GL120'] = '2h glucose from OGTT'
res$trait_long[res$trait=='GL30'] = '30min glucose from OGTT'
res$trait_long[res$trait=='GL60'] = '60min glucose from OGTT'
res$trait_long[res$trait=='GL0'] = 'Fasting glucose from OGTT'
res$trait_long[res$trait=='glu_2h_biopsy'] = '2h glucose closest to biopsy date'
res$trait_long[res$trait=='Glu_AUC_0to30'] = 'Glucose AUC 0 to 30min'
res$trait_long[res$trait=='HIP'] = 'Hip circumference'
res$trait_long[res$trait=='Ins_AUC_0to30'] = 'Insulin AUC 0 to 30min'
res$trait_long[res$trait=='InsGenIn'] = 'Insulinogenic index'
res$trait_long[res$trait=='InsSec30'] = 'InsulinAUC glucoseAUC ratio 0 to 30min'
res$trait_long[res$trait=='matsuda_3pt'] = 'Matsuda index 0, 30, 120min'
res$trait_long[res$trait=='matsuda_4pt'] = 'Matsuda index 0, 30, 60, 120min'
res$trait_long[res$trait=='p_insu'] = 'Fasting plasma insulin'
res$trait_long[res$trait=='p_insu_30'] = '30min plasma insulin from OGTT'
res$trait_long[res$trait=='p_insu_60'] = '60min plasma insulin from OGTT'
res$trait_long[res$trait=='p_insu_120'] = '120min plasma insulin from OGTT'
res$trait_long[res$trait=='RFM'] = 'Relative fat mass'
res$trait_long[res$trait=='S_ALAT'] = 'Alanine aminotransferase(ALT)'
res$trait_long[res$trait=='S_GT'] = 'Glutamyltransferase(GGT)'
res$trait_long[res$trait=='S_hs_CRP'] = 'C-Reactive protein'
res$trait_long[res$trait=='S_Insu'] = 'Fasting serum insulin'
res$trait_long[res$trait=='S_Insu_30'] = '30min serum insulin'
res$trait_long[res$trait=='S_LipoA1'] = 'Apolipoprotein A1(A1)'
res$trait_long[res$trait=='S_Uraat'] = 'Serum uric acid'
res$trait_long[res$trait=='whr'] = 'Waist hip ratio'
res$trait_long[res$trait=='glu_fast_biopsy'] = 'Fasting glucose'
res$trait_long[res$trait=='S_LipoB'] = 'Apolipoprotein B(ApoB)'
res$trait_long[res$trait=='ApoB_A1_ratio'] = 'ApoB A1 ratio'
res$trait_long[res$trait=='whradjbmi'] = 'BMI adjusted WHR'
res$trait_long[res$trait=='B_GHb_A1C'] = 'glycated Hemoglobin A1c'
res$trait_long[res$trait=='fS_Kol_LDL_c'] = 'LDL cholesterol'
res$trait_long[res$trait=='B_HbA1c'] = 'Hemoglobin A1c'
res$trait_long[res$trait=='sbp'] = 'Systolic blood pressure'
res$trait_long[res$trait=='CpepGenIn'] = 'C peptidogenic index'
res$trait_long[res$trait=='dbp'] = 'Diastolic blood pressure'
res$trait_long[res$trait=='HEIGHT'] = 'Height'
res$trait_long[res$trait=='fS_Krea'] = 'Creatinine'
res$trait_long[res$trait=='fS_Kol'] = 'Total cholesterol'

res$trait_long = factor(res$trait_long,levels = c("T2D","Matsuda index 0, 30, 120min",
                                                  "Matsuda index 0, 30, 60, 120min","HOMA","Fasting serum insulin",
                                                  "30min serum insulin","Insulin AUC 0 to 30min",
                                                  "Fasting serum C peptide","Fasting serum C peptide 30min","60min plasma insulin from OGTT",  
                                                  "Fasting plasma insulin","120min plasma insulin from OGTT","30min plasma insulin from OGTT",
                                                  "Insulinogenic index","InsulinAUC glucoseAUC ratio 0 to 30min","Disposition index",'C peptidogenic index',
                                                  "glycated Hemoglobin A1c","Hemoglobin A1c","Fasting glucose","Fasting glucose from OGTT","2h glucose from OGTT","2h glucose closest to biopsy date", 
                                                  "30min glucose from OGTT","60min glucose from OGTT","Glucose AUC 0 to 30min",
                                                  "BMI","Relative fat mass","Waist hip ratio","Hip circumference","Waist","Weight","BMI adjusted WHR",
                                                  "Triglycerides","HDL cholesterol","LDL cholesterol","Total cholesterol","ApoB A1 ratio","Apolipoprotein A1(A1)",
                                                  "Apolipoprotein B(ApoB)","Glutamyltransferase(GGT)","Alanine aminotransferase(ALT)","C-Reactive protein",
                                                  "Serum uric acid",'Creatinine',"Systolic blood pressure","Diastolic blood pressure",'Height'))


dat = res

pdf('figs/muscle_adi_final_model.pdf',width = cm(5),height = cm(4.5))
ggplot(dat,aes(x=trait_long, y=pct, fill=Tissue,color = Tissue)) +
  geom_bar(stat='identity',position=position_dodge(preserve = "single"))+ 
  facet_wrap(~ type,ncol = 1,scales="free_y",drop = F)+
  labs(x="",y=expression("Percent of tests (FDR"<="1%)")) + 
  scale_fill_brewer(palette = "Dark2",drop=FALSE,name = "Tissue") +
  scale_color_brewer(palette = "Dark2",drop=FALSE,name = "Tissue") +
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text=element_text(size=15),
        legend.text=element_text(size=20),
        strip.text=element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(hjust = 0.5,size=20),
        axis.text.x=element_text(angle=-90, hjust=0),
        legend.position="bottom")
dev.off()


####Figure 2.3.13-14 muscle adipose comparison scatterplot for insulin and BMI
################################################################################
##same figures with different data
##first merge muscle and adipose results for the same type of molecular trait together 
##then draw plots
#f="/net/snowwhite/home/aujackso/Tissue_Li/Trait_mRNA.UPDATEDbase/mus.BASE.BASETISSUE/results/tissuefiber_qts.bmi.results.tab"
#f='/net/snowwhite/home/aujackso/Tissue_Li/Trait_miRNA/muscle/results/tissuefiber_qts.bmi.results.tab'
#f='/net/snowwhite/home/aujackso/Tissue_Li/Trait_DNAme/muscle/results/tissuefiber_qts.S_Insu.results.tab.gz'
f='/net/snowwhite/home/aujackso/Tissue_Li/Trait_DNAme/muscle/results/tissuefiber_qts.bmi.results.tab.gz'
mus=data.frame(fread(f, sep="\t", header=T, stringsAsFactors=F,check.names = F),check.names = F)
colnames(mus)[-1]=paste0('mus_',colnames(mus)[-1])
#f='/net/snowwhite/home/aujackso/Tissue_Li/Trait_mRNA.UPDATEDbase/adi.BASEPAIVI/results/paivi_qts.bmi.results.tab'
#f='/net/snowwhite/home/aujackso/Tissue_Li/Trait_miRNA/adipose.PAIVI/results/paivi_qts.bmi.results.tab'
#f='/net/snowwhite/home/aujackso/Tissue_Li/Trait_DNAme/adipose.PAIVI/results/paivi_qts.S_Insu.results.tab'
f='/net/snowwhite/home/aujackso/Tissue_Li/Trait_DNAme/adipose.PAIVI/results/paivi_qts.bmi.results.tab'
adi=data.frame(fread(f, sep="\t", header=T, stringsAsFactors=F,check.names = F),check.names = F)
colnames(adi)[-1]=paste0('adi_',colnames(adi)[-1])

df=merge(mus,adi,by='response')
#bon=0.05/30483 
#bon=0.05/974
bon=0.05/694089
##30483, 974 and 694089 are the total number of mRNAs, miRNAs and DNAme sites tested in two tissues together

df$mus_sig = ifelse(df$mus_p.value<=bon,1,0)
df$adi_sig = ifelse(df$adi_p.value<=bon,1,0)

df$sign_mlog10p_mus=ifelse(df$mus_estimate >0, -log10(df$mus_p.value),log10(df$mus_p.value))
df$sign_mlog10p_adi=ifelse(df$adi_estimate >0, -log10(df$adi_p.value),log10(df$adi_p.value))
df$group[df$mus_sig == 1 & df$adi_sig == 1] = 4
df$group[df$mus_sig == 0 & df$adi_sig == 1] = 2
df$group[df$mus_sig == 1 & df$adi_sig == 0] = 1
df$group[df$mus_sig == 0 & df$adi_sig == 0] = 3
df$group = factor(df$group)
lim_value = ceiling(max(abs(df$sign_mlog10p_mus),abs(df$sign_mlog10p_adi)))

#png('figs/sub_mus_adi_bmi_mrna.png')
#png('figs/sub_mus_adi_insu_mirna.png')
#png('figs/sub_mus_adi_bmi_mirna.png')
#png('figs/sub_mus_adi_insu_dname.png')
png('figs/sub_mus_adi_bmi_dname.png')
ggcust(df,aes(x=sign_mlog10p_mus, y=sign_mlog10p_adi,colour = group)) +
  geom_point(alpha=1,size=2)+ 
  geom_abline(intercept=0, slope = 1, color = "blue",size=0.5,alpha=0.5)+
  geom_abline(intercept=0, slope = -1, color = "blue",size=0.5,alpha=0.5)+
  geom_hline(yintercept=0, color = "blue",alpha=0.5) +
  geom_vline(xintercept=0, color = "blue",alpha=0.5) +
  labs(x="Muscle signed -log10(P)", 
       y="Adipose signed -log10(P)") + 
  xlim(-lim_value,lim_value) +
  ylim(-lim_value,lim_value) +
  scale_color_manual(values = c("#1B9E77","#7570B3","#d2cfd4","#D95F02"),
                     breaks=c("1", "2", "3","4"),
                     name=expression("FDR"<="1%"),
                     labels=c("Muscle", "Adipose","Neither","Both"))
dev.off()


####Figure 2.3.15 EIF4EBP1 ~ physiological traits
####Figure 2.3.16 INHBB ~ physiological traits
################################################################################
#As these two figures are the same figures with different data, the codes for 
#the two figures are the same
#data file
#~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/figs/EIF4EPB1_phy.RData
#~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/figs/INHBB_phy.RData

tmp=subset(res,Tissue=='Adipose') #INHBB order by adi
#tmp=subset(res,Tissue=='Muscle') #EIF order by muscle
tmp=tmp[order(tmp$p.value,decreasing = F),]
level_order=tmp$trait_long
res$trait_long = factor(res$trait_long,levels = level_order)

#pdf('figs/mus_adi_mrna_INHBB_1.pdf',height=cm(3),width = cm(6))
pdf('figs/mus_adi_mrna_EIF4EPB1_1.pdf',height=cm(3),width = cm(6))
ggcust(data=res, 
       aes(x = trait_long, y = estimate, ymin = CI_low, ymax = CI_high, colour=Tissue)) +
  geom_hline(yintercept = 0, colour = 'black', alpha = 0.5) +
  geom_point(position = position_dodge(width=0.5)) +
  geom_errorbar(width=0,position = position_dodge(width=0.5))  +
  #coord_flip() +
  xlab(NULL) + ylab("Regression coefficient [95% Confidence interval]") +
  scale_color_brewer(palette = col_palette,drop=FALSE)+
  theme(axis.text.x=element_text(angle=-90, hjust=0),
        legend.position = "bottom")
dev.off() 


#Codes to generate boxplot like figure 2.3.9 A
#rs516946_cg12439423
#Rscript ~/fusion_mus_adi_cmp/scripts/gene_snp_plot_gt_label.R -t ~/fusion_mus_adi_cmp/figs/paper_figs/QTM/rs516946_cg12439423_mus -y cg12439423 -b ~/muscle_dname/eQTL_201906/cis/data/moltraits.bed.gz -v 8:41519248  -p ~/muscle_trans_snk/data_201906/genotypes.vcf.gz -c ~/muscle_dname/eQTL_201906/cis/data/covariates.txt.gz -q 4 -w 4
#-t is the output prefix
#-y is the name of the quantitative trait
#-b is the expression level file of the quantitative trait
#-v is the SNP
#-p is the genotype file
#-c is the covariate file
#-q and -w is the length and witch of the output figure

#Codes to generate boxplot like figure 2.3.9 D
#INHBB  cg14231073
#Rscript ~/fusion_mus_adi_cmp/figs/paper_figs/scripts/QTM_plot.R -t ~/fusion_mus_adi_cmp/figs/paper_figs/QTM/QTM_cg14231073_INHBB_mus -x cg14231073 -a ~/muscle_dname/eQTL_201906/cis/data/moltraits.bed.gz -y ENSG00000163083.5 -b ~/muscle_trans_snk/eQTL_201906/cis/data/moltraits.bed.gz -v 2:121347612 -p ~/muscle_trans_snk/data_201906/genotypes.vcf.gz -c ~/fusion_mus_adi_cmp/figs/paper_figs/QTM/data/qtm_mus_covariates.txt.gz -q 8 -w 8
#-t is the output prefix
#-x is one of the two quantitative traits
#-a is the expression level file of quantitative  trait X
#-y is the other one of the two quantitative traits
#-b is the expression level file of quantitative  trait b
#-v is the SNP
#-p is the genotype file
#-c is the covariate file
#-q and -w is the length and witch of the output figure

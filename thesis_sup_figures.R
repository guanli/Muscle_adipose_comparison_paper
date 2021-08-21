####Figure 2.7.1 PEER factor figure
################################################################################
f = '~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/muscle_mrna_peer.tsv'
res = read.table(f,header = T,stringsAsFactors = F,sep = '\t')
res$Tissue = 'Muscle'
res$Type = 'mRNA'
dat = data.frame(rbind(res),stringsAsFactors = F)

f = '~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/adipose_mrna_peer.tsv'
res = read.table(f,header = T,stringsAsFactors = F,sep = '\t')
res$Tissue = 'Adipose'
res$Type = 'mRNA'
dat = data.frame(rbind(dat,res),stringsAsFactors = F)

f = '~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/muscle_mirna_peer.tsv'
res = read.table(f,header = T,stringsAsFactors = F,sep = '\t')
res$Tissue = 'Muscle'
res$Type = 'miRNA'
dat = data.frame(rbind(dat,res),stringsAsFactors = F)

f = '~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/adipose_mirna_peer.tsv'
res = read.table(f,header = T,stringsAsFactors = F,sep = '\t')
res$Tissue = 'Adipose'
res$Type = 'miRNA'
dat = data.frame(rbind(dat,res),stringsAsFactors = F)

f = '~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/muscle_dname_peer.tsv'
res = read.table(f,header = T,stringsAsFactors = F,sep = '\t')
res$Tissue = 'Muscle'
res$Type = 'DNAme'
dat = data.frame(rbind(dat,res),stringsAsFactors = F)

f = '~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/adipose_dname_peer.tsv'
res = read.table(f,header = T,stringsAsFactors = F,sep = '\t')
res$Tissue = 'Adipose'
res$Type = 'DNAme'

dat = data.frame(rbind(dat,res),stringsAsFactors = F)
dat$Tissue = factor(dat$Tissue,levels=c('Muscle','Adipose'))
dat$Type = factor(dat$Type, levels = c('mRNA','miRNA','DNAme'))

pdf(file='figs/peer.pdf', height=cm(3.5), width=cm(6))
ggcust(dat,aes(x=num_peer,y=num_eGene,color=Tissue)) +
  geom_point(aes(color=Tissue))+
  geom_line(alpha=0.5,aes(color=Tissue))+
  facet_wrap(~ Tissue + Type ,nrow = 2,scales='free')+
  labs(x="PEER factors", y="Number of molecular traits") +
  scale_colour_brewer(palette = col_palette,drop=FALSE)+
  theme(axis.title=element_text(size=30))
dev.off()

#### Figure 2.7.2 Cumulative fraction of reads comparsion between mRNA and miRNA
################################################################################
f='/net/snowwhite/home/FUSION/Tissue/freeze5/muscle/gene/freeze5.muscle.analysis.counts.all_genes.csv'
#f='~/adipose_snk/data_201906/freeze5.adipose.analysis.counts.all_genes_280.tab.gz'
mrna_rc=data.frame(fread(f, header=T, stringsAsFactors=F)) 
rownames(mrna_rc)=mrna_rc[,1]
mrna_rc=mrna_rc[,-1]
load('~/fusion_mus_adi_cmp/low_exp_thresh_final/keep_traits.RData')
mrna_rc=subset(mrna_rc,rownames(mrna_rc) %in% keep_traits$mus_mrna_qtl)
#mrna_rc=subset(mrna_rc,rownames(mrna_rc) %in% keep_traits$adi_mrna_qtl)
mrna_rc_ori=mrna_rc
mrna_rc_prop=sweep(mrna_rc, 2, colSums(mrna_rc), FUN = '/')
mrna_rc_prop=as.matrix(mrna_rc_prop)
mrna_rc_df=apply(mrna_rc_prop, 2, function(x) sort(x,decreasing = T))
mrna_rc_df= apply(mrna_rc_df,2,cumsum)
mrna_rc_df= data.frame(mrna_rc_df,stringsAsFactors = F)

cumfun=function(x)min(which(x > 0.9))
mrna_rc_num= apply(mrna_rc_df,2,cumfun)
mean(mrna_rc_num) 

mrna_rc_df$num=seq(1,nrow(mrna_rc_df))
df=gather(mrna_rc_df,condi,prop,M12001:MP14793)
#df=gather(mrna_rc_df,condi,prop,A12001:AP14804)

png('figs/cumsum_mrna.png')
#png('figs/cumsum_mrna_adi.png')
ggcust(data=df, aes(x=num,y=prop))+
  geom_line(aes(colour=condi),alpha=0.5) + 
  geom_hline(yintercept=0.9, color = "blue",alpha=0.5) +   
  labs(title = 'mRNA',x="Gene rank", y="Cumulative fraction of reads")+
  scale_y_continuous( limits = c(0,1),breaks = seq(0,1,0.2) )+
  scale_x_continuous(breaks = seq(0,30000,5000),
                     labels = c("0","5K","10K","15K","20K","25K","30K"))+
  #theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none")
dev.off()

f='~/mus_exeRpt/data_201906/exceRpt_miRNA_ReadCounts_no_dup_pos_290.tab.gz'
#f='~/adi_exeRpt/data_201906/exceRpt_miRNA_ReadCounts_no_dup_pos_263.tab.gz'
mirna_rc=data.frame(fread(f, sep="\t", header=T, stringsAsFactors=F)) 
rownames(mirna_rc)=mirna_rc[,1]
mirna_rc=mirna_rc[,-1]
mirna_rc=subset(mirna_rc,rownames(mirna_rc) %in% keep_traits$mus_mirna_qtl)
#mirna_rc=subset(mirna_rc,rownames(mirna_rc) %in% keep_traits$adi_mirna_qtl)
mirna_rc_ori=mirna_rc
mirna_rc_prop=sweep(mirna_rc, 2, colSums(mirna_rc), FUN = '/')
mirna_rc_prop=as.matrix(mirna_rc_prop)
mirna_rc_df=apply(mirna_rc_prop, 2, function(x) sort(x,decreasing = T))
mirna_rc_df= apply(mirna_rc_df,2,cumsum)
mirna_rc_df= data.frame(mirna_rc_df,stringsAsFactors = F)
mirna_rc_num= apply(mirna_rc_df,2,cumfun)

mirna_rc_df$num=seq(1,nrow(mirna_rc_df))
df=gather(mirna_rc_df,condi,prop,M12006_mi:M12167_mi)
#df=gather(mirna_rc_df,condi,prop,A12003_mi:AP14804_mi)

png('figs/cumsum_mirna.png')
#png('figs/cumsum_mirna_adi.png')
ggcust(data=df, aes(x=num,y=prop))+
  geom_line(aes(colour=condi),alpha=0.5) + 
  geom_hline(yintercept=0.9, color = "blue",alpha=0.5) +   
  labs(title = 'miRNA',x="Gene rank", y="Cumulative fraction of reads")+
  scale_y_continuous( limits = c(0,1),breaks = seq(0,1,0.2) )+
  scale_x_continuous(breaks = seq(0,800,100))+
  #theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")
dev.off()


####Figure 2.7.3 Number of predicted target mRNAs and miRNA log10 mean ~ read count
################################################################################
load('~/fusion_mus_adi_cmp/TargetScan/mus_mirna.RData')
load('~/fusion_mus_adi_cmp/TargetScan/adi_mirna.RData')

pdf('figs/mus_mirna_target_rc_TarBase.pdf')
ggcust(mus_mirna,aes(x=log10mean,y=TarBase_num_targets))+
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se=T)+
  labs(title = 'Muscle miRNA',x="Log10 mean read count", 
       y="Number of target mRNAs")+
  scale_x_continuous(limits = c(0,7),breaks = seq(0,7,1))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()  

pdf('figs/adi_mirna_target_rc_TarBase.pdf')
ggcust(adi_mirna,aes(x=log10mean,y=TarBase_num_targets))+
  geom_point(alpha=0.5) +
  geom_smooth(method = "lm", se=T)+
  labs(title = 'Adipose miRNA',x="Log10 mean read count", 
       y="Number of target mRNAs")+
  scale_x_continuous(limits = c(0,7),breaks = seq(0,7,1))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()  


####Figure 2.7.12 muscle fissue and fiber plots
################################################################################
f='/net/snowwhite/home/FUSION/Tissue/freeze5/muscle/muscle_info_for_trait_assoc_with_tissueDESeq2_fibertype.tsv'
res=read.table(f,stringsAsFactors = F,header = T)
df<- res[,c('analysis_id','AdiposeSubcut','LymphocyteEBVtx','MuscleSkeletal','BloodWhole')]
colnames(df)[-1]=c('Adipose','Lymphocyte','Muscle','Blood')
df$Skin=1-rowSums(df[,2:5])
df = df[order(df$Muscle, df$Blood, df$Adipose),]
df <- melt(df)

tmp=subset(df,variable=='Muscle')
tmp=tmp[order(tmp$value,decreasing = T),]
positions <- tmp$analysis_id

df$variable <- factor(df$variable , levels = c('Lymphocyte','Blood','Skin',"Adipose","Muscle"))
pdf("figs/sub_muscle_tissue_type.pdf")
ggcust(data=df, aes(x=analysis_id, y=value, fill=variable ))+
  geom_bar(stat="identity") +
  ylab('Tissue type proportions') +
  xlab('Sample')+
  theme(axis.text.x=element_blank()) +
  scale_x_discrete(limits = positions)+
  scale_fill_brewer(palette = 'Set1',name='')
dev.off() 

##fiber type
df<- res[,c('analysis_id','type1.fiber.proxy','type2A.fiber.proxy','type2X.fiber.proxy')]
colnames(df)[-1]=c('Type1','Type2A','Type2X')
df <- melt(df)

tmp=subset(df,variable=='Type1')
tmp=tmp[order(tmp$value,decreasing = T),]
positions <- tmp$analysis_id

df$variable <- factor(df$variable , levels = c('Type1','Type2A','Type2X'))
pdf("figs/sub_muscle_fiber_type.pdf")
ggcust(data=df, aes(x=analysis_id, y=value, fill=variable ))+
  geom_bar(stat="identity") +
  ylab('Fiber type proportions') +
  xlab('Sample')+
  theme(axis.text.x=element_blank()) +
  scale_x_discrete(limits = positions)+
  scale_fill_brewer(palette = 'Set1',name='')
dev.off() 


#### Figure 2.7.13  Adipose tissue type small version
################################################################################
f='/net/snowwhite/home/FUSION/Tissue/freeze5/adipose/DESeq2_tissue_proportions_296mRNA.txt'
res=read.table(f,stringsAsFactors = F,header = T)
f = '/net/snowwhite/home/FUSION/Tissue/freeze5/adipose/freeze5_adipose_info.analysis.tin.txt'
info = read.table(f, sep="\t", header=T,stringsAsFactors=F,check.names = F)
id=info[info$use_me>0,'sample']
res=subset(res,sample %in% id)
df<- res
colnames(df)[2:4]=c('Adipocyte','Tcell','MVEC')
df <- melt(df)

tmp=subset(df,variable=='Adipocyte')
tmp=tmp[order(tmp$value,decreasing = T),]
positions <- tmp$sample

df$variable <- factor(df$variable , levels = c('Adipocyte','Tcell','MVEC','Macrophage','Blood'))
pdf("figs/sub_adipose_small.pdf")
ggcust(data=df, aes(x=sample, y=value, fill=variable ))+
  geom_bar(stat="identity") +
  ylab('Tissue/cell type proportions') +
  xlab('Sample')+
  theme(axis.text.x=element_blank()) +
  scale_x_discrete(limits = positions)+
  scale_fill_brewer(palette = 'Set1',name='')
dev.off() 

####Figure 2.7.14 Muscle physiological trait results: 3 models
####Figure 2.7.15 Adipose physiological trait results: 4 models
################################################################################
load('~/fusion_mus_adi_cmp/trait_asso/res/res_fdr01_bon.RData')
res$pct = (res$number/res$totmolec)*100
res=subset(res,res$thresh=='FDR01')

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


dat = subset(res,res$tissue == 'Muscle')
dat$model = factor(dat$model,levels = c("Base","TissueFiber","SV"))
pdf('figs/muscle_3model.pdf',width = cm(5),height = cm(4.5))
ggplot(dat,aes(x=trait_long, y=pct, fill=model,color = model)) +
  geom_bar(stat='identity',position=position_dodge(preserve = "single"))+ 
  facet_wrap(~ type,ncol = 1,scales="free_y",drop = F)+
  labs(x="",y=expression("Percent of tests (FDR"<="1%)"),
       title="Associations of physiological traits with mRNA and miRNA expression and DNA methylation
       in skeletal muscle") + 
  scale_fill_brewer(palette = "Set1",drop=FALSE,name = "Model") +
  scale_color_brewer(palette = "Set1",drop=FALSE,name = "Model") +
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

dat = subset(res,res$tissue == 'Adipose')
dat$model = factor(dat$model,levels = c("Base","Component_5","Component_17","SV"))
pdf('figs/adipose_4model.pdf',width = cm(5),height = cm(4.5))
ggplot(dat,aes(x=trait_long, y=pct, fill=model,color = model)) +
  geom_bar(stat='identity',position=position_dodge(preserve = "single"))+ 
  facet_wrap(~ type,ncol = 1,scales="free_y",drop = F)+
  labs(x="",y=expression("Percent of tests (FDR"<="1%)"),
       title="Associations of physiological traits with mRNA and miRNA expression and DNA methylation
       in subcutaneous adipose") + 
  scale_fill_brewer(palette = "Set1",drop=FALSE,name = "Model") +
  scale_color_brewer(palette = "Set1",drop=FALSE,name = "Model") +
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

####Figure 2.7.16-19 pairwise scatterplot comparing fasting insulin/BMI- mRNA associations between differnt models 
################################################################################
#f='~/fusion_mus_adi_cmp/trait_asso/adipose_4_model_20201020.tab.gz'
#f='~/fusion_mus_adi_cmp/trait_asso/adipose_4_model_20201020_bmi.tab.gz'
#f='~/fusion_mus_adi_cmp/trait_asso/muscle_3_model_insu_20201020.tab.gz'
f='~/fusion_mus_adi_cmp/trait_asso/muscle_3_model_bmi_20201020.tab.gz'
df=data.frame(fread(f, sep="\t", header=T, stringsAsFactors=F,check.names = F),check.names = F)
#tmp=df[,c('sign_mlog10p_base','sign_mlog10p_tf','sign_mlog10p_sv')]
#colnames(tmp)=c('Base','TissueFiber','SV')
tmp=df[,c('sign_mlog10p_base','sign_mlog10p_small','sign_mlog10p_paivi','sign_mlog10p_sv')]
colnames(tmp)=c('Base','Component_5','Component_17','SV')

pm=ggpairs(tmp)+theme_bw()
pm2 <- pm
for(i in 2:pm$nrow) {
  for(j in 1:(i-1)) {
    pm2[i,j] <- pm[i,j] +
      #scale_x_continuous(limits = c(-24,24)) +  #mus#change this number by looking at summary(tmp)
      #scale_y_continuous(limits = c(-24,24))+#mus
      scale_x_continuous(limits = c(-28,28))+
      scale_y_continuous(limits = c(-28,28))+
      geom_abline(intercept=0, slope = 1, color = "blue",size=0.5,alpha=0.5)+
      geom_abline(intercept=0, slope = -1, color = "blue",size=0.5,alpha=0.5)+
      geom_hline(yintercept=0, color = "blue",alpha=0.5) +
      geom_vline(xintercept=0, color = "blue",alpha=0.5) 
  }
}
#png('figs/sub_adi_4model_scatterplot_insu.png')
#png('figs/sub_adi_4model_scatterplot_bmi.png')
#png('figs/sub_mus_3model_scatterplot_insu.png')
png('figs/sub_mus_3model_scatterplot_bmi.png')
print(pm2)
dev.off()


#### Figure 2.7.20-21 mRNAs/miRNAs/DNAme sites associated with physiological traits in muscle and adipose: 
#### without or with additional adjustment of fasting serum insulin or BMI
################################################################################

##these two figures are the same figures with different data
load('~/fusion_mus_adi_cmp/trait_asso/res/res_fdr01_bon_bmi_insu_adj_noadj.RData')
res$type = factor(res$type,levels = c('mRNA','miRNA','DNAme'))
res$tissue = factor(res$tissue,levels = c('Muscle','Adipose'))
res$model = factor(res$model,levels = c("TissueFiber","Component_17",
                                        "Fasting_insulin_adjusted","BMI_adjusted"))
res$pct= (res$number/res$totmolec)*100
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

df=subset(res,thresh=='FDR01' & tissue=='Muscle')
df$model = factor(df$model,levels = c("TissueFiber","Fasting_insulin_adjusted","BMI_adjusted"))

pdf('~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/figs/muscle_insu_bmi_adj.pdf',
    width = cm(5),height = cm(4.5))
ggplot(df,aes(x=trait_long, y=pct, fill=model,color = model)) +
  geom_bar(stat='identity',position=position_dodge(preserve = "single"))+ 
  facet_wrap(~ type,ncol = 1,scales="free_y",drop = F)+
  labs(x="",y=expression("Percent of tests (FDR"<="1%)"),
       title="Associations of physiological traits with molecular trait levels in skeletal muscle\n
       without and with additional adjustment of fasting serum insulin or BMI") + 
  scale_fill_brewer(palette = "Set1",drop=FALSE,name = "Model") +
  scale_color_brewer(palette = "Set1",drop=FALSE,name = "Model") +
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

df=subset(res,thresh=='FDR01' & tissue=='Adipose')
df$model = factor(df$model,levels = c("Component_17","Fasting_insulin_adjusted","BMI_adjusted"))

pdf('~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/figs/adipose_insu_bmi_adj.pdf',
    width = cm(5),height = cm(4.5))
ggplot(df,aes(x=trait_long, y=pct, fill=model,color = model)) +
  geom_bar(stat='identity',position=position_dodge(preserve = "single"))+ 
  facet_wrap(~ type,ncol = 1,scales="free_y",drop = F)+
  labs(x="",y=expression("Percent of tests (FDR"<="1%)"),
       title="Associations of physiological traits with molecular trait levels in subcutaneous adipose\n
       without and with additional adjustment of fasting serum insulin or BMI") + 
  scale_fill_brewer(palette = "Set1",drop=FALSE,name = "Model") +
  scale_color_brewer(palette = "Set1",drop=FALSE,name = "Model") +
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

#### Figure 2.7.22-25 fasting insulin or BMI ~ mRNA associations: without and without addjustment comparison
################################################################################

#These four figures are the same figure with different data

f='~/Tissue_Li/Trait_mRNA/muscle/results/tissuefiber_S_Insu.results.tab'
#f='~/Tissue_Li/Trait_mRNA/muscle/results/tissuefiber_bmi.results.tab'
no_adj=data.frame(fread(f, sep="\t", header=T, stringsAsFactors=F,
                        check.names = F),check.names = F)
#f='~/Tissue_Li/Trait_mRNA/muscle_adj_bmi/results/tissuefiber_S_Insu.results.tab'
f='~/Tissue_Li/Trait_mRNA/muscle_adj_insu/results/tissuefiber_bmi.results.tab'
bmi_adj=data.frame(fread(f, sep="\t", header=T, stringsAsFactors=F,
                         check.names = F),check.names = F)
colnames(no_adj)[-1]=paste0('no_adj_',colnames(no_adj)[-1])
colnames(bmi_adj)[-1]=paste0('bmi_adj_',colnames(bmi_adj)[-1])

df=merge(no_adj,bmi_adj,by='response')

df$sign_mlog10p_no_adj_p.value=ifelse(df$no_adj_estimate > 0,-log10(df$no_adj_p.value),log10(df$no_adj_p.value) )
df$sign_mlog10p_bmi_adj_p.value=ifelse(df$bmi_adj_estimate > 0,-log10(df$bmi_adj_p.value),log10(df$bmi_adj_p.value) )

lim_value = ceiling(max(abs(df$sign_mlog10p_no_adj_p.value),abs(df$sign_mlog10p_bmi_adj_p.value)))

#png('~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/figs/mus_s_insu_bmiadj.png')
png('~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/figs/mus_bmi_insuadj.png')
ggcust(df,aes(x=sign_mlog10p_no_adj_p.value, y=sign_mlog10p_bmi_adj_p.value)) +
  geom_point(alpha=0.6,size=2)+ 
  geom_abline(intercept=0, slope = 1, color = "blue",size=0.5,alpha=0.5)+
  geom_abline(intercept=0, slope = -1, color = "blue",size=0.5,alpha=0.5)+
  geom_hline(yintercept=0, color = "blue",alpha=0.5) +
  geom_vline(xintercept=0, color = "blue",alpha=0.5) +
  labs(x="No additional adjustment: signed -log10(P)", 
       #y="Additional adjustment of BMI: signed -log10(P)",
       y="Additional adjustment of fasting serum insulin\nsigned -log10(P)",
       #title="Muscle\nFasting serum insulin-mRNA associations",
       title="Muscle BMI-mRNA associations") + 
  xlim(-lim_value,lim_value) +
  ylim(-lim_value,lim_value) 
dev.off()

colnames(df)[-1]=paste0('mus_',colnames(df)[-1])

#f='~/Tissue_Li/Trait_mRNA/adipose/results/paivi_S_Insu.results.tab'
f='~/Tissue_Li/Trait_mRNA/adipose/results/paivi_bmi.results.tab'
no_adj=data.frame(fread(f, sep="\t", header=T, stringsAsFactors=F,
                        check.names = F),check.names = F)
#f='~/Tissue_Li/Trait_mRNA/adipose_adj_bmi/results/paivi_S_Insu.results.tab'
f='~/Tissue_Li/Trait_mRNA/adipose_adj_insu/results/paivi_bmi.results.tab'
bmi_adj=data.frame(fread(f, sep="\t", header=T, stringsAsFactors=F,
                         check.names = F),check.names = F)

colnames(no_adj)[-1]=paste0('no_adj_',colnames(no_adj)[-1])
colnames(bmi_adj)[-1]=paste0('bmi_adj_',colnames(bmi_adj)[-1])

df=merge(no_adj,bmi_adj,by='response')

df$sign_mlog10p_no_adj_p.value=ifelse(df$no_adj_estimate > 0,-log10(df$no_adj_p.value),log10(df$no_adj_p.value) )
df$sign_mlog10p_bmi_adj_p.value=ifelse(df$bmi_adj_estimate > 0,-log10(df$bmi_adj_p.value),log10(df$bmi_adj_p.value) )

lim_value = ceiling(max(abs(df$sign_mlog10p_no_adj_p.value),abs(df$sign_mlog10p_bmi_adj_p.value)))

#png('~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/figs/adi_s_insu_bmiadj.png')
png('~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/figs/adi_bmi_insuadj.png')
ggcust(df,aes(x=sign_mlog10p_no_adj_p.value, y=sign_mlog10p_bmi_adj_p.value)) +
  geom_point(alpha=0.6,size=2)+ 
  geom_abline(intercept=0, slope = 1, color = "blue",size=0.5,alpha=0.5)+
  geom_abline(intercept=0, slope = -1, color = "blue",size=0.5,alpha=0.5)+
  geom_hline(yintercept=0, color = "blue",alpha=0.5) +
  geom_vline(xintercept=0, color = "blue",alpha=0.5) +
  labs(x="No additional adjustment: signed -log10(P)", 
       #y="Additional adjustment of BMI: signed -log10(P)",
       y="Additional adjustment of fasting serum insulin\nsigned -log10(P)",
       #title='Adipose\nFasting serum insulin-mRNA associations',
       title='Adipose BMI-mRNA associations') + 
  xlim(-20,20) +
  ylim(-20,20) 
dev.off()

#### Figure 2.7.26 hsa-miR-122-5p ~ physiological trait associations
################################################################################

f='~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/has_miR_122_5p.tab.gz'
df=data.frame(fread(f, sep="\t", header=T, stringsAsFactors=F))
df$term[2]='T2D_NGT'

int=c("bmi","DIo","fS_C_pept","fS_C_pept_30","fS_Kol_HDL","fS_Trigly",                    
      "GL120","GL30","GL60","glu_2h_biopsy","Glu_AUC_0to30","HIP" ,"HOMA",
      "Ins_AUC_0to30","InsGenIn","InsSec30","matsuda_3pt","matsuda_4pt",                  
      "p_insu","p_insu_120","p_insu_30","p_insu_60","RFM","S_ALAT",                       
      "S_GT","S_hs_CRP","S_Insu","S_Insu_30","S_LipoA1","S_Uraat",                    
      "T2D_NGT","WAIST","WEIGHT","whr","GL0","glu_fast_biopsy","S_LipoB",
      "ApoB_A1_ratio","whradjbmi","B_GHb_A1C","fS_Kol_LDL_c","B_HbA1c",                     
      "sbp","fS_Krea","HEIGHT","fS_Kol","dbp","CpepGenIn")   

df=subset(df,term %in% int)

df$sign_mlog10p_mus = ifelse(df$mus_estimate >0, -log10(df$mus_p.value),log10(df$mus_p.value))
df$sign_mlog10p_adi = ifelse(df$adi_estimate >0, -log10(df$adi_p.value),log10(df$adi_p.value))

#muscle adipose associations plots
pdf('~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/1225p_mus_adi_cmp_repel.pdf')
ggplot(df,aes(x=sign_mlog10p_mus, y=sign_mlog10p_adi)) +
  geom_point(size=1,alpha = 0.8)+ 
  geom_label_repel(aes(label=trait_long),color = 'black',size = 2,max.overlaps=40)+
  geom_abline(intercept=0, slope = 1, color = "red")+
  geom_hline(yintercept=0, color = "blue")+
  geom_vline(xintercept=0, color = "blue")+
  labs(x="Muscle signed -log10(P)", y="Adipose signed -log10(P)") + 
  theme_bw() + 
  theme(strip.text=element_text(size=15,face = 'bold')) + 
  scale_y_continuous(breaks = seq(-2,12,2),limits = c(-2,12)) +
  scale_x_continuous(breaks = seq(-2,12,2),limits = c(-2,12))
dev.off()

####Figure 2.7.27 Association of hsa-miR-122-5p ~ ALT colored by HBB
################################################################################
f='~/fusion_mus_adi_cmp/trait_asso/HBB_asso/mus_hsa_miR_122_5p.tab.gz'
mus_hsa_miR_122_5p=data.frame(fread(f, sep="\t", header=T, stringsAsFactors=F))
f='~/fusion_mus_adi_cmp/trait_asso/HBB_asso/adi_hsa_miR_122_5p.tab.gz'
adi_hsa_miR_122_5p=data.frame(fread(f, sep="\t", header=T, stringsAsFactors=F))

pdf('~/fusion_mus_adi_cmp/paper_figs/paper_main_sub_figures_202011/figs/1225p_mus_adi_hbb.pdf',width=cm(7.2),height = cm(3.2))
df = mus_hsa_miR_122_5p
p1=ggcust(df, aes(x=S_ALAT,y=hsa.miR.122.5p, colour=HBB)) + 
  geom_point(size=2) +
  scale_colour_gradient2(low= muted("blue"),mid='white', high = muted('red')) +
  labs(title = 'Muscle',
       x="Alanine aminotransferase (ALT)", 
       y="hsa-miR-122-5p") 
df = adi_hsa_miR_122_5p
p2=ggcust(df, aes(x=S_ALAT,y=hsa.miR.122.5p, colour=HBB)) + 
  geom_point(size=2) +
  scale_colour_gradient2(low= muted("blue"),mid='white', high = muted('red')) +
  labs(title = 'Adipose',
       x="Alanine aminotransferase (ALT)", 
       y="hsa-miR-122-5p") 
grid.arrange(p1,p2,nrow=1)
dev.off()



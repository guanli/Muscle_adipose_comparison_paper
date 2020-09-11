# Set_low_exp_thresh
3 steps:
1. plot the relationship for each bin between the proportion of molecular traits having QTL and mean read count of that bin
adipose miRNA
new 1) can do also for DNAme rank the DNAme sites according to varianceRscript
2) expression level file can have more genes than in QTL results as extra genes are filtered out
~/common_scripts/low_exp_BHFDR_each_bin.R -i ~/adipose_dname/eQTL_201906/cis/pct10_qtltools_nominal.tsv.gz -e /net/snowwhite/home/FUSION/Tissue/Methylation/EPIC/processed_data/methyl_v003/epic-beta-A.csv.gz -n 200 -b 2 -g probe -s variance -o ~/fusion_mus_adi_cmp/low_exp_thresh/gam_adipose_dname_beta_variance &

old for mRNA and miRNA only Rscript ~/common_scripts/low_exp_BHFDR_each_bin.R 
-i ~/adi_exeRpt/eQTL_low_exp_201906/cis/qtltools_nominal.tsv.gz
-e ~/adi_exeRpt/data_201906/exceRpt_miRNA_ReadCounts_no_dup_pos_263.tab.gz 
-n 20 -b 2 -g miRNA -o ~/fusion_mus_adi_cmp/low_exp_thresh_final/adipose_mirna_rc &

2. Manually Look at starting from which bin, there is no longer 0 proportion of molecular traits having QTL
f = '~/fusion_mus_adi_cmp/low_exp_thresh_final/adipose_mirna_rc_gam_pred_prop.tab.gz'
res = read.table(f, header = T, sep = '\t', stringsAsFactors = F)
table(res$bin)#1-19 bin 99 genes last 94 genes
dat = aggregate(res[,"sig"], by=list(res$bin), FUN=mean) 
subset(dat,dat[,2]>0) #from bin 9

3. Use different cut off point to evaluate its effects on the detection rate
Rscript ~/common_scripts/mol_traits_diff_low_exp_cutoff.R 
-i ~/adi_exeRpt/eQTL_low_exp_201906/cis/qtltools_nominal.tsv.gz
-e ~/adi_exeRpt/data_201906/exceRpt_miRNA_ReadCounts_no_dup_pos_263.tab.gz 
-p ~/fusion_mus_adi_cmp/low_exp_thresh_final/adipose_mirna_rc_gam_pred_prop.tab.gz
-n 9,10,11,12,13 -b 2 -g miRNA -t 0.05,0.01 
-o ~/fusion_mus_adi_cmp/low_exp_thresh_final/adipose_mirna_rc


The last step can also be used to evaluate its effects on the detection rate for tratis associaionts
Rscript ~/common_scripts/mol_traits_diff_low_exp_cutoff.R 
-i ~/muscle_trans_snk/eQTL_low_exp_201906/tissuefiber_qts.S_Insu.results.tab.gz
-e /net/snowwhite/home/FUSION/Tissue/freeze5/muscle/gene/freeze5.muscle.analysis.counts.all_genes.csv.gz
-p ~/fusion_mus_adi_cmp/low_exp_thresh_final/muscle_mrna_rc_gam_pred_prop.tab.gz
-n 1,35,40,45,50,55 -b 2 -g gene -t 0.05,0.01 
-o ~/fusion_mus_adi_cmp/low_exp_thresh_final/muscle_fast_insulin_rc

This pipleline is used to determined the threshold for filtering out lowly expressed genes or for filtering out DNA methylation (DNAme) sites with low levels of variations in beta-values. It first portrays the relationship between the eQTL detection rate and the mean of read count or variance of beta-values by <br />
  1. Order genes/DNAme sites by their mean read count or variance of beta-values.<br />
  2. Group genes into bins of equal number of genes/DNAme sites.<br />
  3. Calculate the proporting of genes with QTL using the Benjamini-hochberg FDR (BH-FDR) procedure. <br />
  4. Use the generate additive model to fit a smooth spline between the eQTL detection rates and mean read count or variance of beta-values. <br />

Example command of this step for genes (mRNAs or miRNAs)<br />
  `~/common_scripts/low_exp_BHFDR_each_bin.R -i ~/adi_exeRpt/eQTL_low_exp_201906/cis/qtltools_nominal.tsv.gz -e ~/adi_exeRpt/data_201906/exceRpt_miRNA_ReadCounts_no_dup_pos_263.tab.gz -n 20 -b 2 -g miRNA -s mean -o ~/fusion_mus_adi_cmp/low_exp_thresh_final/adipose_mirna_rc`<br />

Example command of this step for genes (mRNAs or DNAme sites)<br />
  `~/common_scripts/low_exp_BHFDR_each_bin.R -i ~/adipose_dname/eQTL_201906/cis/pct10_qtltools_nominal.tsv.gz -e /net/snowwhite/home/FUSION/Tissue/Methylation/EPIC/processed_data/methyl_v003/epic-beta-A.csv.gz -n 200 -b 2 -g probe -s variance -o ~/fusion_mus_adi_cmp/low_exp_thresh/gam_adipose_dname_beta_variance`<br />

Command line arguments <br /> 
  `-i` The cis-QTL associations of all SNP-gene/DNAme sites pairs (QTLoutput format) <br />
  `-e` The gene read count or DNAme sites beta-values <br />
  `-n` Number of genes or DNAme sites within each bin. Currently the bin sizes I am using are 20 for miRNA, 200 for DNAme sites and 100 for mRNA. <br />
  `-b` Tcolumn in exp_level file where the number begins <br />
  `-g` The column name of genes/DNAme sites names <br />
  `-s` Use mean or variance <br />
  `-o` Output prefix <br />


Evaluate the effects of using different cut-off points on the QTL detection rate by including genes or DNAme sites from the first bin where the proportion of genes or DNAme sites with QTLs were no longer zero, until about half of the genes or DNAme sites were included.

Example command of this step <br />
  `Rscript ~/common_scripts/mol_traits_diff_low_exp_cutoff.R -i ~/adi_exeRpt/eQTL_low_exp_201906/cis/qtltools_nominal.tsv.gz -e ~/adi_exeRpt/data_201906/exceRpt_miRNA_ReadCounts_no_dup_pos_263.tab.gz -p ~/fusion_mus_adi_cmp/low_exp_thresh_final/adipose_mirna_rc_gam_pred_prop.tab.gz -n 9,10,11,12,13 -b 2 -g miRNA -t 0.05,0.01 -o ~/fusion_mus_adi_cmp/low_exp_thresh_final/adipose_mirna_rc`<br /> 

Command line arguments <br /> 
  `-p` Output from the first step <br />
  `-n` The starting bin from which to include genes or DNAme sites <br />
  `-t` BH-FDR threshold options <br />


The second step can also be used to evaluate its effects on the detection rate for tratis associaionts with the example command as shown below.<br />
`Rscript ~/common_scripts/mol_traits_diff_low_exp_cutoff.R -i ~/muscle_trans_snk/eQTL_low_exp_201906/tissuefiber_qts.S_Insu.results.tab.gz -e /net/snowwhite/home/FUSION/Tissue/freeze5/muscle/gene/freeze5.muscle.analysis.counts.all_genes.csv.gz -p ~/fusion_mus_adi_cmp/low_exp_thresh_final/muscle_mrna_rc_gam_pred_prop.tab.gz -n 1,35,40,45,50,55 -b 2 -g gene -t 0.05,0.01  -o ~/fusion_mus_adi_cmp/low_exp_thresh_final/muscle_fast_insulin_rc`

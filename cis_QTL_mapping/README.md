
This pipeline can perform the following two analyses.
 - cis-eQTL nominal pass  <br />
  `snakemake qtltools_proximal__nominal_concat --snakefile Snakefile --configfile config_analysis.json`  <br />
 - cis-eQTL permutation pass: An approximate permutation analysis to get a pvalue adjusted for the number of tested SNPs for each gene  <br />
  `snakemake qtltools_proximal__permute_fdr --snakefile Snakefile --configfile config_analysis.json`

Before running the pipeline
 - Create a folder called data and put gene expression, covairate and genotype files into this folder.
 - If running the permutation pass Creat a folder called scripts, place there the runFDR_cis.R downloaded from (https://qtltools.github.io/qtltools/) and name it qtltools-runFDR_cis.R. This script is what QTLtools uses to correct for the number of tested genes through Storey & Tibshirani False Discovery Rate procedure (ST_FDR).
 
Parameters in the config_analysis.json file
 - `NPERM`: Number of permutations to run, a large number is recommended for final QTL mapping, e.g. 10000. <br />
 - `PARAM_QTLTOOLS_NOMINAL_CUTOFF`: Store the associations that pass this threshold, if want to store all of the results, use 1.

Result files that will appear after this pipeline finished. 
 - cis-eQTL nominal pass  <br />
   - `qtltools_nominal.tsv.gz` cis-eQTL associations <br />
 - cis-eQTL permutation pass  <br />
   - `qtltools_permute-significant.tsv.gz`  <br /> 
      - After applying the ST_FDR procedure, the genes are considered to be significant 
   - `qtltools_permute-thresholds.txt.gz`  <br /> 
      - Gene-specific significant threshold from the permutation tests
   - `qtltools_permute.tsv.gz` <br /> 
      - The cis-eQTL mapping results for the lead SNP (the SNP with the most significant pvalues) of each gene 
 
 


This pipeline finds the number of PEER factors that optimizes QTL detection rate. 

This piplefile 1) performs LD-pruning to the input genotype file, and then 2) use the pruned list of genetic variants, the input gene expression file and the input covariate file to conduct QTL mapping using each of the PEER factors specified and 3) generate QTL mapping results into ../peer_runs/factor_k/ folder. Therefore, by comparing the number of genes having eQTL or the number of significant SNP-gene associations, the number of PEER factors that optimizes QTL detection rate can be determined. 

Run the Snakemake pipeline with the following command.
  `snakemake --snakefile peer_opt_cis_perm_based --configfile peer_opt_trans.json -j 12 --latency-wait 60 > snk_jobs.txt 2>&1 &`
  
Before running this pipeline
  - Create a folder called data/ and put original genotype file (not pruned) genotypes.vcf.gz, gene expression level file moltraits_trans.bed.gz there
  - Create a folder called scripts and place Li_peer-optimize_perm_plot.R, qtltools-runFDR_cis.R and realpath files there. 

Options in the configuration file
  - `INPUT_PREP_PHENO_QTLTOOLS_START_COL`: the column number that the first gene expression value column is, e.g. 2
  - `INPUT_PREP_PHENO_QTLTOOLS_FEATURE_COL`: The colomn name of the gene names, e.g. "gene"
  - `PEER_FACTORS_OPT`: the list of the number of PEER factors to try, e.g. ["0","1","2","3","4","5","6"]
  - `NJOBS_NOMINAL` : the number of cores to run the Snakemake pipeline
  - `PLINK_LD_PRUNE` : option for plink-1.9 --vcf --indep-pairwise, used in the LD pruning step.
  
Output results in each PEER factor folder
  - qtltools_permute-significant-summary.tsv
  - qtltools_permute-significant.tsv.gz
  - qtltools_permute-thresholds.txt.gz
  - qtltools_permute.tsv.gz
  
Output in the current folder
  - peer_factors-qtltools_perm-summary.tsv: the number of genes with QTL from result using different number of PEER factors
  - eGene_cis_peer_factors-qtltools_perm-summary.pdf: a line plot displaying the above results. 


This pipeline finds the number of PEER factors that optimizes QTL detection rate. 

This piplefile first performs LD-pruning to the input genotype file, and then use the pruned list of genetic variants, the input gene expression file and the input covariate file to conduct QTL mapping using each of the PEER factors specified and generate QTL mapping results into ../peer_runs/factor_k/ folder.

Before running this pipeline
- Create a folder called data/ and put original genotype file (not pruned) genotypes.vcf.gz, gene expression level file moltraits_trans.bed.gz there
- Create a folder called scripts and place Li_peer-optimize_perm_plot.R, qtltools-runFDR_cis.R and realpath files there. 

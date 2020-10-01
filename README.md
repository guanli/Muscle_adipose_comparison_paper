# Muscle_adipose_comparison_paper
This repository contains code used large-scale analyses for Guan et al. FUSION muscle-adipose comparison paper, including code to 
Generate PEER factors --> PEER_generation/
Find the number of PEER factors that optimize QTL detection rate --> PEER_optimization/
Adipose tissue deconvolution --> adipose_decon_4_ct/
Set a reasonable threshold to exclude lowly expressed genes --> set_low_exp_thresh
signle-variant cis QTL mapping --> cis_QTL_mapping/
Multiple-variant cis QT mapping --> dap/

The analysis pipelines were built using Snakemake and thus there is a Snakefile and a corresponding configuration file in each folder alongwith a configuration file. 

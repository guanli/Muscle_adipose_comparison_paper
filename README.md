# Muscle_adipose_comparison_paper
This repository contains analysis pipelines used for large-scale analyses for the FUSION muscle-adipose comparison paper, including code to 
- Generate PEER factors --> PEER_generation/
- Find the number of PEER factors that optimizes QTL detection rate --> PEER_optimization/
- Adipose tissue deconvolution --> adipose_decon_4_ct/
- Set a reasonable threshold to exclude lowly expressed genes --> set_low_exp_thresh
- signle-variant cis QTL mapping --> cis_QTL_mapping/
- Multiple-variant cis QT mapping --> dap/

The analysis pipelines were built using Snakemake. Each folder has a Snakefile and a corresponding configuration file. 

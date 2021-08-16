# Muscle_adipose_comparison_paper
This repository contains analysis pipelines used for large-scale analyses for the FUSION muscle-adipose comparison paper, including code to 
- Generate PEER factors --> PEER_generation/
- Find the number of PEER factors that optimized QTL detection rate --> PEER_optimization/
- Adipose tissue deconvolution using the 5-component approach --> adipose_decon_4_ct/
- Find the threshold to exclude lowly expressed genes while maximizing the power to detect molecular trait-genotype associations
- Marginal cis-QTL mapping --> cis_QTL_mapping/
- Multiple indepdent cis-QTL analysis --> dap/

The analysis pipelines were built using Snakemake. Each folder has a Snakefile and a corresponding configuration file. 

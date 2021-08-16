This pipeline generates PEER output for the number of PEER factors specified in peer_gnr_config.json.

Before running the pipeline
- Create a folder called data and put the gene expression level file into this folder
- Create a folder called peer_runs and subfolders e.g. factor_1, factor_10. These folders are where the PEER output is going to be stored. 
- Create a folder called scripts and place the peer-run.R into this folder

The snakefile first inverse normalize the input file (if inverse_norm is TRUE) and then run PEER software. 

Command to run the pipeline
  `snakemake --snakefile peer_gnr final` 

Parameters in the Peer_gnr_config.json file 
1. `PEER_FACTORS`  how many PEER factors to run
2. `PARAM_PEER` 
  - iterations 1000 --account_mean FALSE: recommended settings for PEER factors , usually don't change
  - inverse_norm TRUE: if the gene expression file you have is already inverse-normalized, you can set to FALSE
  - feature_ids_colnam: the column name of the gene name 
3. `START_COL` the column: where the first numeric value is

Five result file will appear in each result folder, including 
  - moltraits-peer_factors.tsv.gz
  - moltraits-peer_precision.tsv.gz
  - moltraits-peer_residuals-invnorm.tsv.gz
  - moltraits-peer_residuals.tsv.gz
  - moltraits-peer_weights.tsv.gz

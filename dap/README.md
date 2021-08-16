
This pipeline uses DAP (Wen et al. Am J Hum Genet. 2016, PMID 27236919) to perform cis-eQTL fine mapping.

It 1) converts the cis-eQTL input files (gene expression file, genotype file and covariate file) from the formats required by QTLtools to the formats required by DAP; 2) Run DAP to perform cis-eQTL fine mapping ; 3) Run get_cs.pl (a script from Prof.Xiaoquan Wen) to get the 95% credible sets for each gene

Command to run the pipeline snakemake `--snakefile dap_enloc final`

DAP output will be in dap_res/, where there is the \*.fm.out for each each. <br />
95% credible sets will be in cs/, where there is the \*.cs for each each that has 95% credible set.<br />

Example directory is /net/snowwhite/home/guanli/muscle_trans_snk/eQTL_201906/cis_dap.

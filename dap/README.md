
This pipeline uses DAP (Wen et al. Am J Hum Genet. 2016, PMID 27236919) to perform cis-eQTL fine mapping and then uses output from DAP to perform colocalization analysis using Fastenloc (Wen et al. Plos Genetics 2017, PMID: 28278150; Pividori et al. BioRxiv. 2020).

It 1) converts the cis-eQTL input files(gene expression file, genotype file and covariate file) from the formats required by QTLtools to the formats required by DAP; 2) Run DAP to perform cis-eQTL fine mapping and get the 95% credible sets for each gene; 3) Run fastenloc to perform colocalization analysis. 

Command to run the pipeline snakemake --snakefile dap_enloc final

DAP output will be in dap_res/, where there is the \*.fm.out for each each. <br />
95% credible sets will be in cs/, where there is the \*.cs for each each that has 95% credible set.<br />
colocalization output will be in the current folder, include \*.enloc.sig.out and \*.enloc.snp.out

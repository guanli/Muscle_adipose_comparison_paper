This pipeline uses UNMIX function in the Deseq2 package to perform tissue deconvolution of Bulk RNA-seq data. 

`Rscript deseq2_decon_final.R -m freeze5.adipose.tpm.tsv.gz -r ref.txt -c class.txt -v shift_value.txt -o kerrin_4_type_blood_tabassum_adi_296_mrna/` <br />

Input files 
- Gene expression values of each sample in your RNA-seq data 
- Gene expression values of each sample in the reference data
- A file specifying the tissue/cell type assignment of the reference file 
- A file with a range of shift_values (a key parameter to optimize recommanded by Dr.Love, developer of the UNMIX funtion).

File format (all tab deliminated) <br />
- Gene expression file to deconvolute*<br />
```
gene  A001  A002  A003 
ENSG1 30.1  140.2 86.5
ENSG2 7.8 1.3 7.8 
ENSG3 91.4  66.7  180.9
```
<br />

*Reference gene expression file* 
```
gene  Ref001  Ref002  Ref003 
ENSG1 30.1  140.2 86.5
ENSG2 7.8 1.3 7.8 
ENSG3 91.4  66.7  180.9
```
<br />

*Reference tissue/cell type label*: "1" means is this tissue/cell type, "2" means not. If this example, the first 3 reference samples are adipocyte, the middle 3 samples are myocype, and the last 3 samples are endothelial cells. 
```
adipocyte 1 1 1 2 2 2 2 2 2
myocyte 2 2 2 1 1 1 2 2 2
endothelial 2 2 2 2 2 2 1 1 1
```
<br />

*Shiftvalues*: This following example file is what I acctually used in my analysis
```
1e-04
5e-04
0.001
0.005
0.01
0.05
0.1
0.2
0.3
0.4
0.5
0.6
0.7
0.8
0.9
1
5
```

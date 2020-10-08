
**Example command:** <br />
 - Assume the conditional SNPs have same  effects across studies  <br />
    -- Assume the SNP of interest have same effects acorss studies; pvalue threshold for the step-wise selection procedure is based on ACAT pvalues (pvalues that corrected for the number of tested SNPs)  <br />
 `./yax meta --sumstats {prefix1,prefix2} --stepwise --tests hom --pvalue 0.05 --backward --prefix {output-prefix}` <br />
    -- Assume the SNP of interest have different effects acorss studies (calculate the pvalues under both het and alt assumptions ) <br /> `./yax meta --sumstats {prefix1,prefix2} --stepwise --tests het,alt --pvalue 0.05 --backward --prefix {output-prefix}` <br />
 
 - Assume the conditional SNPs have different effects across studies  <br />
    -- Assume the SNP of interest have same effects acorss studies; pvalue threshold for the step-wise selection procedure is based on marginal pvalues (raw pvalues , not corrected for the number of tested SNPs)  <br />
 `./yax meta --sumstats {prefix1,prefix2} --stepwise --tests hom --het --marginal --pvalue 2.5E-6 --backward --prefix {output-prefix}` <br />
    -- Assume the SNP of interest have different effects acorss studies; in the step-wise selection procedure, do not perform backward selection (do not drop SNPs falling the pvalue threshold in the joint model)<br />
 `./yax meta --sumstats {prefix1,prefix2} --stepwise --tests het,alt --het --pvalue 0.05 --prefix {output-prefix}` <br />
    
**Software concordance.** Regression slopes and standard errors from YAX multiple-variant meta-analysis (where all studies have unrelated samples) are equivalent to the R regression model `lm(trait ~ genotypes + study*covariates, weight = 1/study_mse )`, where `study_mse` is the mean squared error from the null model (no genotypes) fit within each study. 

## Command line arguments
```diff
- Section incomplete!
```
A partial list of options is given below.  Please run `./yax meta --help` to see a complete list of command line flags and options. 
 - **Analysis options**
 	  - `--tests=[hom,het,alt]` : Assumptions under which the pvalues for the SNPs of interestes are estimated. Comma-seperated options. Will estimate under all of the assumptions specified, i.e. [hom,het] with provide pvalues assuming incomplete! here.
	  - `--het` : if specified, assume the conditional SNPs have heterogeneous effects across studies; otherwise assume homogeneous effects. 
	  - `--rsq` : maximum multiple R2 threshold, consider only SNPs with multiple R2 less than this threshold to avoid collinearity.
	  - `--marginal` : if specified, use the raw (unadjusted for the number of tested SNPs) in the stepwise selection procedure; otherwise, use the ACAT pavlues (adjusted for the number of tested SNPs).
	  - `--pvalue`: pvalue threshold for the stepwise selection procedure. 
	  - `--backward`: if specified, perform forward and backward selection in the stepwise selection procedure; otherwise, only perform forward selection, i.e. do not drop SNPs failling the pvalue threshold in the joint model of all selected SNPs.

#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

######################## Required Packages #####################################
#example useage:
#Rscript /net/snowwhite/home/guanli/muscle_trans_snk/peer_all_gene/scripts/peer-run.R 
#--file data/moltraits.bed.gz --sample_column_start 8 
#--hidden_factors 1 --output_prefix feer_01_ori 
#--iterations 1000 --account_mean FALSE --inverse_norm TRUE 
#-ftr_col gene -smp_ex NULL -gene_ex NULL
################################################################################

######################## Required Packages #####################################
suppressMessages(library(peer))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
################################################################################

######################## Functions #############################################
`%notin%` = function(x,y) !(x %in% y)
invnorm_random = function(x) {
  qnorm((rank(x, na.last="keep", ties.method="random") - 0.5) / sum(!is.na(x)))
}

peer_model = function(mtrx, k, iterations, base, other_dat=NULL, 
                      covariates_mtx=NULL, account_mean=T, inverse_norm_before=T,feature_ids) {
  # runs peer
  # mtrx should be a matrix of samples (rows) by feature_ids (cols)
  # other_dat should be a dataframe of samples (cols) by feature_ids (rows)
  # covariates_mtx should be a matrix of samples (rows) by feature_ids (cols)
  
  if (inverse_norm_before) {
    cat("Inverse normalizing each feature id (col) across samples (row).\n")
    # apply invnorm_random to each gene (column) across samples (rows)
    mtrx = as.matrix(apply(mtrx, 2, invnorm_random))
    cat("\nstandard deviation cols:\n")
    print(head(apply(mtrx, 2, sd)))
    cat("\nstandard deviation rows:\n")
    print(head(apply(mtrx, 1, sd)))
    cat("\n")
  }
  print(paste("PEER input matrix shape:", paste(dim(mtrx), collapse=",") ))
  
  # init peer
  model = PEER()
  
  # Set observed data
  #PEER_setPhenoMean(model, mtrx)
  if (!is.null(PEER_setPhenoMean(model,mtrx))) { stop("Error setting pheno mean"); }
  # set number of hidden factors
  #PEER_setNk(model, k)
  if (!is.null(PEER_setNk(model,k))) { stop("Error setting number of factors"); }
  if (k > 0) {
    factor_names = paste("factor", seq(1,k), sep="_")
    factor_names_full = factor_names
  } else {
    factor_names = c()
    factor_names_full = c()
  }
  print(paste("PEER_factors=", (k), sep=""))
  
  # set number of max iterations
  #PEER_setNmax_iterations(model, iterations)
  if (!is.null(PEER_setNmax_iterations(model,iterations))) { stop("Error setting max iterations"); }
  # For most use cases, including the mean effect is likely to be a good choice.
  # If y is inverse normalized including a mean effect does not make sense.
  if (account_mean == T) {
    factor_names_full = c("mean", factor_names_full)
  }
  #PEER_setAdd_mean(model, account_mean)
  if (!is.null(PEER_setAdd_mean(model,account_mean))) { stop("Error adding mean factor"); }
  
  # Measured experimental variables that may contribute to variability in the 
  # data can be included in the inference as covariates
  # (e.g batch, RNA quality, ...)
  if (!is.null(covariates_mtx)) {
    cat("Setting covariates:", colnames(covariates_mtx), "\n")
    factor_names_full = c(colnames(covariates_mtx), factor_names_full)
    #print(factor_names_full)
    # if there is an error it is probly because covariates_mtx is not a matrix
    #PEER_setCovariates(model, covariates_mtx)
    if (!is.null(PEER_setCovariates(model,covariates_mtx))) { stop("Error setting covariates"); }
  }
  
  # Perform inference > results in a model object with posterior distributions 
  # of the variables
  # aka now run peer
  #PEER_update(model)
  if (!is.null(PEER_update(model))) { stop("Model failed to converge"); }
  
  # get the posterior mean of the inferred factors (NxK matrix),
  # their weights (GxK matrix),
  # precision (inverse variance) of the weights (Kx1 matrix),
  # the residual dataset (NxG matrix):
  factors = PEER_getX(model)
  weights = PEER_getW(model)
  precision = PEER_getAlpha(model)
  residuals = PEER_getResiduals(model)
  
  cat('str(factors) is', str(factors), '\n')
  cat('str(weights) is', str(weights), '\n')
  cat('str(precision) is', str(precision), '\n')
  cat('str(residuals) is', str(residuals), '\n')
  
  # add labels
  #print(factor_names_full) # print
  
  if (is.null(feature_ids)) {feature_ids=paste("row", seq(1,nrow(dat)), sep="_")}
  #samples = unlist(lapply(rownames(mtrx), FUN=function(x) { 
  #    strsplit(x, "X")[[1]][2] } ))
  samples = rownames(mtrx)
  
  factors = data.frame(t(factors))
  colnames(factors) = samples
  rownames(factors) = factor_names_full # this may give a bug
  factors$factor = factor_names_full
  # move factor column first
  factors = factors[,c(ncol(factors),1:(ncol(factors)-1))] 
  cat("processed factors\n")
  
  precision_old = precision
  precision = data.frame(precision)
  rownames(precision) = factor_names_full
  colnames(precision) = "precision"
  precision$factor = factor_names_full
  # move factor column first
  precision = precision[,c(ncol(precision),1:(ncol(precision)-1))] 
  cat("processed precision\n")
  
  weights = data.frame(weights)
  rownames(weights) = feature_ids
  colnames(weights) = factor_names_full
  weights$feature_id = feature_ids
  weights = weights[,c(ncol(weights),1:(ncol(weights)-1))]   
  cat("processed weights\n")
  
  residuals = data.frame(t(residuals), check.names=FALSE)
  colnames(residuals) = samples
  rownames(residuals) = feature_ids
  # inverse normalize residuals
  residuals_inv = data.frame(t(apply(as.matrix(residuals), 1, invnorm_random)),
                             check.names=FALSE)
  residuals$feature_id = feature_ids
  # move probe column first
  residuals = residuals[,c(ncol(residuals),1:(ncol(residuals)-1))]
  residuals_inv$feature_id = feature_ids
  residuals_inv = residuals_inv[,c(ncol(residuals_inv),1:(ncol(residuals_inv)-1))]
  colnames(residuals_inv) = colnames(residuals)
  
  cat('Before saving, str(factors) is', str(factors), '\n')
  cat('Before saving, str(weights) is', str(weights), '\n')
  cat('Before saving, str(precision) is', str(precision), '\n')
  cat('Before saving, str(residuals) is', str(residuals), '\n')

  # save our matricies 
  sep="\t"
  gzfh = gzfile(paste(base, "-peer_factors.tsv.gz", sep = ""), "w")
  write.table(factors, gzfh, row.names=F, col.names=T, quote=F, sep=sep)
  close(gzfh)
  
  gzfh = gzfile(paste(base, "-peer_weights.tsv.gz", sep = ""), "w")
  write.table(weights, gzfh, row.names=F, col.names=T, quote=F, sep=sep)
  close(gzfh)
  
  gzfh = gzfile(paste(base, "-peer_precision.tsv.gz", sep = ""), "w")
  write.table(precision, gzfh, row.names=F, col.names=T, quote=F, sep=sep)
  close(gzfh)
  
  gzfh = gzfile(paste(base, "-peer_residuals.tsv.gz", sep = ""), "w")
  write.table(residuals, gzfh, row.names=F, col.names=T, quote=F, sep=sep)
  close(gzfh)
  
  gzfh = gzfile(paste(base, "-peer_residuals-invnorm.tsv.gz", sep = ""), "w")
  write.table(residuals_inv, gzfh, row.names=F, col.names=T, quote=F, sep=sep)
  close(gzfh)
}

run_peer = function() {
  suppressPackageStartupMessages(library(optparse))
  
  optionList = list(
    make_option(c("-f", "--file"), 
                type="character", 
                help=paste0(
                  "Input molecular trait data in tsv file. ",
                  "Molecular trait features should be rows ",
                  "(i.e., each row a different trait).",
                  "It is assumed that the sample ids will be sequentially in order ",
                  "after intro columns ",
                  "(e.g., chr, start, end, feature_id, sample_1, sample_2, ...)."
                )
    ),
    
    make_option(c("-l", "--sample_column_start"), 
                type="numeric", default=2, 
                help=paste0(
                  "Column in which the sample ids begin (1 based). For instance ",
                  "(chr, start, end, feature_id, sample_1, sample_2, ...) ",
                  "would be 5, and ",
                  "(feature_id, sample_1, sample_2, ...) would be 2. [default %default]"
                )
    ),
    
    make_option(c("-k", "--hidden_factors"), 
                type="numeric", default=15, 
                help="The number of hidden factors. [default %default]"
    ),
    
    make_option(c("-i", "--iterations"), 
                type="numeric", default=1000, 
                help="Maximum number of iterations. [default %default]"
    ),
    
    make_option(c("-o", "--output_prefix"), 
                type="character", default="moltraits", 
                help="Prefix for output files. [default %default]"
    ),
    
    make_option(c("-c", "--covariate_file"), 
                type="character", default=NULL, 
                help=paste0(
                  "Input covariates file. ",
                  "Covariates should be rows ",
                  "(i.e., each row a different covariate). ",
                  "It is assumed that the sample ids start at position 2 ",
                  "(e.g., covariate_id, sample_1, sample_2, ...). ",
                  "[default %default]"
                )
    ),
    
    make_option(c("-m", "--account_mean"), 
                type="logical", default=TRUE, 
                help=paste0(
                  "Add additional covariate to account for mean y signal. ",
                  "[default %default]"
                )
    ),
    
    make_option(c("-n", "--inverse_norm"), 
                type="logical", default=FALSE, 
                help=paste0(
                  "Inverse normalize y features before running PEER ",
                  "(so each y feature has a normal distribution across samples). ",
                  "Note if true, you do not need to include a mean covariate. ",
                  "[default %default]"
                )
    ),
    
    make_option(c("-s", "--seed"), 
                type="numeric", default=1938, 
                help="Random number seed. [default %default]"
    ),
    
    make_option(c("-a", "--feature_ids_colnam"), 
                type="character", default="gene", 
                help="Colname of feature ID. [default %default]"
    ),
    
    make_option(c("-y", "--exclude_sample_file"), 
                type="character", default=NULL, 
                help=paste0(
                  "Samples need to be excluded; a column of sample id without colname ",
                  "[default %default]"
                )
    ),
    
    make_option(c("-z", "--exclude_gene_file"), 
                type="character", default=NULL, 
                help=paste0(
                  "Genes need to be excluded; a column of gene id without colname ",
                  "[default %default]"
                )
    )
  )
  
  parser = OptionParser(
    usage='%prog', 
    option_list=optionList,
    description = paste0(
      'Runs PEER'
    )
  )
  
  # a hack to fix a bug in optparse that won't let you use positional args
  # if you also have non-boolean optional args:
  getOptionStrings = function(parserObj) {
    optionStrings = character()
    for(item in parserObj@options) {
      optionStrings = append(optionStrings, 
                             c(item@short_flag, item@long_flag))
    }
    optionStrings
  }
  optStrings = getOptionStrings(parser)
  arguments = parse_args(parser, positional_arguments=TRUE)
  
  # set the seed
  set.seed(arguments$options$seed)
  
  # read in the matrix
  f=arguments$options$file
  if (f == 'stdin') {
    dat = read.table(file=file('stdin'), header=TRUE, 
                     stringsAsFactors=FALSE, check.names=FALSE)
  } else {
    dat = data.frame(fread(paste("zcat", f), header=TRUE, 
                           stringsAsFactors=FALSE, showProgress=FALSE, check.names=FALSE), 
                     check.names=FALSE)
  }
  # fix case if input file has #chr in header
  colnames(dat) = gsub("X.ch", "ch", colnames(dat))
  colnames(dat) = gsub("#ch", "ch", colnames(dat))
  
  # subset columns
  smpl_start = arguments$options$sample_column_start
  sample_ids = colnames(dat)[smpl_start:ncol(dat)]
  other_cols = colnames(dat)[1:smpl_start-1]
  
  if (!is.null(arguments$options$exclude_sample_file)){ 
    sample_ex = read.table(arguments$options$exclude_sample_file, sep="\t", header=FALSE, stringsAsFactors=FALSE,check.names=FALSE)$V1
    cat('Excluded samples are ', print(sample_ex), '\n')
    sample_ids=sample_ids[!sample_ids%in%sample_ex]
    }
  
  cat('length(sample_ids) is ', str(sample_ids), '\n')
  # feature id
  feature_ids_colname = as.character(arguments$options$feature_ids_colnam)
  print(feature_ids_colname)
  print(dat[1:6,1:10])
  rownames(dat) = dat[,feature_ids_colname]
  
  
  #exclude genes
  if (!is.null(arguments$options$exclude_gene_file)) { 
   exclude_gene = read.table(arguments$options$exclude_gene_file, sep="\t", header=FALSE, stringsAsFactors=FALSE,check.names=FALSE)$V1
   cat('Excluded genes are ', print(exclude_gene), '\n')
   dat = subset(dat, !rownames(dat) %in% exclude_gene)
  }
  
  cat('length(sample_ids) is ', str(sample_ids), '\n')
  cat('After subseting genes dim dat is', dim(dat), '\n')
  dat_mtx = dat[,sample_ids]
  cat('After subseting dat is', dim(dat_mtx), '\n')
  feature_ids = rownames(dat_mtx)
  cat('After subseting dat, feature id is', head(feature_ids), '\n')
  dat_mtx = t(dat_mtx)
  dat_other = dat[,other_cols]
  
  
  # read in covariates
  if (!is.null(arguments$options$covariate_file)) {
    f2=arguments$options$covariate_file
    covariates_mtx = data.frame(fread(paste("zcat", f2), sep="\t", 
                                      header=TRUE, stringsAsFactors=FALSE, 
                                      showProgress=FALSE, check.names=FALSE), check.names=FALSE)
    # set the row names to covariate ids
    rownames(covariates_mtx) = covariates_mtx[[colnames(covariates_mtx[1])]]
    covariates_mtx[colnames(covariates_mtx[1])] = NULL
    if (!all(sample_ids %in% colnames(covariates_mtx))) {
      cat("ERROR: sample ids dont match between covariates and mol trait.\n")
      print(sample_ids[!(sample_ids %in% colnames(covariates_mtx))])
      stop()
    }
    covariates_mtx = t(as.matrix(covariates_mtx[,sample_ids]))
  } else {
    covariates_mtx=NULL
  }
  
  cat('dat_mtx is ', dat_mtx[1:3,1:3], '\n')
  # run peer
  peer_model(dat_mtx,
             arguments$options$hidden_factors,
             arguments$options$iterations,
             arguments$options$output_prefix,
             dat_other,
             covariates_mtx,
             arguments$options$account_mean,
             arguments$options$inverse_norm,
             feature_ids)
}

run_peer()
################################################################################


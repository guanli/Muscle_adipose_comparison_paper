#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

######################## Required Packages #####################################
suppressMessages(library(ggplot2))
################################################################################

######################## Functions #############################################
run = function() {
  suppressPackageStartupMessages(library(optparse))
  
  optionList = list(
    make_option(c("-f", "--file"), 
      type="character", 
      help=paste0(
        "Input tab delimited file with the following columns: ",
        "(num_peer\tnum_eGene)."
      )
    ),
    
    make_option(c("-y", "--y_axis_label"), 
      type="character", default="Number of molecular traits with a QTL", 
      help="Label for y axis. [default %default]"
    ),
    
    make_option(c("-o", "--out_file"), 
      type="character", default="peer-optimize_perm", 
      help="Name of output file (will be in pdf format)"
    )
  )
  
  parser = OptionParser(
      usage='%prog', 
      option_list=optionList,
      description = paste0(
          'Plots summary of QTLs per PEER factor'
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
    
  # read in the matrix
  f=arguments$options$file
  if (f == 'stdin') {
    dat = read.table(file=file('stdin'), sep='\t', header=TRUE, 
      stringsAsFactors=FALSE, check.names=FALSE)
  } else {
    dat = read.table(file=f, header=TRUE, sep='\t')
    dat$num_eGene = dat$num_eGene-1
  }
    
  base=arguments$options$out_file
  y_lab=arguments$options$y_axis_label

  pdf(file=paste(base, ".pdf", sep = ""), height=4.5, width=5)
    plt = ggplot(dat, aes(x=num_peer, y=num_eGene)) 
    plt = plt + theme_bw() 
    plt = plt + geom_point(alpha=0.4) + geom_line(alpha=0.7)
    plt = plt + labs(x="PEER factors", y=y_lab)
    print(plt)
  dev.off()
  return(0)  
}
################################################################################


run()

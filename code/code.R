args <- commandArgs(trailingOnly = TRUE)

valid.args <- args %in% c("no-output", "no-figs", "no-examples", 
                          "incl-appendix", "appendix-only", "out-to-file", "test")
if (!all(valid.args)) {
  stop("Command line argument ",
       min(which(!valid.args)),
       " not recognized.")
}

if ("out-to-file" %in% args) {
  sink("./code_output.txt")
}

test <- "test" %in% args
do_main <- !("appendix-only" %in% args)
do_appendix <- ("incl-appendix" %in% args) || ("appendix-only" %in% args)
do_examples <- !(("no-examples" %in% args) || ("appendix-only" %in% args))

library(tidyverse)
library(ashr)
library(ebnm)
library(EbayesThresh)
library(microbenchmark)
library(xtable)


# Output:

if (!("no-output" %in% args)) {
  setwd("./make_output/")
  
  # Output needed for figures in main text:
  if (do_main) {
    source("./ashr_grid_smn.R")
    source("./ashr_grid_symm.R")
    source("./time_comps.R")
    source("./wOBA.R")
  }
  
  # Output needed for figures in appendices:
  if (do_appendix) {
    source("./optmethods.R")
    source("./ebayesthresh.R")
    source("./ebnm_v_ashr.R")
    source("./optgrids.R")
  }
  
  setwd("../")
}


# Figures:

if (!("no-figs" %in% args)) {
  cat("\nGenerating figures...\n\n")
  
  setwd("./make_figs/")
  
  if (do_main) {
    source("./ashr_grid_figs.R")
    source("./timecomps_fig.R")
  }
  
  if (do_appendix) {
    source("./optmethods_figs.R")
    source("./ebayesthresh_fig.R")
    source("./ebnm_v_ashr_fig.R")
    source("./optgrids_fig.R")
  }
  
  setwd("../")
}


# Examples:

if (do_examples) {
  setwd("./examples/")
  
  source("./sims_example.R")
  source("./wOBA_example.R")
  
  setwd("../")
}


cat("\n\n")
sessionInfo()
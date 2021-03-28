args <- commandArgs(trailingOnly = TRUE)

valid.args <- args %in% c("no-output", "no-figs", "no-examples", "incl-appendix", "out-to-file", "test")
if (!all(valid.args)) {
  stop("Command line argument ",
       min(which(!valid.args)),
       " not recognized.")
}

if ("out-to-file" %in% args) {
  sink("./code_output.txt")
}

test <- "test" %in% args

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
  source("./ashr_grid_smn.R")
  source("./ashr_grid_symm.R")
  source("./time_comps.R")
  source("./wOBA.R")
  
  # Output needed for figures in appendices:
  if ("incl-appendix" %in% args) {
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
  
  source("./ashr_grid_figs.R")
  source("./timecomps_fig.R")
  
  if ("incl-appendix" %in% args) {
    source("./optmethods_figs.R")
    source("./ebayesthresh_fig.R")
    source("./ebnm_v_ashr_fig.R")
    source("./optgrids_fig.R")
  }
  
  setwd("../")
}


# Examples:

if (!("no-examples" %in% args)) {
  setwd("./examples/")
  
  source("./sims_example.R")
  source("./wOBA_example.R")
  
  setwd("../")
}


cat("\n\n")
sessionInfo()
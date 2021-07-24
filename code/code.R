args <- commandArgs(trailingOnly = TRUE)

valid.args <- args %in% c(
  "no-output", 
  "no-figs", 
  "incl-appendix", 
  "appendix-only", 
  "out-to-file", 
  "test"
)
if (!all(valid.args)) {
  stop("Command line argument ",
       min(which(!valid.args)),
       " not recognized.")
}

if ("out-to-file" %in% args) {
  if (!("incl-appendix" %in% args)) {
    sink("./code_output_main.txt")
  } else if ("appendix-only" %in% args) {
    sink("./code_output_appendix.txt")
  } else {
    sink("./code_output_all.txt")
  }
}

test <- "test" %in% args
do_main <- !("appendix-only" %in% args)
do_appendix <- ("incl-appendix" %in% args) || ("appendix-only" %in% args)


library(tidyverse)
library(ashr)
library(ebnm)
library(EbayesThresh)
library(REBayes)
library(microbenchmark)
library(xtable)
library(gt)
library(ggpubr)


# Output:

if (!("no-output" %in% args)) {
  cat("\nGenerating output...\n\n")
  
  setwd("./make_output/")
  
  # Output needed for figures in main text:
  if (do_main) {
    source("./time_comps.R")
    source("./gridtests.R")
  }
  
  # Output needed for figures in appendices:
  if (do_appendix) {
    source("./optmethods.R")
    source("./ebayesthresh.R")
    source("./ebnm_v_ashr.R")
    source("./rebayes.R")
    source("./wOBA.R")
  }
  
  setwd("../")
}


# Figures:

if (!("no-figs" %in% args)) {
  cat("\nGenerating figures...\n\n")
  
  setwd("./make_figs/")
  
  if (do_main) {
    source("./smnKLdiv_fig.R")
    source("./timecomps_fig.R")
    source("./gridtests.R")
  }
  
  if (do_appendix) {
    source("./optmethods_figs.R")
    source("./ebayesthresh_fig.R")
    source("./ebnm_v_ashr_fig.R")
    source("./rebayes_fig.R")
  }
  
  setwd("../")
}


# Examples:

setwd("./examples/")

if (do_main) {
  source("./sims_study.R")
  source("./eight_schools.R")
  source("./gtex.R")
}

if (do_appendix) {
  source("./wOBA_example.R")
}

setwd("../")


cat("\n\n")
sessionInfo()

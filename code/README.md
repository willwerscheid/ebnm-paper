# Code

1. Install the necessary packages:

install.packages(c("tidyverse", "ashr", "EbayesThresh", "microbenchmark", "xtable", "devtools"))
devtools::install_github("stephenslab/ebnm")

2. Run code.R from this directory. Command-line options include:
  no-output: omit output
  no-figs: omit figures
  no-examples: omit examples
  incl-appendix: also generate output/figures from the appendix
  test: changes parameters so that the code runs really quickly

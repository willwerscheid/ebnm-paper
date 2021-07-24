# Code

1. Install the necessary packages:

install.packages(c("tidyverse", "ashr", "EbayesThresh", "REBayes", "microbenchmark", "xtable", "gt", "ggpubr", "devtools"))
devtools::install_github("stephenslab/ebnm")
devtools::install_github("willwerscheid/flashier")

2. (Appendix only): Install MOSEK and acquire MOSEK license (http://www.mosek.com).

3. Run code.R from this directory. Command-line options include:
  no-output: omit output
  no-figs: omit figures
  incl-appendix: also generate output/figures from the appendix
  appendix-only: only generate output/figures from the appendix
  out-to-file: print output (including progress updates and session Info) to file
  test: changes parameters so that the code runs really quickly

# Code

1. Install the necessary packages:

install.packages(c("tidyverse", "ashr", "EbayesThresh", "REBayes", "microbenchmark", "xtable", "gt", "ggpubr", "devtools", "latex2exp", "webshot2"))
devtools::install_github("stephenslab/ebnm")
devtools::install_github("willwerscheid/flashier")

2. (Appendix only): Install MOSEK and acquire MOSEK license (http://www.mosek.com). Then install Rmosek within R as follows:

install.packages("Rmosek")
Rmosek::mosek_attachbuilder("<MSKHOME>/mosek/<MSKVERSION>/tools/platform/<PLATFORM>/bin")
install.rmosek()

3. Run code.R from this directory. Command-line options include:
  no-output: omit output
  no-figs: omit figures
  incl-appendix: also generate output/figures from the appendix
  appendix-only: only generate output/figures from the appendix
  out-to-file: print output (including progress updates and session Info) to file
  test: changes parameters so that the code runs quicker

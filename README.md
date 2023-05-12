Hi! Welcome to the first Peach Gene Coexpression Network (PeachGCN).

In this website, you will find the data and scripts needed to perform a candidate gene analysis (CGA) using the PeachGCN. If you are not familiar with coding languages, don't worry, the PeachGCN can be used for everyone, even those with no experience in coding. All you have to do is follow the simple guidelines presented in this readme.  The scripts presented here are coded in python and R, so you will need to have R and python installed in your computer. You can download them from here: https://www.python.org/downloads/ and https://cran.r-project.org/bin/windows/base/. Also, you will need to set up some variables in those scripts, so you need a coding enviroment. You can use spyder for python and Rstudio for R, they are user-friendly and easy to install: https://www.spyder-ide.org/, https://posit.co/products/open-source/rstudio/.

The CGA is divided in two steps:
  - First, you will extract from the PeachGCN the neighbors (coexpressing genes) of your gene or genes of interest. This step is performed with a super simple python script named neighbors.py.
  - Second, you will perform a enrichment analysis using the neighbors of your gene or genes of interest using clusterProfiler. This step is even more simple, as you will run it using an R script named enrichment_analysis.R.

With enrichment_analysis.R you will use Gene Onthology and Mapman as onthologies to perform your candidate gene analysis. Last version of Gene Onthology can be accessed through biomaRt, but for Mapman you will have to use a downloaded version. You can find MapMan Pathways version 4.2 (Thimm et al., 2004) in this repository. 


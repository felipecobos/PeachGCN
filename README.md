Hi! Welcome to the first Peach Gene Coexpression Network (PeachGCNv1).

In this repository, you will find the data and scripts needed to perform a candidate gene analysis (CGA) using the PeachGCNv1. If you are not familiar with coding languages, don't worry, the PeachGCNv1 can be used for everyone, even those with no experience in coding. All you have to do is follow the simple guidelines presented in this readme.

The first thing you will have to do is to download the PeachGCNv1. For a matter of memory, it is compressed in rar format. So after download it you will need to decompress the file. Once downloaded and decompressed, you will be able to use it. We provided the PeachGCNv1 to all the scientific community as a tool for overcoming the intrinsic limitations of working with crop tree species, prioritize research lines and outline new ones. It can be used in countless ways, but here we provided to our colleges the scripts needed to perform a CGA as it is done in Pérez de los Cobos et al., 2023.

The scripts presented here are coded in python and R, so you will need to have R and python installed in your computer. You can download them from here: https://www.python.org/downloads/ and https://cran.r-project.org/bin/windows/base/. Also, you will need to set up some variables in those scripts, so you will need a coding enviroment. You can use spyder for python and Rstudio for R, they are user-friendly and easy to install: https://www.spyder-ide.org/, https://posit.co/products/open-source/rstudio/.

The CGA is divided in two steps:
  - First, you will extract from the PeachGCNv1 the neighbors (coexpressing genes) of your gene or genes of interest. This step is performed with a super simple python script named Neighbors.py.
  - Second, you will perform a enrichment analysis using the neighbors of your gene or genes of interest using clusterProfiler. This step is even more simple, as you will run it using an R script named Enrichment_analysis.R.

With neighbors.py you will use two python libraries called pandas and itertools. Itertools is a built-in module in python and does not need to be installed. However, pandas need to be installed. You can download it from here: https://sparkbyexamples.com/pandas/install-pandas-on-windows/.

With Enrichment_analysis.R you will use Gene Onthology and Mapman as onthologies to perform your candidate gene analysis. With this script you will use five R packages called biomaRt, ClusterProfiler, dplyr, readxl and ggplot2. All these packages can be installed using the commands at the beginning of the script. Last version of Gene Onthology can be accessed through biomaRt, but for Mapman you will have to use a downloaded version. You can find MapMan Pathways version 4.2 (Thimm et al., 2004) in this repository as X4.2_prunus_persica.xlsx. Additionaly, we provide the annotation of the peach genome using MapMan as MapMan_annotation.xlsx.

References:
Pérez de los Cobos, F., García-Gómez, B. E., Orduña-Rubio, L., Batlle, I., Arús, P., Matus, J. T., & Eduardo, I. (2023). First large-scale peach gene coexpression network: A new tool for predicting gene function. bioRxiv.

Thimm, O., Bläsing, O., Gibon, Y., Nagel, A., Meyer, S., Krüger, P., Selbig, J., Müller, L. A., Rhee, S. Y., & Stitt, M. (2004). mapman: a user-driven tool to display genomics data sets onto diagrams of metabolic pathways and other biological processes. The Plant Journal, 37(6), 914–939. 

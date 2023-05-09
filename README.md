Hi! Welcome to the first Peach Gene Coexpression Network (PeachGCN).

In this website, you will find the data and scripts needed to perform a candidate gene analysis (CGA) using the PeachGCN. If you are not familiar with coding languages, don't worry, the PeachGCN can be used for everyone, even those with no experience in coding. All you have to do is follow the simple guidelines presented in this readme. The CGA in divided in two steps:
  - First, you will extract from the PeachGCN the neighbors (coexpressing genes) of your gene or genes of interest. This step is performed with a super simple python script named neighbors.py.
  - Second, you will perform a enrichment analysis using the neighbors of your gene or genes of interest using clusterProfiler. This step is even more simple, as you will run it using an R script named enrichment_analysis.R.
 
 

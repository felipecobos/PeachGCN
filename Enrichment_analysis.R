#Installing required packages
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("biomaRt")
#BiocManager::install("clusterProfiler")
#install.packages("dplyr")
#install.packages("readxl")
#install.packages("ggplot2")

library(biomaRt)
library(clusterProfiler)
library(dplyr)
library(readxl)
library(ggplot2)

###
#Enrichment analysis with clusterProfiler
###

#Setting variables
output_path = "path/to/folder"
genes_path = "path/to/genes/GCNs"
mapman_path = "path/to/MapMan/onthology"
mapman_name_path = "path/to/MapMan/names"
gene_list = c("PRUPE_5G006200")

m <- useMart("plants_mart", dataset="ppersica_eg_gene", host="plants.ensembl.org")
setwd(output_path)

for (gene in gene_list){
  
  neighboors = read.csv(paste0(GCN_path,"/",gene,"_neighborhood.csv"), header = TRUE)
  
  GO_MF <- getBM(attributes=c("go_id", "ensembl_gene_id", "namespace_1003"), mart=m)
  GO_MF = subset(GO_MF, GO_MF$namespace_1003 == "molecular_function")
  GO_MF = data.frame(go_id=GO_MF$go_id, ensembl_gene_id=GO_MF$ensembl_gene_id)
  
  GO_MF_name = getBM(attributes=c("go_id", "name_1006"), mart=m)
  GO_MF_name = subset(GO_MF_name, GO_MF_name$go_id != "")
  
  GO_MF_network =  getBM(attributes=c("ensembl_gene_id", "go_id", "namespace_1003", "name_1006"), mart=m, filters = "ensembl_gene_id", values = neighboors$genes)
  GO_MF_network = subset(GO_MF_network, GO_MF_network$namespace_1003 == "molecular_function")
  GO_MF_network_unique = unique(GO_MF_network$go_id)
  
  id_reps = list()
  for (GO_MF_term in GO_MF_network_unique) {
    reps = length(subset(GO_MF_network, GO_MF_network$go_id == GO_MF_term)[,1])
    id_reps = paste0("id_reps = c(id_reps, list(",GO_MF_term,"=reps))")
  }
  
  vector_reps = c()
  for (GO_MF_term in GO_MF_network$go_id){
    vector_reps = paste0("c(vector_reps, id_reps[[",GO_MF_term,"]])")
  }
  
  GO_MF_network$reps = vector_reps
  write.csv(GO_MF_network, "./GO_MF_network.csv", row.names = FALSE)
  
  GO_CC <- getBM(attributes=c("go_id", "ensembl_gene_id", "namespace_1003"), mart=m)
  GO_CC = subset(GO_CC, GO_CC$namespace_1003 == "cellular_component")
  GO_CC <- data.frame(go_id=GO_CC$go_id, ensembl_gene_id=GO_CC$ensembl_gene_id)
  
  GO_CC_name = getBM(attributes=c("go_id", "name_1006"), mart=m)
  GO_CC_name = subset(GO_CC_name, GO_CC_name$go_id != "")
  
  GO_CC_network =  getBM(attributes=c("ensembl_gene_id", "go_id", "namespace_1003", "name_1006"), mart=m, filters = "ensembl_gene_id", values = neighboors$genes)
  GO_CC_network = subset(GO_CC_network, GO_CC_network$namespace_1003 == "cellular_component")
  GO_CC_network_unique = unique(GO_CC_network$go_id)
  
  id_reps = list()
  for (GO_CC_term in GO_CC_network_unique) {
    reps = length(subset(GO_CC_network, GO_CC_network$go_id == GO_CC_term)[,1])
    id_reps = pasteo("c(id_reps, list(",GO_CC_term,"=reps))")
  }
  
  vector_reps = c()
  for (GO_CC_term in GO_CC_network$go_id){
    vector_reps = paste0("c(vector_reps, id_reps[[",GO_CC_term,"]])")
  }
  
  GO_CC_network$reps = vector_reps
  write.csv(GO_CC_network, "./GO_CC_network.csv", row.names = FALSE)
  
  
  GO_BP <- getBM(attributes=c("go_id", "ensembl_gene_id", "namespace_1003"), mart=m)
  GO_BP = subset(GO_BP, GO_BP$namespace_1003 == "biological_process")
  GO_BP <- data.frame(go_id=GO_BP$go_id, ensembl_gene_id=GO_BP$ensembl_gene_id)
  
  GO_BP_name = getBM(attributes=c("go_id", "name_1006"), mart=m)
  GO_BP_name = subset(GO_BP_name, GO_BP_name$go_id != "")
  
  GO_BP_network =  getBM(attributes=c("ensembl_gene_id", "go_id", "namespace_1003", "name_1006"), mart=m, filters = "ensembl_gene_id", values = neighboors$genes)
  GO_BP_network = subset(GO_BP_network, GO_BP_network$namespace_1003 == "biological_process")
  GO_BP_network_unique = unique(GO_BP_network$go_id)
  
  id_reps = list()
  for (GO_BP_term in GO_BP_network_unique) {
    reps = length(subset(GO_BP_network, GO_BP_network$go_id == GO_BP_term)[,1])
    id_reps = paste0("id_reps = c(id_reps, list(",GO_BP_term,"=reps))")
  }
  
  vector_reps = c()
  for (GO_BP_term in GO_BP_network$go_id){
    vector_reps = paste0("c(vector_reps, id_reps[[",GO_BP_term,"]])")
  }
  
  GO_BP_network$reps = vector_reps
  write.csv(GO_BP_network, "./GO_BP_network.csv", row.names = FALSE)
  
  Mapman <- read.csv(mapman_path,"/Mapman_annotations.xlsx")
  Mapman <- data.frame(bincode = Mapman$bincode, ensembl_gene_id = Mapman$ensembl_gene_id)
  Mapman_name <- read_excel(mapman_name_path,"/X4.2_prunus_persica.xlsx")
  Mapman_name <- Mapman_name[!duplicated(Mapman_name),]
  colnames(Mapman_name)=c("bincode", "Description")
  
  colnames(neighboors)=c("ensembl_gene_id")
  Mapman_network = merge(neighboors, Mapman)
  Mapman_network = merge(Mapman_network, Mapman_name)
  Mapman_network = Mapman_network[,c(2,1,3)]
  Mapman_network_unique = unique(Mapman_network$bincode)
  
  id_reps = list()
  for (Mapman_term in Mapman_network_unique) {
    reps = length(subset(Mapman_network, Mapman_network$bincode == Mapman_term)[,1])
    expression = paste0("id_reps = c(id_reps, list(\"",Mapman_term,"\"=reps))")
    eval(parse(text=expression))
  }
  
  vector_reps = c()
  for (Mapman_term in Mapman_network$bincode){
    expression=paste0("vector_reps = c(vector_reps, id_reps[[\"",Mapman_term,"\"]])")
    eval(parse(text=expression))
  }
  
  Mapman_network$reps = vector_reps
  Mapman_network = subset(Mapman_network, Mapman_network$Description != "not assigned.not annotated")
  Mapman_network = subset(Mapman_network, Mapman_network$Description != "not assigned.annotated")
  write.csv(Mapman_network, "./Mapman_network.csv", row.names = FALSE)
  colnames(neighboors)=c("genes")
  
  genes_anotados = c(length(unique(GO_BP_network$ensembl_gene_id)), length(unique(GO_MF_network$ensembl_gene_id)), length(unique(GO_CC_network$ensembl_gene_id)), length(unique(Mapman_network$ensembl_gene_id)))
  genes_anotados = data.frame(Genes_Anotados = genes_anotados, Dataset = c("GObp", "GOmf", "GOcc", "Mapman"))
  

  entrez_annotation <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"), filters = "ensembl_gene_id", values = neighboors$genes,  mart=m)
  metadata <- getBM(attributes=c("ensembl_gene_id", "description"), filters = "ensembl_gene_id", values = neighboors$genes,  mart=m)
  write.csv(metadata, "./Network_metadata_15072021.csv", row.names = FALSE)
  
  GO_MF_enriched = enricher(neighboors$genes, pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 2, maxGSSize = 100000, qvalueCutoff = 0.05, TERM2GENE = GO_MF, TERM2NAME = GO_MF_name)
  GO_MF_enriched_table <- data.frame(GO_MF_enriched)
  GO_MF_enriched_table$Dataset = rep("GOmf", nrow(GO_MF_enriched_table))
  #GO_MF_enriched_table = subset(GO_MF_enriched_table, GO_MF_enriched_table$Count > 1)
  write.csv(GO_MF_enriched_table, "./Tables/GO_MF_enriched.csv", row.names = FALSE)
  
  GO_CC_enriched = enricher(neighboors$genes, pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 2, maxGSSize = 100000, qvalueCutoff = 0.05, TERM2GENE = GO_CC, TERM2NAME = GO_CC_name)
  GO_CC_enriched_table <- data.frame(GO_CC_enriched)
  GO_CC_enriched_table$Dataset = rep("GOcc", nrow(GO_CC_enriched_table))
  #GO_CC_enriched_table = subset(GO_CC_enriched_table, GO_CC_enriched_table$Count > 1)
  write.csv(GO_CC_enriched_table, "./Tables/GO_CC_enriched.csv", row.names = FALSE)

  GO_BP_enriched = enricher(neighboors$genes, pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 2, maxGSSize = 100000, qvalueCutoff = 0.05, TERM2GENE = GO_BP, TERM2NAME = GO_BP_name)
  GO_BP_enriched_table <- data.frame(GO_BP_enriched)
  GO_BP_enriched_table$Dataset = rep("GObp", nrow(GO_BP_enriched_table))
  #GO_BP_enriched_table = subset(GO_BP_enriched_table, GO_BP_enriched_table$Count > 1)
  write.csv(GO_BP_enriched_table, "./Tables/GO_BP_enriched.csv", row.names = FALSE)

  Mapman_enriched <- enricher(neighboors$genes, pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 2, maxGSSize = 100000, qvalueCutoff = 0.05, TERM2GENE = Mapman, TERM2NAME = Mapman_name)
  Mapman_enriched_table <- data.frame(Mapman_enriched)
  Mapman_enriched_table$Dataset = rep("MapMan", nrow(Mapman_enriched_table))
  write.csv(Mapman_enriched_table, "./Tables/Mapman_enriched.csv", row.names = FALSE)
  
  enriched_table <- rbind(GO_BP_enriched_table, GO_MF_enriched_table, GO_CC_enriched_table, Mapman_enriched_table)
  }

enriched_table_plot = subset(enriched_table, enriched_table$Description != "")
enriched_table_plot$Dataset[enriched_table_plot$Dataset == "goBP"] = "GObp"
enriched_table_plot$Dataset[enriched_table_plot$Dataset == "goCC"] = "GOcc"
enriched_table_plot$Dataset[enriched_table_plot$Dataset == "goMF"] = "GOmf"
enriched_table_plot$Dataset[enriched_table_plot$Dataset == "Mapman"] = "MapMan"
enriched_table_plot = subset(enriched_table_plot, enriched_table_plot$Dataset == "GObp" | enriched_table_plot$Dataset == "GOcc" | enriched_table_plot$Dataset == "GOmf")
enriched_table_plot$Dataset = factor(enriched_table_plot$Dataset, levels = c("GObp", "GOmf", "GOcc", "MapMan"))
row.names(enriched_table_plot) = NULL

i=1
for (Description in enriched_table_plot$Description){
  
  if (nchar(Description) >= 60){
    Description = paste(substr(Description, 1, 60), "...", sep="", collapse = NULL)
    print(Description)
    
    enriched_table_plot$Description[i] = Description
  }
  i=i+1
}

ggplot(enriched_table_plot, aes(y = Count, x = reorder(Description, -Count), fill = Dataset)) + geom_segment(aes(y = 0, yend = Count, reorder(Description, -Count), xend = Description, color = Dataset), size = 1) + geom_point(aes(color=Dataset, size = -qvalue)) + theme_classic() + theme(axis.text.x = element_text(size = 8, angle = 270, vjust = 1, hjust = 0)) + xlab("Term Description") + ylab("Number of genes") + scale_color_manual(values = c("#5AA4CD", "#DE4A6A", "#73AC31")) + scale_size_continuous(limits=c(-0.05, -0),breaks=c(0,-0.01,-0.02,-0.03,-0.04,-0.05), labels=c("0", "0.01", "0.02", "0.03", "0.04", "0.05"), name = "qvalue")

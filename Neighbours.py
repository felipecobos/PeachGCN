#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 09:21:35 2020

@author: felipe
"""

"""
Importing libraries
"""
import pandas as pd
import itertools as it

"""
Declaring functions
"""
def calculate_neighbours(gene, network):
    neighbours = []
    gene_df = pd.DataFrame(network[gene])
    for neighbour in list(gene_df.index.values) :
        if gene_df.loc[neighbour][gene] > 0:
            neighbours.append(neighbour)
    
    return neighbours 


def unique(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    
    return unique_list



"""
Importing the PeachGCN
"""
#Example: Binary_matrix = pd.read_table("/home/felipe/Binary_matrix_300_nozeros_filtered_08062021.csv", sep=",", header=0, index_col=0)
Binary_matrix = pd.read_table("", sep=",", header=0, index_col=0)


"""
Declaring variables
"""
#Example: candidate_genes = ["PRUPE_4G261900", "PRUPE_4G262200"]
candidate_genes = ["", ""]

#Example: ouput_path = "home/felipe/"
ouput_path = ""


"""
This loop checks if the candidate gene is in the PeachGCN and exports a csv with the candidate gene's neighborhood
"""
for candidate_gene in candidate_genes:
    
    if not candidate_gene in Binary_matrix.columns :
        print(candidate_gene + " is not in the GCN\n")
        continue
    
    candidate_gene_df = pd.DataFrame(Binary_matrix[candidate_gene])
    
    gene_neighborhood = calculate_neighbours(candidate_gene, Binary_matrix)
    
    genes = [candidate_gene] + gene_neighborhood
    gene_neighborhood_matrix = pd.DataFrame(index = genes, columns = genes)
    gene_neighborhood_matrix = gene_neighborhood_matrix.fillna(0)
    
    combinaciones = list(it.combinations(genes, 2))
    
    for couple in combinaciones:
        
        if Binary_matrix.loc[couple[0]][couple[1]] == 1:
            
            gene_neighborhood_matrix.loc[couple[0]][couple[1]] = 1
            gene_neighborhood_matrix.loc[couple[1]][couple[0]] = 1
    
    output = ouput_path + candidate_gene + "_neighborhood_matrix.csv"
    gene_neighborhood_matrix.to_csv(output, index=True)
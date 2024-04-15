# Non Negative Matrix Factorization to find disease subtypes

This repository contains R script to perform Non Negative Matrix Factorization (NMF). The input for this script requires .csv file where rows contain the genomic events and the columns contain the patients. The matrix contains discrete values such as 0 (event absent) and 1 (event present). The ouput consists of a pdf file showing the subtypes at a chosen membership threshold. 

- Input: /consensus_clustering/input_data/WES_iteration_145.txt
- Output: Membership file, event frequency calculation, fisher test  calculation, pdf for final clustering

# How to run the project
After cloning the repository, cd into consensus clustering folder and run the NMF_Consensus_Clustering.R in R Studio. 

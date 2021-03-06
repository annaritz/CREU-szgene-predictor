# CREU-szgene-predictor

This is the code base for the following paper:

M. Bern, A. King, D. A. Applewhite, and A. Ritz. [Network-based prediction of polygenic disease genes involved in cell motility](https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-019-2834-1). _BMC_Bioinformatics_ 20: 313 (2019). 

It was an outcome of a [Collaborative REU at Reed College supported by CRA-W](https://cra.org/cra-w/creu/), and was presented at the Fifth International Workshop on Computational Network Biology: Modeling, Analysis and Control ([CNB-MAC 2018](https://cnbmac.org/)) in Washington, D.C., in 2018. 

## Network-Based Prediction of Polygenic Disease Genes Involved in Cell Motility
For a disease and a phenotype of interest, this method identifies genes that are functionally associated with both the disease and the phenotype. This method incorporates the functional interactome GIANT to assess functional similarity between genes. We assign scores to genes based on their functional association with user-curated genes known to be involved with the disease and then assign scores to genes based on their functional association with user-curated genes known to be involved with the phenotype. The output is a ranked list of genes based on their involvement with the disease.

## Installation Instructions
* CREU-szgene-predictor was tested on Python 3.7 and requires the following python packages:
    * NetworkX 2.2
    * scipy
    * matplotlib

## Input Files
* ‘-g ‘OR ‘--interaction_graph’: A tab-delimited file that represents the undirected, weighted interactome. Each line has to have at least 3 columns: gene 1, gene 2, and weight. The weight represents the confidence that the two genes are functionally interacting based on a naive bates classifier. This data can be pulled from GIANT.
* '-b' OR ’--biological_process_positives': A file that gives the entrez ID numbers of the genes involved with the cellular phenotypes in question.
* '-d'OR '--disease_positives': A file that gives the entrez ID numbers of the genes involved with the disease in question.
* '-n' OR '--negatives': A file that gives the entrez ID numbers of the genes NOT involved with the disease in question.

## Main Arguments

* '--sinksource_method': Run the sinksource+ algorithm instead of sinksource. This attaches a negatively labeled node to every single other unlabeled node, which reduces the value of nodes that have very few edges.
* '-c' OR '--sinksource_constant': When run with sinksource+, change the value of the edges connecting all the nodes to the sink negative.
* '--single': Run single experiments for disease and biological process and computes score. Also generates figures of results and formats a LaTeX table.
* '-out': Changes the outfile prefix

Other arguments can be found with -help.

## Output Files
* prefix_disease_output.txt: scores of unlabeled genes based on association with disease genes
* prefix_process_output.txt: scores of unlabeled genes based on association with cellular phenotype genes
* prefix_combined_output.txt: scores of unlabeled genes based on association with disease genes and cellular phenotype genes
* prefix_disease_stats.txt: time between iterations when calculating scores of unlabeled genes based on association with disease genes
* prefix_process_stats.txt: time between iterations when calculating scores of unlabeled genes based on association with cellular phenotype genes

## Toy Example

python3 predict.py --single --sinksource_method -c 1 -o schizo_motility

does a single run with the sinksource+ method with a constant of 1. It changes the outfile prefix to -o schizo_motility

## GIANT Network
The GIANT networks are tissue specific functional interactomes built from a naive bayes classifier. The classifier uses protein interaction, expression, and coexpression data to determine the probability that a tissue specific functional interaction occurs between 2 genes.

The HumanBase (GIANT) brain-specific network can be downloaded from: http://hb.flatironinstitute.org/download or http://giant-v2.princeton.edu

## Cell Motility Positive Set

The set of cell motility positives was created by downloading genes from the KEGG database. The genes were downloaded from five cell motility pathways: the cell adhesion molecule (CAM) pathway, focal adhesion kinase (FAK) pathway, ErbB signaling pathway, regulation of actin cytoskeleton pathway, and tight junction pathway.

Kanehisa, Furumichi, M., Tanabe, M., Sato, Y., and Morishima, K.; KEGG: new perspectives on genomes, pathways, diseases and drugs. Nucleic Acids Res. 45, D353-D361 (2017). 

Kanehisa, M., Sato, Y., Kawashima, M., Furumichi, M., and Tanabe, M.; KEGG as a reference resource for gene and protein annotation. Nucleic Acids Res. 44, D457-D462 (2016). 

Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30 (2000). 



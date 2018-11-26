# CREU-szgene-predictor

## Network-Based Prediction of Polygenic Disease Genes Involved in Cell Motility

For a disease and a phenotype of interest, this method identifies genes that are functionally associated with both the disease and the phenotype. This method incorporates the functional interactome GIANT to assess functional similarity between genes. We assign scores to genes based on their functional association with user-curated genes known to be involved with the disease and then assign scores to genes based on their functional association with user-curated genes known to be involved with the phenotype. The output is a ranked list of genes based on their involvement with the disease.
* HumanBase (GIANT) brain-specific network downloaded from: http://hb.flatironinstitute.org/download or http://giant-v2.princeton.edu
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

Other arguments can be found with -help.


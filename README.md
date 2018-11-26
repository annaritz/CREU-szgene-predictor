# CREU-szgene-predictor

## Network-Based Prediction of Polygenic Disease Genes Involved in Cell Motility

For a disease and a phenotype of interest, this method identifies genes that are functionally associated with both the disease and the phenotype. This method incorporates the functional interactome GIANT to assess functional similarity between genes. We assign scores to genes based on their functional association with curated genes known to be involved with the disease and then assign scores to genes based on their functional association with curated genes known to be involved with the phenotype. The output is a ranked list of genes based on their functional involvement with the disease and the phenotype.
* HumanBase (GIANT) brain-specific network downloaded from: http://hb.flatironinstitute.org/download or http://giant-v2.princeton.edu
## Installation Instructions
* CREU-szgene-predictor was tested on Python 3.7 and requires the following python packages:
	* NetworkX 2.2
	* scipy
	* matplotlib

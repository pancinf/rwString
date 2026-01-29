# General Overview
This code conducts a random walk with restart (RWR) analysis based on the diffusr R-package (https://github.com/dirmeier/diffusr) on the STRING V12 database for gene prioritization.

# Required inputs:

## In the "raw" directory:
* 9606.protein.links.v12.0.txt.gz
* 9606.protein.info.v12.0.txt.gz

Both can be downloaded from STRING DB.

## In the "seeds" directory:
* phenoSeeds.txt, containing a one column list of gene symbols. Those genes should represent known causal genes.

Note that "pheno" is a placeholder for your phenotype and is also used as the only required argument for the rwString.R script (so choose the same one!). Just choose a simple term here that describes the studied phenotype/disease best.

# How to use
Rscript --vanilla rwString.R --pheno myPheno --backProb 0.8

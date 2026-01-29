# General Overview
This code conducts a random walk analysis from the diffusr R-package (https://github.com/dirmeier/diffusr) on the STRING V12 database.

# Required inputs:

## In the "raw" directory:
* 9606.protein.links.v12.0.txt.gz
* 9606.protein.info.v12.0.txt.gz

Both can be downloaded from STRING DB.

## In the "seeds" directory:
* phenoSeeds.txt, containing a one column list of gene symbols. Those genes should represent known causal genes.

Note that "pheno" is a placeholder for your phenotype and is also used as the only required argument for the rwString.R script. Just choose a simple term here that describes the studied phenotype/disease best.

# How to use
Rscript --vanilla rwString.R --pheno pancreatitis --backProb 0.8

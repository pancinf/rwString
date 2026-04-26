# General Overview
This code conducts a random walk with restart (RWR) analysis based on the diffusr R-package (https://github.com/dirmeier/diffusr) on the STRING V12 database for gene prioritization.

#  Required R-packages (R v. 4.2.0)
*  optparse v. 1.7.5
*  data.table v. 1.15.4
*  igraph v. 1.4.2
*  diffusr v. 0.1.4 (https://github.com/dirmeier/diffusr)

# Required inputs:

## In the "raw" directory:
* 9606.protein.links.v12.0.txt.gz
* 9606.protein.info.v12.0.txt.gz

Both can be downloaded from STRING DB.

## In the "seeds" directory:
* phenoSeeds.txt, containing a one column list of gene symbols. Those genes should represent known causal genes.

Note that "pheno" is a placeholder for your phenotype and is also used as the only required argument for the rwString.R script (so choose the same one!). Just choose a simple term here that describes the studied phenotype/disease best.

# How to use
* Rscript --vanilla rwString.R --pheno myPheno --backProb 0.8 --norm FALSE

-> This version runs the RWR in permutation mode to compute p-values (slow)
* Rscript --vanilla rwString.R --pheno myPheno --backProb 0.8 --norm TRUE

-> This version runs the RWR in hub-gene-correction mode to compute probabilities (fast)

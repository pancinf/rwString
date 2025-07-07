#This script performs a random walk analysis on the weighted STRING V12 database. Furthermore a permutation test is performed (1000 iterations).

###
###
###Libraries
library(optparse)
library(data.table)
library(igraph)
library(diffusr)

###
###
###Parameters
options(scipen=999)
set.seed(1234)
nPerm <- 1000

###
###
###External arguments
option_list = list(
  make_option(c("-p","--pheno"), type="character", default=NULL, 
              help="Phenotype", metavar="character"),

  make_option(c("-b","--backProb"), type="numeric", default=NULL,
              help="Restart probability", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$pheno) | is.null(opt$backProb)){
  print_help(opt_parser)
  stop("STOP. Please provide a phenotype", call.=FALSE)
}

###
###
###Main

##
##Prepro STRING data
seedGenes <- fread(paste0("../data/seeds/",opt$pheno,"Seeds.txt"), header = FALSE)
seedGenes <- seedGenes$V1
stringInfo <- fread(("../data/raw/9606.protein.info.v12.0.txt.gz"), select = c(1,2))
colnames(stringInfo) <- c("protein","symbol")
stringDb <- fread("../data/raw/9606.protein.links.v12.0.txt.gz")
colnames(stringDb)[3] <- "weight"
stringDb <- as_adjacency_matrix(graph_from_data_frame(stringDb, directed = FALSE),sparse = FALSE,attr = "weight")

##
##Perform random walk
knownSeeds <- stringInfo$protein[stringInfo$symbol %in% seedGenes]
knownSeeds <- sort(knownSeeds)
baseProb <- rep(0,nrow(stringDb))
baseProb[rownames(stringDb) %in% knownSeeds] <- 1/length(knownSeeds)
finalTable <- data.frame(protein = rownames(stringDb),rwStringProb = random.walk(p0 = baseProb,graph = stringDb,correct.for.hubs = FALSE, r = opt$backProb, thresh = 1e-6)$p.inf)
finalTable <- merge(stringInfo,finalTable, by = "protein")
finalTable$rankRwString <- rank(x = -finalTable$rwStringProb, ties.method = "max")

##
##Permutation test
finalTable$lowerThanBackground <- nPerm
usedSeedSets <- c(paste0(sort(knownSeeds),collapse = ", "))
for(u in 1:nPerm){
  novelSeedSet <- FALSE
  
  #Do not consider duplicated seed sets
  while(novelSeedSet == FALSE){
    backSeeds <- sample(x = stringInfo$protein[!(stringInfo$protein %in% knownSeeds)],size = length(knownSeeds),replace = FALSE)
    if(!(paste0(sort(backSeeds),collapse = ", ") %in% usedSeedSets)){
      novelSeedSet <- TRUE  
    }
  }
  usedSeedSets <- c(usedSeedSets,paste0(sort(backSeeds),collapse = ", "))
  baseProb <- rep(0,nrow(stringDb))
  baseProb[rownames(stringDb) %in% backSeeds] <- 1/length(backSeeds)
  rwStringProbBack <- random.walk(p0 = baseProb,graph = stringDb,correct.for.hubs = FALSE,r = opt$backProb, thresh = 1e-6)$p.inf
  finalTable$lowerThanBackground[finalTable$rwStringProb > rwStringProbBack] <- finalTable$lowerThanBackground[finalTable$rwStringProb > rwStringProbBack] - 1
}
finalTable$pValueRwString <- finalTable$lowerThanBackground/nPerm
finalTable$lowerThanBackground <- NULL

##
##Write to output

#rw output
write.table(x = finalTable, file = paste0("../data/output/", opt$pheno,"_backProb",opt$backProb,"_rwString.txt"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)

#used seed sets
write.table(x = usedSeedSets, file = paste0("../data/output/", opt$pheno,"_backProb",opt$backProb,"_rwStringUsedSeedSets.txt"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)

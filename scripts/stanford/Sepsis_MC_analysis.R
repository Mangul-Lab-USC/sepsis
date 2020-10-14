


####################  Run multi-cohort analysis  #############################
# install.packages("rmeta")
# source("http://bioconductor.org/biocLite.R")
# biocLite("multtest")
library(rmeta)
library(multtest)

## Set the path where the file was unzipped
# should contain both scripts and the .RDS file
path <- "/sepsis"
source(file="./Sepsis_MC_analysis_functions.R")

discoveryGEMs <- readRDS(file=paste0(path, "./Sepsis_discovery_GEMs.RDS"))

## All of the discovery cohorts are stored in the list:
# These datasets have already been subsetted according to the criteria listed in the paper (Sweeney et al.). 
names(discoveryGEMs)

## Each item in the list has four objects:
# $expr containts processed probe-level microarray data
# $pheno contains phenotype data, in the same order (by rows) as the columns of $expr
# $keys contains the GEO-assigned probe:gene mappings for the given data set's chip type
# $class contains a binary vector with a 0 for controls (non-infected SIRS/trauma) and a 1 for cases (sepsis)


discoveryResults.Mort = runMetaAnalysis(discoveryGEMs)
discoveryResults$GSEnames <- names(discoveryGEMs)

## View summary statistics
View(discoveryResults.Mort$pooled.ES)

#####################  RESULTS  #########################

threshold = 0.05
pValHet = 0
ESfold = 1.3
scaled <- F

thresh <- analyzeNonLOO(discoveryResults.Mort, threshold, 
                        pValHet, allowedOut=4, ESfold=ESfold)


### regression data -- run separately, see regression file.
out.firthBeckerRand <- readRDS("sepMortRegress_firthBeckerRand.RDS")
# View(out.firthBeckerRand$pooled.ES[order(out.firthBeckerRand$pooled.ES$gene.beta.pchisq), ])
threshReg <- list(pos=rownames(subset(out.firthBeckerRand$pooled.ES, gene.beta.pchisq<0.01 & gene.beta>0)),
                  neg=rownames(subset(out.firthBeckerRand$pooled.ES, gene.beta.pchisq<0.01 & gene.beta<0)) )


thresh$pos <- union(thresh$pos, threshReg$pos)
thresh$neg <- union(thresh$neg, threshReg$neg)

discovery.genes <- convertDiscoveryListToGenes(discoveryResults.Mort, 
                                               thresh$pos, thresh$neg, 
                                               scale.pos=F)  
totalN <- sum(unlist(lapply(discovery.genes, function(x) length(x$class))))


forwardTB <- forwardSearchWeighted(discovery.genes, 
                             thresh$pos, thresh$neg, 
                             yes.pos=c("MPO", "CTSG", "DEFA4"), yes.neg=NULL, 
                             forwardThresh=1, backThresh=totalN/50, scaled=scaled)
getDiscoveryAUCs(discoveryResults.Mort, forwardTB[[1]], forwardTB[[2]], scaled=F)
pos.genes <- forwardTB[[1]]
neg.genes <- forwardTB[[2]]
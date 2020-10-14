

##############################################################
##  This code takes the user through an example of multivariate 
##  meta-regression implemented by Tim Sweeney.
##
##  It will use as an example the MetaIntegrator test object. 
##############################################################

library(MetaIntegrator)
data(tinyMetaObject)
dataList <- tinyMetaObject$originalData


### let's say we wanted to predict group as a function of age, sex, and genes. 
### so the model looks like: glm(group ~ age + sex + gene)
### first we have to create a data.frame within each GSE list named $regDF
### $regDF has a first column named 'response', which here will be the 'group' data. 
### any other columns in this data frame are taken as covariates in meta-regression.

lapply(dataList, function(GSE) colnames(GSE$pheno))

## clean up phenotype names to match
dataList$Whole.Blood.Study.1$pheno$age.years <- dataList$Whole.Blood.Study.1$pheno$age

## notes that in one dataset, 'sex' has only 1 level, (female); in another, there are 3 (M/F/U).
## This will cause failure since categories don't match. so generate random two-level data for the test object. 
dataList$Whole.Blood.Study.1$pheno$sex <- as.numeric(runif(nrow(dataList$Whole.Blood.Study.1$pheno))<0.5)
dataList$PBMC.Study.1$pheno$sex <- as.numeric(runif(nrow(dataList$PBMC.Study.1$pheno))<0.5)

## add the 'regDF' object
dataList <- lapply(dataList, function(GSE) {
    GSE$regDF <- data.frame(response=GSE$pheno$group,
                            age.years=GSE$pheno$age.years,
                            sex=GSE$pheno$sex)
    GSE
})


##############   application of meta-regression   ######################
source("/projects/Sepsis_Tim/scripts/regression_general_metaAnalysis.R")

## a single function exists that does everything
# out.beckerRand <- regressGenMetaAnalysis(dataList) 

## however, it makes more sense to call the two subfunctions individually:
## convertDiscoveryGEMsToGenes() is fairly self-explanatory; converts from 
## probes -> genes, and then ComBat-normalizes between datasets (to allow for 
## same scale in computing betas).
dataListGenes <- convertDiscoveryGEMsToGenes(dataList)

## regressGenMetaAnalysisInner() is the main function. There are four options:
## firth = T/F: Default F. firth regression allows for logistic regression even in the case of perfect separators. 
## becker = T/F: Default T. closed-form solution, includes covariance matrix for all predictors in analysis. 
## intercept = T/F: Default F. DOES NOT calculate regression through intercept if F; instead, simply leaves intercept out of meta-analysis. 
## randomEffects = T/F: Default T. If becker==T, use a fixed or random effect model? 

out.beckerRand <- regressGenMetaAnalysisInner(dataListGenes, firth=F, becker=T, randomEffects=T)

View(out.beckerRand$pooled.ES)

regGenes <- list(pos=rownames(subset(out.beckerRand$pooled.ES, gene.beta.pchisq<0.01 & gene.beta>0)),
                 neg=rownames(subset(out.beckerRand$pooled.ES, gene.beta.pchisq<0.01 & gene.beta<0)) )


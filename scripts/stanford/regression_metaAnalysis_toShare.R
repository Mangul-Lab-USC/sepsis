
## Written by Tim Sweeney, April 2016

###############  data prep for meta-regression  #########################
combatNormalizeGEMsList <- function(gems){
    ## Need to remove any genes which are present but invariant, as prevent Combat. 
    ## may be able to remove if can figure a way to allow Combat to deal with missingness
    gems <- lapply(gems, function(gem) {
        tmp <- apply(gem$genes, 1, function(x) length(unique(x)))
        gem$genes <- gem$genes[tmp!=1, ]
        gem
    })
    
    allgenes <- Reduce(intersect, lapply(gems, function(GEM) rownames(GEM$genes)))
    allgenes <- allgenes[allgenes != ""]
    BigMat <- Reduce(cbind, lapply(gems, function(GEM) GEM$genes[allgenes, ]))
    BigIndex <- Reduce(c, lapply(1:length(gems), function(i) rep(i, dim(gems[[i]]$genes)[2])))     
    
    cat("\nComBat normalizing all GEMs... only common genes included. ")
    require(sva)
    BigMatCB <- ComBat(dat=BigMat, batch=BigIndex)
    
    
    cat("\nCo-normalized gene matrices re-stored separately in study$genes.")
    for(i in 1:length(gems)){
        gems[[i]]$genes <- BigMatCB[, BigIndex == i]
    }
    return(gems)
}


convertGEMtoGeneMtx <- function(GEM){
    uniqueGenes <- GEM$keys[!duplicated(GEM$keys)]
    uniqueGenes <- uniqueGenes[!is.na(uniqueGenes)]
    uniqueGenes <- uniqueGenes[!(uniqueGenes == "")]
    
    singleProbes <- names(table(GEM$keys))[table(GEM$keys)==1]
    multipleProbes <- uniqueGenes[!(uniqueGenes %in% singleProbes)]
    
    singletons <- GEM$keys[GEM$keys %in% singleProbes]
    multiple <- GEM$keys[GEM$keys %in% multipleProbes]
    
    ## get single probes data
    GEM.genes.single <- GEM$expr[names(singletons), ]
    rownames(GEM.genes.single) <- unname(singletons)
    
    if(!length(multipleProbes)>1){
        cat("\nNo multiple-probe genes present\n")
        return(GEM.genes.single[uniqueGenes, ])
        
    } else {
        ## get multiple probes data (simple means)
        GEM.genes.multiple <- t(sapply(multipleProbes, function(gene){
            probes <- names(multiple[multiple %in% gene]);
            colMeans(GEM$expr[probes, , drop=F], na.rm=T);
        }))
        rownames(GEM.genes.multiple) <- multipleProbes
        
        GEM.genes <- rbind(GEM.genes.single, GEM.genes.multiple)
        
        return(GEM.genes[uniqueGenes, ])
    }
}


################  2-var meta MAIN FUNCTION  ##########################
regressGenMetaAnalysisInner <- function(gems, intercept=F, randomEffects=T) {
    ### GEMs passed in are prior converted to genes and quant-normalized
    ### Note: GEMs must have $regDF instead of $class; here $regDF
    ### functions to hold whatever vector is being regressed on (explained by) 
    ### gene expression and thus nrow(GEM$regDF) must equal ncol(GEM$expr).
    ### $regDF is a data.frame with response variable in column 'response' and 
    ### covariates in the other columns
    library(parallel)
    checkGEMregDFs(gems)
    ## Regression effects, combined using inverse variance weighted model
    cat("\nComputing single regression betas...")
    beta.cov.list <- lapply(gems, multiRegressionBecker)
    pooled.ES <- beckerWuGLS(beta.cov.list, intercept=intercept, randomEffects=randomEffects)
    pooled.ES <- pooled.ES[order(pooled.ES$gene.beta.pchisq), ]
    return(list(gems = gems, 
                pooled.ES = pooled.ES,
                GSEnames = names(gems)))
}

checkGEMregDFs <- function(gemList){
    levelMat <- sapply(gemList, function(study){
        levels <- unlist(apply(study$regDF, 2, function(x) length(unique(x))))
        levels[ levels > 4 ] <- 1
        levels
    })
    invisible(lapply(1:nrow(levelMat), function(var){
        varLevels <- levelMat[var, ]
        if(length(unique(varLevels)) > 1){
            warning("Some apparently categorical variables have different numbers of levels:\n", rownames(levelMat)[var], varLevels)
        }
    }))
}

check.metaGEM.study.regress <- function(study){
    stopifnot("genes" %in% names(study)) 
    stopifnot("regDF" %in% names(study))
    stopifnot(typeof(study$regDF) == "list")
    stopifnot(nrow(study$regDF) == ncol(study$expr))
    #stopifnot(any(is.na(study$regDF)))
}



multiRegressionBecker <- function(study){
    check.metaGEM.study.regress(study)
    
    if(length(levels(factor(study$regDF$response)))==2){
        coefs <- apply(study$genes, 1, function(gene){
            tmp <- glm(response~., data=data.frame(study$regDF, gene=gene), family="binomial")
            list(betas.gene = tmp$coefficients, 
                 cov.gene = vcov(tmp))
        })
    } else {
        coefs <- apply(study$genes, 1, function(gene){
            tmp <- glm(response~., data=data.frame(study$regDF, gene=gene), family="gaussian")
            list(betas.gene = tmp$coefficients, 
                 cov.gene = vcov(tmp))
        })
    }
    
    coefs
}

beckerWuGLS <- function(beta.cov.list, intercept=F, randomEffects=T){
    library(MASS)
    library(Matrix)
    ## param for all 
    
    if(intercept){
        betaVec <- beta.cov.list[[1]][[1]]$betas.gene
    } else {
        betaVec <- beta.cov.list[[1]][[1]]$betas.gene[-1]
    }
    P.1 <- length(betaVec)
    k <- length(beta.cov.list)
    W <- Reduce(rbind, lapply(1:k, function(i) diag(1, P.1, P.1)))
    
    if(!randomEffects){
        coefs <- t(sapply(1:length(beta.cov.list[[1]]), function(i){
            if(intercept){
                b <- as.matrix(unlist(lapply(beta.cov.list, function(study) study[[i]]$betas.gene)))
                sigma <- as.matrix(bdiag(lapply(beta.cov.list, function(study) study[[i]]$cov.gene)))
            } else {
                b <- as.matrix(unlist(lapply(beta.cov.list, function(study) study[[i]]$betas.gene[-1])))
                sigma <- as.matrix(bdiag(lapply(beta.cov.list, function(study) study[[i]]$cov.gene[-1,-1])))
            }
            
            ## Becker & Wu, Stat Science 2007, pp9
            cov.beta.hat.star <- ginv(t(W) %*% ginv(sigma) %*% W)
            beta.hat.star <- cov.beta.hat.star %*% t(W) %*% ginv(sigma) %*% b
            
            ## standard for GLS beta p-values
            beta.pchisq <- 1 - pchisq(beta.hat.star^2/diag(cov.beta.hat.star), df=1)
            
            ## to output all values
            out <- data.frame(beta=beta.hat.star, beta.sd=sqrt(diag(cov.beta.hat.star)), beta.pchisq=beta.pchisq)
            names.out <- c(t(outer(names(betaVec), colnames(out), paste, sep = ".")))
            out <- as.vector(t(out))
            names(out) <- names.out
            out
        }))
    } else if (randomEffects){
        coefs <- t(sapply(1:length(beta.cov.list[[1]]), function(i){
            if(intercept){
                betaList <- lapply(beta.cov.list, function(study) study[[i]]$betas.gene)
                sigmaList <- lapply(beta.cov.list, function(study) study[[i]]$cov.gene)
            } else {
                betaList <- lapply(beta.cov.list, function(study) study[[i]]$betas.gene[-1])
                sigmaList <- lapply(beta.cov.list, function(study) study[[i]]$cov.gene[-1,-1])
            }
            
            ## Becker & Wu, Stat Science 2007, pp9
            ## Chen, Manning, Dupuis, Biometrics 2012
            ### start as in fixed effect
            b <- as.matrix(unlist(betaList))
            sigma <- as.matrix(bdiag(sigmaList))
            cov.beta.hat.star <- ginv(t(W) %*% ginv(sigma) %*% W) ## same as Psi
            beta.hat.star <- cov.beta.hat.star %*% t(W) %*% ginv(sigma) %*% b
            
            ### calculate heterogeneity, t.hat
            phi <- ginv(cov.beta.hat.star) - Reduce(`+`, lapply(sigmaList, function(sigmai) ginv(sigmai) %*% cov.beta.hat.star %*% ginv(sigmai) ))
            A <- Reduce(`+`, lapply(1:length(sigmaList), function(j) ginv(sigmaList[[j]]) %*% (betaList[[j]] - beta.hat.star) %*% t((betaList[[j]] - beta.hat.star)) ))
            A <- A - diag(x=length(sigmaList)-1, P.1)
            t.hat <- (ginv(phi) %*% A + t(A) %*% ginv(phi))/2
            
            ### Chen Supplment. Make t.hat positive semi-definite.
            e <- eigen(t.hat)
            t.hat.psd <- e$vectors %*% replace(diag(e$values), diag(e$values)<0, 0) %*% ginv(e$vectors)
            
            ### add t.hat back to block diag vcov mtx
            Omega <- as.matrix(bdiag(lapply(sigmaList, function(sigmai) sigmai + t.hat.psd)))
            
            ## recalc beta/cov with Omega
            cov.beta.hat.R <- ginv(t(W) %*% ginv(Omega) %*% W) 
            beta.hat.R <- cov.beta.hat.R %*% t(W) %*% ginv(Omega) %*% b
            
            ## standard for GLS beta p-values
            beta.pchisq <- 1 - pchisq(beta.hat.R^2/diag(cov.beta.hat.R), df=1)
            
            ## to output all values
            out <- data.frame(beta=beta.hat.R, beta.sd=sqrt(diag(cov.beta.hat.R)), beta.pchisq=beta.pchisq)
            names.out <- c(t(outer(names(betaVec), colnames(out), paste, sep = ".")))
            out <- as.vector(t(out))
            names(out) <- names.out
            out
        }))
    }
    rownames(coefs) <- names(beta.cov.list[[1]])
    gene.beta.FDRchisq = p.adjust(coefs[, "gene.beta.pchisq"], method="BH")
    data.frame(coefs, gene.beta.FDRchisq=gene.beta.FDRchisq)
}

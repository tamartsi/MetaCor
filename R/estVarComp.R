estVarComp <-
function(Y, W, covMatList, IDs = NULL, start = NULL, AIREML.tol = 1e-6, maxIter = 100, dropZeros = TRUE, verbose = TRUE){
		
    
    stopifnot(length(Y) == nrow(W))
    if (is.null(names(Y))) names(Y) <- rownames(W)
        
    k <- ncol(W)
	if (is.null(IDs)) IDs <- rownames(W)
	stopifnot(all(is.element(IDs, rownames(W))))

    # sample size
    n <- length(IDs)
    
    ### subset the data
    W <- W[match(IDs, rownames(W)),]
    Y <- Y[match(IDs, rownames(W))]

      
    # if covMatList is a matrix, convert to a list
    if(class(covMatList) == "matrix"){
        covMatList <- list(A = covMatList)
    }
    
    ## add error group
    group.names <- "E"
    g <- 1

    # number of covariance structure matrices
    m <- length(covMatList)
    
    # if covMatList doesn't have names, assign them
    if(is.null(names(covMatList))){
        names(covMatList) <- paste("A",1:m,sep="")
    }
    
    # check for starting values
    if(!is.null(start)){
            if(length(start) != (m+g)){
                stop("length of start must equal the length of covMatList + number of groups for gaussian")
            }

    }
    
    # subest covariance structure matrices
    for(i in 1:m){
        if(!all(IDs %in% colnames(covMatList[[i]]))){
            stop(paste("All of the included Samples must be in matrix ", i, " of covMatList"))
        }
        # subset matrix
        keepMat <- colnames(covMatList[[i]]) %in% IDs
        covMatList[[i]] <- covMatList[[i]][keepMat,keepMat]
        
        # check that names match
        if(!all(colnames(covMatList[[i]]) == IDs)){
            stop("Column and Row names of matrix ", i, " of covMatList must match the IDs \n (given either directtly as IDs or as rownames of the design matrix W)")
        }
    }
    
    if(verbose) message("Using AIREML Procedure...")
    if(verbose) message(paste("Sample Size: ", n))
    
      # estimate variance components
     if(verbose) message("Computing Variance Component Estimates...")
     if(verbose) message(paste(paste("Sigma^2_",c(names(covMatList),group.names),sep="", collapse="     "), "log-lik", "RSS", sep="     "))
        
      # estimate variance components
      out <- .runAIREMLgaussian(Y=Y, W=W, start=start, m=m, g=g, n=n, k=k, covMatList=covMatList, group.idx=group.idx, AIREML.tol=AIREML.tol, dropZeros=dropZeros, maxIter=maxIter, verbose=verbose)
        
    
    # get estimates and covariance of estimates
    varComp <- out$sigma2.k
    
    
    names(varComp) <- paste("V_",c(names(covMatList),group.names),sep="")
    varCompCov <- matrix(NA, nrow=(m+g), ncol=(m+g))
    colnames(varCompCov) <- paste("V_",c(names(covMatList),group.names),sep="")
    rownames(varCompCov) <- paste("V_",c(names(covMatList),group.names),sep="")
    
    
    if(dropZeros){
        varCompCov[!out$zeroFLAG, !out$zeroFLAG] <- solve(out$AI)
    }else{
        varCompCov <- solve(out$AI)
    }
    
    # compute Cholesky decomposition of Covariance matrix
    if(verbose) message("Computing Cholesky Decomposition of Inverse Covariance Matrix of Phenotype...")
    # Covariance Matrix
     Sigma <- Reduce("+", mapply("*", covMatList, varComp[1:m], SIMPLIFY=FALSE))
     if(g == 1){
         diag(Sigma) <- diag(Sigma) + varComp[m+1]
     }else{
         diagV <- rep(0,n)
         for(i in 1:g){
             diagV[group.idx[[i]]] <- varComp[m+i]
         }
         diag(Sigma) <- diag(Sigma) + diagV
     }
    
    
    # Inverse
    SigmaInv <- chol2inv(chol(Sigma))
    # Cholesky Decomposition
    cholSigmaInv <- t(chol(SigmaInv))
    colnames(cholSigmaInv) <- colnames(covMatList[[1]])
    rownames(cholSigmaInv) <- rownames(covMatList[[1]])
    
    return(list(varComp = varComp, varCompCov = varCompCov, cholSigmaInv = cholSigmaInv, beta = as.vector(out$beta), workingY=as.vector(Y), eta=as.vector(out$eta),  converged = out$converged, zeroFLAG = out$zeroFLAG, logLikR = out$logLikR, logLik = out$logLik, dispersion = out$RSS))
    
    
}

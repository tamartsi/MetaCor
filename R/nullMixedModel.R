nullMixedModel <-
function(Y,
						   W,
                           VC,
                           IDs = NULL,
                           verbose = TRUE){
    
    cholSigmaInv <- VC$cholSigmaInv
    varComp <- VC$varComp

	stopifnot(length(Y) == nrow(W))
    if (is.null(names(Y))) names(Y) <- rownames(W)
        
    k <- ncol(W)
	if (is.null(IDs)) IDs <- rownames(W)
	stopifnot(all(is.element(IDs, rownames(W))))

    k <- ncol(W)

    
    # which samples to remove from cholSigmaInv
    if(!all(IDs %in% colnames(cholSigmaInv))){
        stop("All of the included Samples must be in the cholSigmaInv matrix")
    }
    chol.idx <- which(!(colnames(cholSigmaInv) %in% IDs))
    cholSigmaInv <- .subsetCholSigmaInv(cholSigmaInv, chol.idx)
    
    
    # sample size 
    n <- length(IDs)
    if(verbose) message("Fitting Model with ", n, " Samples")
    
    # calculate fixed effects estimates (beta)
    CW <- crossprod(cholSigmaInv,W)
    CY <- crossprod(cholSigmaInv,Y)
    XtSigInvX <- crossprod(CW)
    XtSigInvXInv <- chol2inv(chol(XtSigInvX))
    beta <- crossprod(XtSigInvXInv, crossprod(CW,CY))
    
    
    # marginal residuals
    fits <- tcrossprod(W, t(beta))
    residM <- as.vector(Y - fits)
    
    # conditional residuals  
    m <- length(varComp)
    residtmp <- crossprod(residM,cholSigmaInv)
    
    residC <- as.vector(varComp[m]*tcrossprod(cholSigmaInv, residtmp))
    
    
    # Variance Covariance of betas
    RSS <- sum(residtmp^2)/(n-k)
    Vbeta <- RSS*XtSigInvXInv
    
    # test statistics and p-values
    SE <- sqrt(diag(Vbeta))
    Stat <- (beta/SE)^2
    pval <- pchisq(Stat, df=1, lower.tail=FALSE)
    
    # find likelihood and AIC
    logLik <- as.numeric(-0.5*n*log(2*pi*RSS) + sum(log(diag(cholSigmaInv))) - 0.5*tcrossprod(residtmp)/RSS)
    AIC <- 2*(k+m) - 2*logLik
    logLikR <- as.numeric( logLik + 0.5*k*log(2*pi*RSS) - 0.5*log(det(XtSigInvX)) )  
    
    # prepare results
    dimnames(Vbeta) <- list(colnames(W), colnames(W))
    
    fixef <- data.frame(Est = beta, SE = SE, Stat = Stat, pval = pval)
    rownames(fixef) <- colnames(W)
    
    return(list(fixef = fixef,
                varComp = varComp,
                resid.marginal = residM,
                resid.conditional = residC,
                logLik = logLik,
                AIC = AIC,
                logLikR = logLikR,
                model.matrix = W,
                Vbeta = Vbeta,
                RSS =RSS))
}

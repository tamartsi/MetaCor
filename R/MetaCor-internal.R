.runAIREMLgaussian <-
function(Y, W, start, m, g, n, k, covMatList, group.idx, AIREML.tol, dropZeros, maxIter, verbose){
    
    # initial values
    sigma2.p <- var(Y)
    AIREML.tol <- AIREML.tol*sigma2.p  # set convergence tolerance dependent on trait
    if(is.null(start)){
        sigma2.k <- rep((1/(m+1))*sigma2.p, (m+g))
    }else{
        sigma2.k <- as.vector(start)
    }
    
    reps <- 0
    repeat({
        reps <- reps+1
        
        zeroFLAG <- sigma2.k < AIREML.tol # which elements have converged to "0"
        sigma2.k[zeroFLAG] <- 0 # set these to 0
        
        # phenotype covariance matrix
        Vre <- Reduce("+", mapply("*", covMatList, sigma2.k[1:m], SIMPLIFY=FALSE))
        
        V <- Vre
        if(g == 1){
            diag(V) <- diag(V) + sigma2.k[m+1]
        }else{
            diagV <- rep(0,n)
            for(i in 1:g){
                diagV[group.idx[[i]]] <- sigma2.k[m+i]
            }
            diag(V) <- diag(V) + diagV
        }
        
        # cholesky decomposition
        cholV <- chol(V)
        # inverse
        Vinv <- chol2inv(cholV)
        VinvW <- crossprod(Vinv,W)
        cholWtVinvW <- chol(crossprod(W, VinvW))
        WtVinvWInv <- chol2inv(cholWtVinvW)
        beta <- crossprod(WtVinvWInv, crossprod(VinvW,Y))
        fits <- tcrossprod(W, t(beta))
        residM <- as.vector(Y - fits)
        VinvR <- crossprod(Vinv, residM)
        RVinvR <- crossprod(residM, VinvR)
        # residual sum of squares
        RSS <- as.numeric(RVinvR/(n-k))
        # log likelihood
        logLik <- as.numeric( -0.5*n*log(2*pi*RSS) - sum(log(diag(cholV))) - 0.5*RVinvR/RSS )
        logLikR <- as.numeric( logLik + 0.5*k*log(2*pi*RSS) - sum(log(diag(cholWtVinvW))) )
        
        # print current estimates
        if(verbose) print(c(sigma2.k, logLikR, RSS))
        
        # projection matrix
        P <- Vinv - tcrossprod(tcrossprod(VinvW,WtVinvWInv),VinvW)
        
        # vector for later use
        PY <- crossprod(P,Y)
        
        if(reps > 1){
            # Average Information and Scores
            AI <- matrix(NA, nrow=(m+g), ncol=(m+g))
            score <- rep(NA,(m+g))
            for(i in 1:m){
                PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
                score[i] <- -0.5*(sum(P*covMatList[[i]]) - crossprod(Y, PAPY)) # tr(PA) - YPAPY
                AI[i,i] <- 0.5*crossprod(PY, crossprod(covMatList[[i]],PAPY)) # YPAPAPY
                if((i+1) <= m){
                    for(j in (i+1):m){
                        AI[i,j] <- 0.5*crossprod(PY, crossprod(covMatList[[j]],PAPY)) # YPDPAPY
                        AI[j,i] <- AI[i,j]
                    }
                }
                if(g == 1){
                    AI[i,(m+1)] <- 0.5*crossprod(PY, PAPY) # YPIPAPY
                    AI[(m+1),i] <- AI[i,(m+1)]
                }else{
                    for(j in 1:g){
                        AI[i,m+j] <- 0.5*crossprod(PY[group.idx[[j]]], PAPY[group.idx[[j]]]) # YP(I_group)PAPY
                        AI[m+j,i] <- AI[i,m+j]
                    }
                }
            }
            if(g == 1){
                score[m+1] <- -0.5*(sum(diag(P)) - crossprod(PY)) # tr(P) - YPIPY
                AI[(m+1),(m+1)] <- 0.5*crossprod(PY,crossprod(P,PY)) # YPIPIPY
            }else{
                for(i in 1:g){
                    PIPY <- crossprod(P[group.idx[[i]], ],PY[group.idx[[i]]])
                    score[m+i] <- -0.5*(sum(diag(P)[group.idx[[i]]]) - crossprod(PY[group.idx[[i]]])) # tr(P(I_group)) - YP(I_group)PY
                    AI[m+i,m+i] <- 0.5*crossprod(PY[group.idx[[i]]], PIPY[group.idx[[i]]]) # YP(I_group)P(I_group)PY
                    if((i+1) <= g){
                        for(j in (i+1):g){
                            AI[m+i,m+j] <- 0.5*crossprod(PY[group.idx[[j]]], PIPY[group.idx[[j]]]) # YP(I_group2)P(I_group)PY
                            AI[m+j,m+i] <- AI[m+i,m+j]
                        }
                    }
                }
            }
            
            if(dropZeros){
                # remove Zero terms
                AI <- AI[!zeroFLAG,!zeroFLAG]
                score <- score[!zeroFLAG]
            }
            
            # update
            AIinvScore <- solve(AI, score)
            
            if(dropZeros){
                sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + AIinvScore
                sigma2.kplus1[zeroFLAG] <- 0
            }else{
                sigma2.kplus1 <- sigma2.k + AIinvScore
                sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
            }
            
            # step-halving if step too far
            tau <- 1
            while(!all(sigma2.kplus1 >= 0)){
                tau <- 0.5*tau
                if(dropZeros){
                    sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + tau*AIinvScore
                    sigma2.kplus1[zeroFLAG] <- 0
                }else{
                    sigma2.kplus1 <- sigma2.k + tau*AIinvScore
                    sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
                }
            }
            
            # test for convergence
            stat <- sqrt(sum((sigma2.kplus1 - sigma2.k)^2))
            # update estimates
            sigma2.k <- sigma2.kplus1
            if(stat < AIREML.tol){
                converged <- TRUE
                break()
            }
            if(reps == maxIter){
                converged <- FALSE
                warning("Maximum number of iterations reached without convergence!")
                break()
            }
            
        }else{
            # EM step
            sigma2.kplus1 <- rep(NA,(m+g))
            for(i in 1:m){
                PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
                sigma2.kplus1[i] <- (1/n)*(sigma2.k[i]^2*crossprod(Y,PAPY) + n*sigma2.k[i] - sigma2.k[i]^2*sum(P*covMatList[[i]]))
            }
            if(g == 1){
                sigma2.kplus1[m+1] <- (1/n)*(sigma2.k[m+1]^2*crossprod(PY) + n*sigma2.k[m+1] - sigma2.k[m+1]^2*sum(diag(P)))
            }else{
                for(i in 1:g){
                    sigma2.kplus1[m+i] <- (1/n)*(sigma2.k[m+i]^2*crossprod(PY[group.idx[[i]]]) + n*sigma2.k[m+i] - sigma2.k[m+i]^2*sum(diag(P)[group.idx[[i]]]))
                }
            }
            sigma2.k <- sigma2.kplus1
        }
        
    })
    
    # linear predictor
    eta <- fits + crossprod(Vre, VinvR) # X\beta + Zb
    
    return(list(sigma2.k = sigma2.k, AI = AI, converged = converged, zeroFLAG = zeroFLAG, beta=beta, eta=eta, logLikR=logLikR, logLik=logLik, RSS=RSS))
    
}
.subsetCholSigmaInv <-
function(cholSigmaInv, chol.idx) {
    if(length(chol.idx) > 0){
        # subset cholSigmaInv
        SigmaInv <- tcrossprod(cholSigmaInv)
        for(i in sort(chol.idx, decreasing=TRUE)){
            SigmaInv <- SigmaInv[-i,-i] - tcrossprod(SigmaInv[-i,i])/SigmaInv[i,i]
        }
        cholSigmaInv <- t(chol(SigmaInv))
    }
    
    cholSigmaInv
}

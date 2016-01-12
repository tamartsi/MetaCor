estVarCompCI <-
function(x, prop=TRUE){
    if(prop){
        ci <- matrix(NA, nrow=length(x$varComp), ncol=2)
        est <- x$varComp/sum(x$varComp)
        varCompCov <- x$varCompCov
        varCompCov[is.na(varCompCov)] <- 0
        for(i in 1:length(est)){
            deltaH <- rep(-x$varComp[i]/(sum(x$varComp)^2),length(x$varComp))
            deltaH[i] <- deltaH[i] + sum(x$varComp)/(sum(x$varComp)^2)
            varH <- crossprod(deltaH, crossprod(varCompCov, deltaH))
            ci[i,] <- est[i] + sqrt(varH)*qnorm(c(0.025,0.975))
        }
        ci[x$zeroFLAG,] <- NA
        res <- as.data.frame(cbind(est, ci))
        names(res) <- c("Proportion", "Lower 95", "Upper 95")
        
    }else{
        ci <- x$varComp + sqrt(diag(x$varCompCov)) %o% qnorm(c(0.025,0.975))
        res <- as.data.frame(cbind(x$varComp, ci))
        names(res) <- c("Est", "Lower 95", "Upper 95")
    }
    
    print(res)
}

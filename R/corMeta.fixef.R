corMeta.fixef <-
function(beta, var, cov){

		ordered <- check.column.order.input.MetaCor(beta, var, cov)
		if (!ordered$good.input) stop(cat("The column names of beta, var, and cov indicates wrong order of columns or wrong names, got this message: ", "\n", ordered$message, "\n"))

        K <- ncol(beta)

        W.fixed <- 1/var
        Beta.times.W.fixed <- beta*W.fixed
        beta.fixed <- rowSums(Beta.times.W.fixed, na.rm = T)/rowSums(W.fixed, na.rm = T)

        covs <- cbind(var, cov)

        W.covs <- W.fixed^2
        for (i in 1:(K-1)){
                W.covs <- cbind(W.covs, W.fixed[,i]*matrix(W.fixed[,(i+1):K]*2, nrow = nrow(covs)))
        }

        covs.times.W.covs <- covs*W.covs
        var.fixed <- rowSums(covs.times.W.covs, na.rm = T)/rowSums(W.fixed, na.rm = T)^2

        test.stat <- beta.fixed/sqrt(var.fixed)
        pval <- 2*pnorm(abs(test.stat), lower.tail = F)

        return(data.frame(beta = beta.fixed, var = var.fixed, test.stat  = test.stat^2, pval = pval))
}

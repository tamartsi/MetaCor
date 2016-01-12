corMeta.GLS.invert.mats <-
function(beta, var, cov){
	
		ordered <- check.column.order.input.MetaCor(beta, var, cov)
		if (!ordered$good.input) stop(cat("The column names of beta, var, and cov indicates wrong order of columns or wrong names, got this message: ", "\n", ordered$message, "\n"))

	
		if (nrow(beta) == 1) return(corMeta.GLS.1(beta, var, cov))
        S <- ncol(beta)

        cov.mats <- change.cov.mat.representation(var, cov)
        inv.cov.mats <- invert.matrices.jointly(cov.mats)

        inds <- seq(1,  S*S, by = S)

        args1 <- rowSums(inv.cov.mats)

        for (i in 1:length(inds)){
                if (i == 1) betas <- rowSums(beta*inv.cov.mats[,inds[i]:(inds[i] + S-1)]) else{
                betas <- betas + rowSums(beta*inv.cov.mats[,inds[i]:(inds[i] + S-1)])}
        }
        betas <- betas/args1
        vars <- 1/args1

        test.stats <- betas^2/vars
        pvals <- pchisq(test.stats, df = 1, lower.tail = F)

        return(data.frame(beta = betas, var = vars, test.stat = test.stats, pval = pvals))
}

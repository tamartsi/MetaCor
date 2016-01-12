test2levelInteraction <-
function(beta, var, cov){
	K <- ncol(beta)
	if (K != 2) stop(paste("The are", K, "strata, there should be 2"))
	  	
	beta.diff <- beta[,1] - beta[,2]
	
	var <- var[,1] + var[,2] - 2*cov
	
	test.stat <- beta.diff/sqrt(var)
	pval <- 2*pnorm(abs(test.stat), lower.tail = F)

	res <- data.frame(beta = beta.diff, var = var,  test.stat = test.stat^2, pval =pval)
	colnames(res) <- c("beta", "var", "test.stat", "pval")
	return(res)
	
}

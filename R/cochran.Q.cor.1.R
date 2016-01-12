cochran.Q.cor.1 <-
function(beta.1, var.1, cov.1, testType){
	stopifnot(is.element(testType, c("MetaGLS", "MetaCor.fixef")))
	
	S <- length(beta.1)
	beta.1 <- drop(beta.1)
	var.1 <- drop(var.1)
	cov.1 <- drop(cov.1)
	
	na.strat <- which(is.na(beta.1))
	if (length(na.strat) > 0){
		if (length(na.strat) >= S-1) return(NA)
		for (i in 1:length(na.strat)){
			mis.strat.name <- 	names(na.strat)[i]
			mis.strat.name <- strsplit(mis.strat.name, split = ".", fixed = T)[[1]][2]
			beta.1 <- beta.1[-grep(mis.strat.name, names(beta.1))]
			var.1 <- var.1[-grep(mis.strat.name, names(var.1))]
			cov.1 <- cov.1[-grep(mis.strat.name, names(cov.1))]
		}
		
	}
	
	S <- length(beta.1)
	weights <- 1/var.1
	
	W <- diag(1/var.1)
	
	cov.beta <- diag(drop(var.1))
	ind.first <- 1
	for (i in 1:(S-1)){
		cov.beta[i,(i+1):S] <- cov.beta[(i+1):S, i] <- cov.1[ind.first:(ind.first+S-i-1)]
		ind.first <- ind.first+S-i
	}
	
	inv.cov.beta <- solve(cov.beta)
	W.sumInvCov <- diag(weights/sum(inv.cov.beta))
	
	one.mat <- matrix(1, nrow = S, ncol = S)
	
	if (testType == "MetaGLS"){
		A <- W - 2* W.sumInvCov %*% one.mat %*% inv.cov.beta + inv.cov.beta %*% one.mat %*% inv.cov.beta * sum(weights)/(sum(inv.cov.beta)^2)	
		beta.meta <- sum(inv.cov.beta %*% beta.1)/sum(inv.cov.beta)
	} else if (testType == "MetaCor.fixef"){
		A <- W - 1/sum(weights)*weights %*% t(weights)
		beta.meta <- sum(weights*beta.1)/sum(weights)
	}
	
	K <- crossprod(A,cov.beta)
	Q.stat <- sum((beta.1 - beta.meta)^2/var.1)
	res <- davies(Q.stat, eigen(K, only.values = T)$values)$Qq
	return(res)
}

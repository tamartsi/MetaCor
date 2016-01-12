twoStratCovs <-
function(covMatList, varComps, IDsList){
	stopifnot(length(varComps) == length(IDsList))
	n.strat <- length(varComps)
	analysis.names <- names(varComps)
	
	covMatList.names <- names(covMatList)
	
	
	data.list <- vector(mode = "list", length = n.strat)
	names(data.list) <- analysis.names
	two.strata.cov <- data.list[1:(n.strat -1)]
	
	for (i in 1:(n.strat-1)){
	two.strata.cov[[i]] <- data.list[(i+1):n.strat]
	for (j in 1:length(two.strata.cov[[i]])){
		strat.1.ind <- i
		strat.2.ind <- i + j
		
		## build all component correlation matrices
		
		cov.12 <- matrix(0, nrow = length(IDsList[[strat.1.ind]]), ncol = length(IDsList[[strat.2.ind]]))
		
		for (k in 1:length(covMatList)){
			ind <- grep(covMatList.names[k], names(covMatList))
			cor.12 <- covMatList[[ind]][ match(IDsList[[strat.1.ind]], rownames(covMatList[[ind]])), match(IDsList[[strat.2.ind]], colnames(covMatList[[ind]]))]
			
			ind.vc.strat.1 <- grep(covMatList.names[k], names(varComps[[strat.1.ind]])) 
			ind.vc.strat.2 <- grep(covMatList.names[k], names(varComps[[strat.2.ind]])) 
			
			cov.12.k <- cor.12*sqrt(varComps[[strat.1.ind]][ind.vc.strat.1])*sqrt(varComps[[strat.2.ind]][ind.vc.strat.2])
			
			cov.12 <- cov.12 + cov.12.k
		}
				
		two.strata.cov[[i]][[j]] <- cov.12
		
	}
}
return(two.strata.cov)	
}

stratLMMTest <-
function(Y, W, G, covMatList, IDsList, blockSize = 5000, metaCorBlockSize = 50000, testType = "MetaGLS", verbose = TRUE){
	
	if (!is.list(IDsList)) 	stop("participants IDs should be given in a list, each item defines a stratum.")
	if (is.null(names(IDsList))) {
		message("The list of IDs defining the strata should be named, but it has no names. \n New names will be A,B,...")
		names(IDsList) <- LETTERS[1:length(IDsList)]
	}
	if (!is.list(covMatList)) stop("covMatList should be a list")
	if (!is.numeric(G)) stop("genotype matrix G should be numeric, try as.matrix()")
	
	stopifnot(length(Y) == nrow(W))
	if (is.vector(G)) {
		stopifnot(length(Y) == length(G))	
		G <- matrix(G)
	}
	if (is.matrix(G)) stopifnot(length(Y) == nrow(G))
	
	if (is.null(names(Y))) names(Y) <- rownames(W)
	if (is.null(rownames(G))) rownames(G) <- rownames(W)
	
	k <- ncol(W)
	n.strat <- length(IDsList)
	if (n.strat == 1) message("only one strata was defined!")
	
	if (is.null(IDsList)) IDsList <- list(pooled = rownames(W))
	stopifnot(all(is.element(do.call(c, IDsList), rownames(W))))
	
	na.inds.Y <- which(is.na(Y))
	na.inds.W <- which(apply(W, 1, function(x) sum(is.na(x))) > 0)
	na.inds.G <- which(apply(G, 1, function(x) sum(is.na(x))) > 0)
	na.inds <- c(na.inds.Y, na.inds.W, na.inds.G)
	na.inds <- unique(na.inds)
	if (length(na.inds) > 0) {
		message(paste("Some missing outcome, covariates and/or genotypes. Removing", length(na.inds), "individuals with missing data"))
		Y <- Y[-na.inds]
		G <- G[-na.inds, , drop = F]
		W <- W[-na.inds, , drop = F]
		for (i in 1:n.strat){
			IDsList[[i]] <- intersect(IDsList[[i]], rownames(W))
		}

	}
	
	if (verbose) message(paste("There are", length(Y), "individuals in the data"))
	
	if (verbose) message("Defining strata-specific outcome, covariates, genotype, and variance objects...")
	

	strat.names <- names(IDsList)
	
	Y.old <- Y
	W.old <- W
	G.old <- G
	
	data.list <- vector(mode = "list", length = n.strat)
	names(data.list) <- names(IDsList)
	Y <- W <- G <- cholSigmaInv <- var.comps <- data.list
	
	if (verbose) message("Estimating strata-specific variance components...")
	for (i in 1:n.strat){
		Y[[i]] <- Y.old[match(IDsList[[i]], names(Y.old))]
		W[[i]] <- W.old[match(IDsList[[i]], rownames(W.old)), , drop = F]
		G[[i]] <- G.old[match(IDsList[[i]], rownames(G.old)), , drop = F]
	
		temp.varComp <- estVarComp(Y[[i]], W[[i]], covMatList, IDs = IDsList[[i]], verbose = verbose)
		cholSigmaInv[[i]] <- temp.varComp$cholSigmaInv
		var.comps[[i]] <- temp.varComp$varComp[-grep("V_E", names(temp.varComp$varComp))]
	}
	
	if (verbose) message("Calculating covariances between individuals in different strata...")
	if (n.strat > 1) {
			two.strat.covs <- twoStratCovs(covMatList = covMatList, varComps = var.comps, IDsList = IDsList)
			sum.two.strat.covs <- 0
			for (i in 1:length(two.strat.covs)){
				for (j in 1:length(two.strat.covs[[i]])){
					sum.two.strat.covs <- sum.two.strat.covs + sum(abs(two.strat.covs[[i]][[j]]))
				}
			}
			if (sum.two.strat.covs == 0) {
				message("Strata are assumed independent according to covariance between them!")
				strata.indep.flag <- TRUE
			} else strata.indep.flag = FALSE
		}
	
	### define blocks
	
	nloci <- ncol(G.old)
	if (nloci == 1) b.start <- b.end <- 1
	if (nloci > 1) {
		
	}
	b.start <- seq(1, by = blockSize, to = nloci)
	b.end <- seq(min(blockSize, nloci), by = blockSize,  length.out = length(b.start))
	b.end[length(b.end)] <- nloci
	nblocks <- length(b.end)
	
	#### initiallizng the results data.frame
	if (verbose) message("Initializing the results object...")
	
	if (n.strat > 1){
		length.col.cov.names <- n.strat*(n.strat - 1)/2
		col.cov.names <-  rep("", length.col.cov.names)
		if (length.col.cov.names == 1){
			col.cov.names <- paste0("cov.",strat.names[1], ":", strat.names[2])	
		} else{
			ind.1 <- 1
			for (i in 1:(length(strat.names)-1)){
				 col.cov.names[ind.1:(ind.1 + length(strat.names) - i - 1)] <-  paste0("cov.",strat.names[i], ":", strat.names[(i+1):length(strat.names)])	
				 ind.1 <- ind.1 + length(strat.names) - i 
				
				}
			}
	} else col.cov.names <- NULL
		
	col.res.names <- c("snpID",  paste0("Beta.", strat.names), paste0("var.", strat.names), col.cov.names)
	
	
	res <- matrix(NA, ncol = length(col.res.names), nrow = ncol(G.old) , dimnames = list(NULL, col.res.names)) ## change to list(NULL, col.res.names)
	res <- data.frame(res)
	colnames(res) <- col.res.names
	res[,1] <- colnames(G.old)
	
	
	
	
	 Mt <- Xtilde <- Ytilde <- XtX <- beta <- var.beta <- geno <- N.in.strat <- data.list
	if (verbose) message("Calculating individuals strata parameters that are independent of genotypes...")
	
	for (i in 1:n.strat){
		
		C.temp <- cholSigmaInv[[i]]
		CW <- crossprod(C.temp, W[[i]])
		Mt[[i]] <- C.temp - tcrossprod(tcrossprod(C.temp,tcrossprod(chol2inv(chol(crossprod(CW))),CW)),CW)
		Ytilde[[i]] <- crossprod(Mt[[i]],Y[[i]])	
		N.in.strat[[i]] <- length(Ytilde[[i]])
	}
	
	if (verbose) message("Starting association analyses...")
	
	for (b in 1:nblocks){
		if (verbose) message(paste0("Block number ", b, " out of ", nblocks, " blocks"))
		bidx <- b.start[b]:b.end[b]
		
		if (verbose) message(paste0("\t preparing all individual-strata variables...")) 
		
		for (i in 1:n.strat){
			geno[[i]] <- G[[i]][,bidx, drop = F]
		 	freq <- 0.5*colMeans(geno[[i]], na.rm=T)
		    maf <- ifelse(freq < 0.5, freq, 1-freq)
		
			Xtilde[[i]] <- crossprod(Mt[[i]],geno[[i]])
			XtX[[i]] <- colSums(Xtilde[[i]]^2)
			XtX[[i]][which(maf==0)] <- NA
			
		    sY2 <- sum(Ytilde[[i]]^2)
			
			res[bidx,grep(paste0("Beta.", strat.names[i]), colnames(res))] <- as.vector(crossprod(Xtilde[[i]],Ytilde[[i]]))/XtX[[i]]
			res[bidx,grep(paste0("var.", strat.names[i]), colnames(res))] <- ((sY2/XtX[[i]] - (res[bidx,grep(paste0("Beta.", strat.names[i]), colnames(res))])^2)/(N.in.strat[[i]] - k - 1)) 
		}
		
		if (n.strat > 1){
			if (verbose) message("\t calculating covariances between effect estimates of all strata...")
			
			if (strata.indep.flag) 
				res[bidx, grep("cov", colnames(res))] <- 0 else{
						
			two.strata.effect.cov     <- data.list[1:(n.strat -1)]
			
			for (i in 1:(n.strat-1)){
		 			for (j in 1:length(two.strat.covs[[i]])){
						strat.1.ind <- i
						strat.2.ind <- i + j
						two.strata.cov.sub <- two.strat.covs[[i]][[j]]
						mid.1 <- crossprod(Mt[[strat.1.ind]],two.strata.cov.sub)
						mid <- crossprod(t(mid.1), Mt[[strat.2.ind]])
						mid.2 <- crossprod(Xtilde[[i]],mid)/XtX[[i]]
						res[bidx, grep(paste0("cov.", strat.names[strat.1.ind], ":", strat.names[strat.2.ind]), colnames(res))] <- rowSums(mid.2 * t(Xtilde[[i+j]])/XtX[[i+j]])
						
				}
			}}
		
		}
				
		
	}
	
	### add meta analysis and heterogeneity test?
	if (n.strat > 1){
		meta <- MetaCor(res, block.size = metaCorBlockSize, testType = testType, verbose = verbose)
		res$meta.beta <- meta$meta.beta
		res$meta.sd <- meta$meta.sd
		res$meta.test <- meta$meta.test
		res$meta.pval <- meta$meta.pval
		res$hetero.pval <- meta$hetero.pval
		if (n.strat == 2){
			res$hetero.eff <- meta$hetero.eff
			res$hetero.sd <- meta$hetero.sd
			res$hetero.test <- meta$hetero.test
			
		}
		
	} else{ ## just one strata
		res$test <- res[,grep("Beta", colnames(res))]/sqrt(res[,grep("var", colnames(res))])
		res$pval <- pchisq(res$test^2, df = 1, lower.tail = F)
		
		}
	
	return(res)
	
}

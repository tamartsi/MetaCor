MetaCor <-
function(dat, block.size = 50000, testType = "MetaGLS", subgroup.vec = NULL, test.heterogeneity = TRUE, verbose = FALSE){
	
	if (!is.element(testType, c("MetaGLS", "MetaCor.fixef"))){
		message("unfamiliar testType, setting testType to be MetaGLS")
		testType <- "MetaGLS"
	}
	stopifnot(length(grep("snpID", colnames(dat))) == 1)
	
	if (verbose) message("Subsetting the data to matrices of effects, variances, and covariances...")
	
	beta <- as.matrix(dat[,grep("Beta", colnames(dat)), drop = F])
	var <- as.matrix(dat[,grep("var", colnames(dat)), drop = F])
	cov <- as.matrix(dat[,grep("cov", colnames(dat)), drop = F])
	
	if (verbose) message("Checking that all required collumns are in the input and in the right order...")
	
	ordered <- check.column.order.input.MetaCor(beta, var, cov)
	if (!ordered$good.input) stop(cat("The column names of beta, var, and cov indicates wrong order of columns or wrong names, got this message: ", "\n", ordered$message, "\n"))
	
	#### if need to subset by a given vector of groups -- subset the beta, var, and cov matrices accordingly, i.e. find the column numbets corresponding to the given groups and keep only them.  
	if (!is.null(subgroup.vec)){
		if (verbose) message("subsetting the data to the required stata (subgroups)...")
		beta.inds <- sapply(colnames(beta), function(x){is.element(x, paste0("Beta.", subgroup.vec))})
		var.inds  <- sapply(colnames(var), function(x){is.element(x, paste0("var.", subgroup.vec))})
		
		
		cov.first <- sapply(colnames(cov), function(x){
			temp <- strsplit(x, split = ".", fixed = TRUE)[[1]][2]
			temp <- strsplit(temp, split = ":", fixed = TRUE)[[1]][1]
			return(temp)
		})
		
		cov.second <- sapply(colnames(cov), function(x){
			temp <- strsplit(x, split = ".", fixed = TRUE)[[1]][2]
			temp <- strsplit(temp, split = ":", fixed = TRUE)[[1]][2]
			return(temp)
		})

		
		cov.inds.first  <- sapply(cov.first, function(x){is.element(x, subgroup.vec)})
		cov.inds.second  <- sapply(cov.second, function(x){is.element(x, subgroup.vec)})
		cov.inds <- which(cov.inds.first & cov.inds.second)
		
		beta <- beta[, beta.inds, drop = F]
		var <- var[, var.inds, drop = F]
		cov <- cov[, cov.inds, drop = F]
		
	}
	
	if (verbose) message("initializing return values....")
	## initialize return matrices:
	res <- matrix(NA, nrow = nrow(beta), ncol = 4)
	rownames(res) <- rownames(beta)
	colnames(res) <- c("meta.beta", "meta.sd", "meta.test", "meta.pval")
	res <- data.frame(res)
	### if tests of heterogeneity/interaction are needed, columns will be added later  
	
	if (verbose) message("assigning SNPs to blocks...")

	S <- ncol(beta)  ## the number of strata
	
	#### detemine blocks of SNPs (rows of input matrices) 
	nloci <- nrow(beta)
	if (testType == "MetaGLS" & S <= 6){
		if (verbose) message(paste0("Testing using MetaGLS, jointly inverting matrices in blocks of ", block.size, " SNPs..."))
		nblocks <- ceiling(nloci/block.size)
	    # start and end positions for SNP blocks
	    if(nblocks == 1){
	        snp.start <- 1
	        snp.end <- nloci
	    }else{
	        snp.start <- (0:(nblocks-1))*block.size+1
	        snp.end <- c( (1:(nblocks-1))*block.size, nloci )
	    }
	
	} else{
		if (verbose) message("Testing using MetaCor.fixef...")
		snp.start <- 1
		snp.end <- nloci
		}


	
	if (verbose) message("meta-analyzing...")
	
	for (b in 1:length(snp.start)){
		if (verbose) message(paste0("\t block ", b, " out of ", length(snp.start), " blocks..."))
		if (testType == "MetaCor.fixef") {
			res.meta <- corMeta.fixef(beta[snp.start[b]:snp.end[b],, drop = F], var[snp.start[b]:snp.end[b], , drop = F], cov[snp.start[b]:snp.end[b],, drop = F])
			}
		
		if (testType == "MetaGLS" & S <= 6) {			
			res.meta <- corMeta.GLS.invert.mats(beta[snp.start[b]:snp.end[b],, drop = F], var[snp.start[b]:snp.end[b], , drop = T], cov[snp.start[b]:snp.end[b],, drop = F])
			na.inds <- which(is.na(res.meta[,grep("beta", colnames(res.meta))]))
			
			if (length(na.inds) >0){
				if (verbose) message("\t  Testing SNPs with missing effect sizes in at least one strata one by one using MetaGLS...")
				
				for (i.na.inds in 1:length(na.inds)){
					res.meta[na.inds[i.na.inds], ] <- corMeta.GLS.1(beta[(b-1)*block.size + na.inds[i.na.inds],,drop = F], var[(b-1)*block.size + na.inds[i.na.inds],, drop = F], cov[ (b-1)*block.size + na.inds[i.na.inds],, drop = F])
					
				}

			}
			
			}
		
		if (testType == "MetaGLS" & S > 6) {
			if (verbose) message("More than 6 strata, applying MetaGLS on SNPs one by one...")
			res.meta <-  t(apply( cbind(beta, var, cov), 1, function(x){
					corMeta.GLS.1(x[1:S], x[(S+1):(S+S)], x[(2*S+1):(S*(S + 3)/2)])
				}))
				
			}	
		
		res$meta.beta[snp.start[b]:snp.end[b]] <- res.meta[, grep("beta", colnames(res.meta))]
		res$meta.sd[snp.start[b]:snp.end[b]] <- sqrt(res.meta[, grep("var", colnames(res.meta))])
		res$meta.pval[snp.start[b]:snp.end[b]] <- res.meta[, grep("pval", colnames(res.meta))]
		res$meta.test[snp.start[b]:snp.end[b]] <- res.meta[, grep("test.stat", colnames(res.meta))]
		
		}
	
	
	
	###### perform test of heterogeneity (or interaction)
	if (test.heterogeneity){
		if (verbose) message("testing heterogeneity...")
		if (S == 2){
			if (verbose) message("Only two strata, testing for interaction...")
			
			res.hetero <- test2levelInteraction(beta[,,drop = F], var[,,drop = F], cov)
			res$hetero.eff <- res.hetero[, grep("beta", colnames(res.hetero))]
			res$hetero.sd <- sqrt(res.hetero[, grep("var", colnames(res.hetero))])
			res$hetero.test <- res.hetero[, grep("test.stat", colnames(res.hetero))]
			res$hetero.pval <- res.hetero[, grep("pval", colnames(res.hetero))]
			
		}	else { ## S > 2, use Cochran's Q test
			if (verbose) message("More than two strata, testing using Cochran's Q...")
			hetero.res <- t(apply( cbind(beta, var, cov), 1, function(x){
					cochran.Q.cor.1(x[1:S], x[(S+1):(S+S)], x[(2*S+1):(S*(S + 3)/2)], testType = testType)
				}))
			
			hetero.res <- matrix(hetero.res)
			res$hetero.pval <- hetero.res
			
		}	
	} 
	
	res$snpID <- dat[,grep("snpID", colnames(dat))]
	
	return(res)
	
}

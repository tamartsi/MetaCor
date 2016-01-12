corMeta.GLS.1 <-
function(beta.1, var.1, cov.1){
        
        beta.1 <- drop(beta.1)
        var.1 <- drop(var.1)
        cov.1 <- drop(cov.1)
        S <- length(beta.1)

        na.strat <- which(is.na(beta.1))
        if (length(na.strat) > 0){
                if (length(na.strat) == S) return(cbind(beta = NA, var = NA, test.stat = NA, pval = NA))
                if (length(na.strat) == S-1) return(cbind(beta = beta.1[-na.strat], var = var.1[-na.strat], test.stat = beta.1[-na.strat]^2/var.1[-na.strat] , pval = pchisq(beta.1[-na.strat]^2/var.1[-na.strat], df = 1, lower.tail = F) ) )
                for (i in 1:length(na.strat)){
                        mis.strat.name <-       names(na.strat)[i]
                        beta.1 <- beta.1[-grep(mis.strat.name, names(beta.1))]

                        mis.strat.name <- strsplit(mis.strat.name, split = ".", fixed = T)[[1]][2]
                        var.1 <- var.1[-grep(paste0("var.",mis.strat.name), names(var.1))]
                        inds.rm <- c(grep(paste0("cov.", mis.strat.name), names(cov.1)) ,  grep(paste0(":", mis.strat.name), names(cov.1)))
                        cov.1 <- cov.1[-inds.rm]

                }
        }

        S <- length(beta.1)


        cov.beta <- diag(var.1)
        ind.first <- 1
        for (i in 1:(S-1)){
                cov.beta[i,(i+1):S] <- cov.beta[(i+1):S, i] <- cov.1[ind.first:(ind.first+S-i-1)]
                ind.first <- ind.first+S-i
        }

        if (min(eigen(cov.beta, only.values = TRUE)$values) <= 1e-10) return(cbind(beta = NA, var = NA, test.stat = NA, pval = NA))

         arg1 <-  sum(solve(cov.beta, rep(1, nrow(cov.beta))))
         beta.fixed <- sum(solve(cov.beta, beta.1))/arg1
         var.fixed <- 1/arg1


        test.stat <- beta.fixed/sqrt(var.fixed)
        pval <- 2*pnorm(abs(test.stat), lower.tail = F)

		res <- data.frame(beta = beta.fixed, var = var.fixed, test.stat = test.stat^2, pval = pval)
        return(res)
}

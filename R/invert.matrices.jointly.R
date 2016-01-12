invert.matrices.jointly <-
function(cov.mats){
        S <- sqrt(ncol(cov.mats))

        inv.cov.mats <- matrix(0, ncol = ncol(cov.mats), nrow = nrow(cov.mats))

        ind.row <- rep(1:S, each = S)
        ind.col <- rep(1:S, times = S)

        for (i in 1:ncol(cov.mats)){
                rm.col.inds <- union(c((ind.row[i]*S):(ind.row[i]*S - S + 1)), which(ind.col == ind.col[i]))
                inv.cov.mats[,i] <- (-1)^(ind.row[i] + ind.col[i])*det.cov.mat(cov.mats[,-rm.col.inds], S-1)

        }

        det.all <- det.cov.mat(cov.mats, S)

        return(inv.cov.mats/det.all)

}

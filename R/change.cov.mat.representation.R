change.cov.mat.representation <-
function(var, cov){

        stopifnot(nrow(var) == nrow(cov), ncol(var) > 1)
        S = ncol(var)

        cov.mats <- matrix(0, nrow = nrow(var), ncol = S^2)

        row.inds.in.cov.mats <- seq(1, S^2 + 1, by = S)
        first.ind.in.cov <- 1
        last.ind.in.cov <- first.ind.in.cov + (S-2)
        for (i in 1:S){
                if (i == 1) cov.mats[, row.inds.in.cov.mats[i]:(row.inds.in.cov.mats[i+1] -1)] <- cbind(var[,i], cov[,first.ind.in.cov:(last.ind.in.cov)]) else{
                        back.inds <- seq(i, length = i-1, by = S)
                        if (i < S) cov.mats[, row.inds.in.cov.mats[i]:(row.inds.in.cov.mats[i+1] -1)] <- cbind(cov.mats[,back.inds], var[,i], cov[,first.ind.in.cov:(last.ind.in.cov)]) else{
                                cov.mats[, row.inds.in.cov.mats[i]:(row.inds.in.cov.mats[i+1] -1)] <- cbind(cov.mats[,back.inds], var[,i])
                        }}
                first.ind.in.cov <- last.ind.in.cov + 1
                last.ind.in.cov <- first.ind.in.cov + S - i - 2
        }
        return(cov.mats)

}

det.cov.mat <-
function(cov.mats, S){

        stopifnot(ncol(cov.mats) == S^2)

        if (S == 1) return(cov.mats)
        if (S == 2) return(cov.mats[,1]*cov.mats[,4] - cov.mats[,2]*cov.mats[,3])
        if (S > 2){
                dets <- rep(0, nrow(cov.mats))
                sign.1 <- 1
                for (i in 1:S){
                        dets <- dets + sign.1*cov.mats[,i]*det.cov.mat(cov.mats = cov.mats[,-c(union(1:S, seq(i, ncol(cov.mats), by = S)) )], S = S-1)
                        sign.1 <- -sign.1
                }
                return(dets)
        }

}

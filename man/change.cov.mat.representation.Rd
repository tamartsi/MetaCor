\name{change.cov.mat.representation}
\alias{change.cov.mat.representation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Takes variances and covarianves of effect estimates and puts them in matrix that is used to calculate MetaGLS estimators. 
}
\description{
This function is used by corMeta.GLS.inverst.mats, the function that calculates the MetaGLS estimators for many SNPs at the same time. 
It puts variances and covariances of effect estimates from stratified analysis in a single matrix, in which each rows corresponds to  
a stacked covariance matrix of the estimated effects. 
}
\usage{
change.cov.mat.representation(var, cov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var}{
   a p by S matrix of variances of the estimated effects of each of p SNPs on a trait in each of S strata. 
}
  \item{cov}{
a p by (S*(S-1)/2) matrix representing p vectors of covariances between the estimated effects of each of p SNPs on a trait in each of S strata.
}
}
\details{
The given variances and covrariances are put in a single matrix, in which each row corresponds to the covariance matrix of the SNP effects of a single SNP across multiple strata. This function is used by corMeta.GLS.inverst.mats. 
}
\value{
An p  by s^2 matrix in which each row represent a covariance matrix, as the concatanated rows of the matrix. 
}

\author{
Tamar Sofer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
### create matrices of covariances and variances

n.strat <- 5
n.est <- 5
## create a positive definite matrix to be a covariance matrix
sigma <- matrix(rnorm(n.strat^2), n.strat, n.strat)
sigma <- t(sigma) %*% sigma

vars <- matrix(rep(diag(sigma), n.est), nrow = n.est, ncol = n.strat, byrow = T)
colnames(vars) <- paste0("var.s", 1:n.strat)
covs <- matrix(rep(sigma[lower.tri(sigma)], n.est), nrow = n.est, byrow = T)
colnames(covs) <- c(paste0("cov.s1:", paste0("s", 2:n.strat)), paste0("cov.s2:", paste0("s", 3:n.strat)), paste0("cov.s3:", paste0("s", 4:n.strat)) , "cov.s4:s5")

### now change their representation:
change.cov.mat.representation(vars, covs)

}

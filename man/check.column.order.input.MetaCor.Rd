\name{check.column.order.input.MetaCor}
\alias{check.column.order.input.MetaCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check that the input is organized in the right way for MetaCor functions. 
}
\description{
The functions that compute the meta-analysis estimators from MetaCor family (corMeta.fixef, corMeta.GLS.1 and corMeta.GLS.invert.mats) assume a specific order on the input data, since they perform matrix operations. This function checks that the order is correct. 
}
\usage{
check.column.order.input.MetaCor(beta, var, cov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
  a p by  s matrix of the estimated effects of each of p SNPs on a trait in each of s strata. It must have column names of the form "Beta.stratum_name"
}
  \item{var}{
  a p by  s matrix of variances of the estimated effects of each of p SNPs on a trait in each of s strata. It must have column names of the form "var.stratum_name"
}
  \item{cov}{
  a p by  (s*(s-1)/2) matrix of a vector of length p of covairances between the estimated effects of each of p SNPs on a trait in each of s strata.It must have column names of the form "cov.stratum_name_1:stratum_name_2"
}
}
\details{
The order of the columns is automatically correct when this function is called by stratLMMTest. This function is important if the user is performing additional meta-analysis based on effect estimates, they variances and covariances. 
}
\value{
A list with a logical element named good.input (if TRUE, the columns are ordered correctly according to their names, otherwise it is FALSE). 
If good.input is FALSE, additional message is provided with the type of problem. 
}

\author{
Tamar Sofer
}

\examples{
### create a small data set with the right formating:
require(mvtnorm)
n.strat <- 5
n.est <- 5
## create a positive definite matrix to be a covariance matrix
sigma <- matrix(rnorm(n.strat^2), n.strat, n.strat)
sigma <- t(sigma) %*% sigma

effects <- rmvnorm(n.est, c(2,2,3,3,0), sigma)
colnames(effects) <- paste0("Beta.s", 1:n.strat)
vars <- matrix(rep(diag(sigma), n.est), nrow = n.est, ncol = n.strat, byrow = T)
colnames(vars) <- paste0("var.s", 1:n.strat)
covs <- matrix(rep(sigma[upper.tri(sigma)], n.est), nrow = n.est, byrow = T)
colnames(covs) <- c(paste0("cov.s1:", paste0("s", 2:n.strat)), paste0("cov.s2:", paste0("s", 3:n.strat)), paste0("cov.s3:", paste0("s", 4:n.strat)) , "cov.s4:s5")

## order is correct:
check.column.order.input.MetaCor(effects, vars, covs)

## we want to meta-analyzed only strata s1 and s2. Are we choosing them correctly?
check.column.order.input.MetaCor(effects[,1:2], vars[,1:2], covs[,1, drop = F])
## yes!

check.column.order.input.MetaCor(effects[,1:2], vars[,1:2], covs[,2, drop = F])
## no!

}

\name{corMeta.fixef}
\alias{corMeta.fixef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The fixed effects meta-analysis estimator accounting for correlations between effect estimates
}
\description{
Computes the fixed effects meta-analysis estimator while accounting for correlations 
between effect estimates in the variance estimation. 
The computation is done for many SNPs at once.
}
\usage{
corMeta.fixef(beta, var, cov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
 a p by S matrix of estimated effects of each of p SNPs on a trait in each of S strata.
}
  \item{var}{
 a p by S matrix of variances of the estimated effects of each of p SNPs on a trait in each of S strata.
}
  \item{cov}{
  a p by (S*(S-1)/2) matrix representing p vectors of covariances between the estimated effects of each of p SNPs on a trait in each of S strata.
}
}
\details{
This function requires column names for each of the columns. It is used by 
stratLMMTest, which prepares the data with correct column names. 
}
\value{
a matrix with beta, var, test.stat and pval corresponding to the meta-analysis estimates for each of the SNPs, 
their variances, the test statistics and their p-values. 
}

\author{
Tamar Sofer
}
\note{
This function is called by MetaCor and it is advisable to use MetaCor.
It can also be called directly, but then the test of heterogeneity is not performed. 
Also, MetaCor allows for selection of strata while this function does not. 
}

%%
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
covs <- matrix(rep(sigma[lower.tri(sigma)], n.est), nrow = n.est, byrow = T)
colnames(covs) <- c(paste0("cov.s1:", paste0("s", 2:n.strat)), paste0("cov.s2:", paste0("s", 3:n.strat)), paste0("cov.s3:", paste0("s", 4:n.strat)) , "cov.s4:s5")

## test:
corMeta.fixef(effects, vars, covs) 
}

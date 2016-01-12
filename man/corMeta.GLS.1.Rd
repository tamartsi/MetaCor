\name{corMeta.GLS.1}
\alias{corMeta.GLS.1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The efficient meta-analysis estimator accounting for correlations between effect estimates
}
\description{
Computes the fixed effects meta-analysis estimator while accounting for correlations 
between effect estimates in the variance estimation. 
The computation is done for a single SNP at a time.
}
\usage{
corMeta.GLS.1(beta.1, var.1, cov.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta.1}{
A vector of S effect estimates
}
  \item{var.1}{
A vector of S variances of the effect estimates in beta.1
}
  \item{cov.1}{
A vector of S*(S-1)/2 covariances between the effect estimates in beta.1
}
}
\details{
Computes the MetaGLS (MetaCor1 form the MetaCor paper) estimator and test for a single SNP.
}
\value{
A data frame with 1 row, with beta, var, test.stat and pval providing the meta-analysis effect estimate, 
its variance, test statistic, and p-value.
}

\author{
Tamar Sofer
}
\note{
This function is called by MetaCor and it is advisable to use MetaCor.
It can also be called directly, but then the test of heterogeneity is not performed. 
Also, MetaCor allows for selection of strata while this function does not. 
}


\examples{
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


corMeta.GLS.1(effects[1,], vars[1,], covs[1,], testType = "MetaGLS")
}


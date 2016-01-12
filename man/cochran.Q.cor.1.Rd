\name{cochran.Q.cor.1}
\alias{cochran.Q.cor.1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The modified Cochran's Q test for heterogeneity that accounts for correlations between effect estiates
}
\description{
Calculates the modified Cochran's Q test for heterogeneity between effect estimates, while accounting 
for correlations between the effect estimates.
}
\usage{
cochran.Q.cor.1(beta.1, var.1, cov.1, testType)
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
  \item{testType}{
  Which test should the Cochran Q test be based on? MetaCor.fixef (MetaCor1), or MetaGLS (MetaCor2)? preferably, one should 
  use the same testType that was used to compute the meta-analysis estimator. 
}
}

\value{
A p-value of the test of heterogeneity. 
}

\author{
Tamar Sofer
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


cochran.Q.cor.1(effects[1,], vars[1,], covs[1,], testType = "MetaGLS")
cochran.Q.cor.1(effects[1,], vars[1,], covs[1,], testType = "MetaCor.fixef")


}


\name{test2levelInteraction}
\alias{test2levelInteraction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test SNPxstrata interaction when there are only two strata
}
\description{
Estimate and test the difference in SNP effect the trait between two strata. 
}
\usage{
test2levelInteraction(beta, var, cov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
  a p by s matrix of estimated effects of each of p SNPs on the trait in each of 2 strata. 
}
  \item{var}{
 a p by 2 matrix of variances of the estimated effects of each of p SNPs on the trait in each of 2 strata. 
}
  \item{cov}{
a p by 1 matrix of a vector of length p of covairances between the estimated effects of each of p SNPs on the trait in each of 2 strata. 
}
}
\details{
  For each SNP, estimates the difference in SNP effect between strata, and tests the null hypothesis that this effect is zero. 
}
\value{
A matrix with column corresponding to beta = the estimated effect difference (equivalently, the interaction effect), var = the variance
of the estimated beta, test.stat = the test.statistics and pval = the p-value of the test. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tamar Sofer
}
\note{
This function is used  by stratLMMTest when two strata are analyzed. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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

## this will give an error:
# test2levelInteraction(effects, vars, covs)
### error! there are more than 2 strata

test2levelInteraction(effects[,1:2], vars[,1:2], covs[,1])

test2levelInteraction(effects[,4:5], vars[,4:5], covs[,grep("cov.s4:s5", colnames(covs))])

}

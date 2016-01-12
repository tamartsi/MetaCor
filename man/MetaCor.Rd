\name{MetaCor}
\alias{MetaCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine stratum-specific effect estimates in meta-analysis using MetaCor and tests for association and for heterogeneity. 
}
\description{
Given association analysis results and covariances between the stratum-specific estimates, perform meta-analysis using MetaCor. 
}
\usage{
MetaCor(dat, block.size = 50000, testType = "MetaGLS", subgroup.vec = NULL, test.heterogeneity = TRUE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A matrix with rows corresponding to m SNPs, and columns with the following information 
(note that this matrix is prepared by the function stratLMMTest that calls this function)
ID for each SNP, columns with names "Beta.strat_name", "var.strat_name", "cov.strat_name:other_strat_name" corresponding to 
estimated SNP effect, variances of these estimators, and covariances between them, for a few strata (at least 2).
}
  \item{block.size}{
Number of SNPs to meta-analyze together (to speed computation). Used only if MetaGLS is testType. Default is 50,000.
}
  \item{testType}{
Either MetaGLS (MetaCor2 from paper) or MetaCor.fixef (MetaCor1 from paper). If an unfamiliar test name is specified, will be taken to be MetaGLS.
}
  \item{subgroup.vec}{
If one wants to meta-analyze only some of the strata, specify the names in this vector of strings. 
}
  \item{test.heterogeneity}{
Should a test of heterogeneity between stratum-specific effects be performed? default is yes. 
}
  \item{verbose}{
If verbose, print progress updates. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
return a matrix with information for each SNP - 
meta.beta, meta.se, meta.test and meta.pval are the meta-analyzed effect sizes, standard errors, test statistics and p-values. Similar information is given for the test of heterogeneity using hetero.pval, etc.  
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tamar Sofer
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
 ### This function is used by stratLMMTest when combining the effect estimates from all strata. Example of usage outside of stratLMMTest:
  
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
dat <- cbind(data.frame(snpID = 1:n.est), effects, vars, covs)

## test:
MetaCor(dat) ## uses MetaGLS, MetaCor2 from the paper
MetaCor(dat, testType = "MetaCor.fixef") ## MetaCor1 from paper
MetaCor(dat, testType = "MetaCor.fixed") ## if wrong name of test, it uses MetaGLS

## can choose to meta-analyze only part of the strata:
MetaCor(dat, subgroup.vec = c("s1", "s2", "s5"), testType = "MetaCor.fixef")

}


\name{invert.matrices.jointly}
\alias{invert.matrices.jointly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Joint inversion of many covariance matrices
}
\description{
Invert p covariance matrices of dimension S by S,  reprented in a single matrix of dimension p by (s^2) together. 
This is used to speed up MetaGLS computations, that are based on estimated effects across multiple strata, and their covariances. 
}
\usage{
invert.matrices.jointly(cov.mats)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov.mats}{
 A p by s^2 matrix representing p s by s covariance matrices between SNP effects in s strata. 
}
}
\details{
This function is called by corMeta.GLS.invert.mats, which is called by MetaCor when testType is MetaGLS.
}
\value{
 A p by s^2 matrix representing p s by s inversees of the input covariance matrices.
}

\author{
Tamar Sofer
}
\note{
Called by corMeta.GLS.invert.mats
}


\examples{

n.strat <- 5
n.est <- 5
## create a positive definite matrix to be a covariance matrix
sigma.1 <- matrix(rnorm(n.strat^2), n.strat, n.strat)
sigma.1 <- t(sigma.1) %*% sigma.1

vars <- matrix(rep(diag(sigma.1), n.est), nrow = n.est, ncol = n.strat, byrow = T)
colnames(vars) <- paste0("var.s", 1:n.strat)
covs <- matrix(rep(sigma.1[lower.tri(sigma.1)], n.est), nrow = n.est, byrow = T)
colnames(covs) <- c(paste0("cov.s1:", paste0("s", 2:n.strat)), paste0("cov.s2:", paste0("s", 3:n.strat)), paste0("cov.s3:", paste0("s", 4:n.strat)) , "cov.s4:s5")


### now change their representation:
full.covs.1 <- change.cov.mat.representation(vars, covs)



sigma.2 <- matrix(rnorm(n.strat^2), n.strat, n.strat)
sigma.2 <- t(sigma.2) %*% sigma.2

vars <- matrix(rep(diag(sigma.2), n.est), nrow = n.est, ncol = n.strat, byrow = T)
colnames(vars) <- paste0("var.s", 1:n.strat)
covs <- matrix(rep(sigma.2[lower.tri(sigma.2)], n.est), nrow = n.est, byrow = T)
colnames(covs) <- c(paste0("cov.s1:", paste0("s", 2:n.strat)), paste0("cov.s2:", paste0("s", 3:n.strat)), paste0("cov.s3:", paste0("s", 4:n.strat)) , "cov.s4:s5")

### now change their representation:
full.covs.2 <- change.cov.mat.representation(vars, covs)

full.covs <- rbind(full.covs.1, full.covs.2)
full.inv.covs <- invert.matrices.jointly(full.covs)
all(abs(full.inv.covs[10,] - solve(sigma.2)[1:25]) < 1e-10)
}

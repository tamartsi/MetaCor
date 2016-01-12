\name{estVarComp}
\alias{estVarComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate variance components
}
\description{
Estimate the variance components in a mixed model, via the AI-REML procedure. 
}
\usage{
estVarComp(Y, W, covMatList, IDs = NULL, start = NULL, group.var = NULL, AIREML.tol = 1e-06, maxIter = 100, dropZeros = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
An n vector of quantitative outcomes. 
}
  \item{W}{
An n by k matrix of k covariates for each of n participants. 
}
  \item{covMatList}{
A list of covariance matrices modeling the correlations between study participants (e.g. kinship/GRM matrix, household matrix)/
}
  \item{IDs}{
A vectpr of IDs that should be in the model. Should be a subset of the rownames of W. Will be taken to be W if not provided. 
}
  \item{start}{
A potential vector of starting values for the variance components. 
}
  \item{AIREML.tol}{
Tolerance for convergence of the AI-REML procedure. 
}
  \item{maxIter}{
Maximum number of iterations of the AI-REML procedure.
}
  \item{dropZeros}{
If TRUE covariance matrices with respective variance components being equal to zero will be dropped and variance
components will be re-computed to the other covariance matrices. 
}
  \item{verbose}{
If TRUE progress will be reported. 
}
}
\details{
More soon...
}
\value{
A list with estimated variance componets, and the cholesky decomposition of the inverse covariance matrix, and more. 
}

\author{
Matt Conomos
}
\note{
More details soon!
}



\name{det.cov.mat}
\alias{det.cov.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the determinants of multiple covairance matrices jointly
}
\description{
Calculates the determinants of multiple covairance matrices jointly using recurssion. 
This is used by the function invert.matrices.jointly, and ultimately it is used to speed-up the 
calculation of the MetaGLS test statistics for many SNPs together.
}
\usage{
det.cov.mat(cov.mats, S)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov.mats}{
  A p by S^2 matrix of covairance matrices between effect estimates across S strata, corresponding to p SNPs 
}
  \item{S}{
 The number of strata
}
}

\value{
A vector of determinants of the p covariance matrices.
}

\author{
Tamar Sofer
}
\note{
This function is called by invert.matrices.jointly, which is called by corMeta.GLS.invert.mats.
}


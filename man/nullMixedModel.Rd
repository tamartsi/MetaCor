\name{nullMixedModel}
\alias{nullMixedModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fits the null model (model without any genotypes) based on estimated variance components. 
}
\description{
This function is called by stratLMMTest. It fits an linear mixed model 
to the trait Y and covariates W based on computed variance components. 
}
\usage{
nullMixedModel(Y, W, VC, IDs = NULL, group.var = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
An n vector of quantitative outcomes. 
}
  \item{W}{
An n by k matrix of k covariates for each of n participants. 
}
  \item{VC}{
Variance components, estimate by estVarComp.
}
  \item{IDs}{
A vectpr of IDs that should be in the model. Should be a subset of the rownames of W. Will be taken to be W if not provided. 
}
  \item{verbose}{
If TRUE progress will be reported. 
}
}
\details{
Called by stratLMMTest.
}
\value{
A null mixed model list. More details soon. 
}
\author{
Matt Conomos
}
\note{
More details - soon!
}

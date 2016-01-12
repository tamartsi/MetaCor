\name{twoStratCovs}
\alias{twoStratCovs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates covariance matrices between any two strata 
}
\description{
Given a list of participants IDs and variance components in each strata, calculates a list containing the covariance matrices between any two strata defined by a list of participant IDs
}
\usage{
twoStratCovs(covMatList, varComps, IDsList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covMatList}{
  A list of pre-defined covariance matrices between all individuals in the data, based on a random effect model. 
}
  \item{varComps}{
  The variance components estimated based on the list of covariance matrices covMatList, in each of the strata defined by IDsList.
}
  \item{IDsList}{
 A list of IDs defining strata (each item in the list provides IDs for people in a single strata). Overlapping
 individuals between strata are allowed. 
}
}
\details{
This function is used by stratLMMTest function. 
}
\value{
A list of lists of between-strata covariance matrices. The i-th item in the list holds a list 
of covariance matrices between the individuals in the i-th strata and the individuals in the i+1,...,S strata (out of a total of S strata).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Tamar Sofer
}

\examples{
### this function is normally used by stratLMMTest

## upload simulated genotype data data
data(G)
## load a data frame with simulated covariates info
data(W)
## load the simulated outcome
data(Y)
## load the simple kinship matrix for the simulated data and create a list of covariance matrices (containing here only a single, kinship matrix)
data(K)
covMatList = list(kinship = K)

## create a list of IDs to define the strata
all(W$ID == G$ID) #TRUE
IDsList <- list(male = W$ID[which(W$sex == 1)], female = W$ID[which(W$sex == 0)])  


## put information in a list by strata:
n.strat <- length(IDsList)

strat.names <- names(IDsList)
	
	Y.old <- Y
	W.old <- W
	G.old <- G
	
	data.list <- vector(mode = "list", length = n.strat)
	names(data.list) <- names(IDsList)
	Y <- W <- G <- cholSigmaInv <- var.comps <- data.list
  ## Estimating strata-specific variance components...
  for (i in 1:n.strat){
		Y[[i]] <- Y.old[match(IDsList[[i]], names(Y.old))]
		W[[i]] <- W.old[match(IDsList[[i]], rownames(W.old)), , drop = F]
		G[[i]] <- G.old[match(IDsList[[i]], rownames(G.old)), , drop = F]
	
		temp.varComp <- estVarComp(Y[[i]], W[[i]], covMatList, IDs = IDsList[[i]], verbose = verbose)
		cholSigmaInv[[i]] <- temp.varComp$cholSigmaInv
		var.comps[[i]] <- temp.varComp$varComp[-grep("V_E", names(temp.varComp$varComp))]
	}
two.strat.covs <- twoStratCovs(covMatList = covMatList, varComps = var.comps, IDsList = IDsList)

}


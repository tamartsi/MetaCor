\name{MetaCor-package}
\alias{MetaCor-package}
\docType{package}
\title{
An R package to perform genetic association study with stratification, 
when individuals are correlated across the strata. 
}
\description{ Analyzed stratified a stratigied genetic association study with correlated individuals between strata

\details{

Code used to perform stratified SNP-trait association study when individuals are correlated (and possibly overlapping) between the strata. 
Effect estimates from the individual strata are meta-analyzed using MetaCor. 
}
\author{ Tamar Sofer
\packageAuthor{MetaCor}

Maintainer: Tamar Sofer \packageMaintainer{MetaCor}
}

\examples{

#### load the data - note: all data in the package is simulated and is very simple. 
### load a data frame with genotype data
data(G)
G[1:3,]
## note: sex_stratInd are indicators for selection of a maximum set of participants if we want the two sex strata independent.

## load a data frame with covariates info
data(W)
head(W)

## load the outcome
data(Y)
head(Y)

## load the simplified kinship matrix
data(K)
K[1:5,1:5] ## uncorrelated individuals since they are from different families
K[which(W$ID_fam == 1),which(W$ID_fam == 1)] ## these individuals are from the same family and have kinship coefficient 0.5.

## Prepare to fitting a stratified model with strata of sex - create a list of participants IDs in each of the strata:
## check:
all(W$ID == G$ID) #TRUE
IDs.list.gender <- list(male = W$ID[which(W$sex == 1)], female = W$ID[which(W$sex == 0)])
IDs.list.gender.stratInd <- list(male = W$ID[which(W$sex == 1 & W$sex_stratInd == 1)], female = W$ID[which(W$sex == 0 & W$sex_stratInd == 1)])

### fit a stratified model for males and females separately, and combine the results using MetaCor using a single function that does it all:
strat.model <- stratLMMTest(Y, as.matrix(W[,c("intercept", "bmi", "bmi2", "age")]), as.matrix(G[,grep("SNP", colnames(G))]), covMatList = list(kinship = K), IDsList = IDs.list.gender, verbose = T, testType = "MetaGLS")

### When according to the kinship matrix, the two strata are indpendent... (and also trying to use MetaCor1, names here "MetaCor.fixef") 
strat.model <- stratLMMTest(Y, as.matrix(W[,c("intercept", "bmi", "bmi2", "age")]), as.matrix(G[,grep("SNP", colnames(G))]), covMatList = list(kinship = K), IDsList = IDs.list.gender.stratInd, verbose = T, testType = "MetaCor.fixef")

## if the list of IDs provides only a single stratum, a pooled model is fit:
pooled.model   <- stratLMMTest(Y,  as.matrix(W[,c("intercept", "bmi", "bmi2", "age", "sex")]), as.matrix(G[,grep("SNP", colnames(G))]), covMatList = list(kinship = K), IDsList = list(pooled = rownames(W)), verbose = T, testType = "MetaCor.fixef")





}

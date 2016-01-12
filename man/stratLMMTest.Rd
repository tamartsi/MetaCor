\name{stratLMMTest}
\alias{stratLMMTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test genetic variants using a stratified linear mixed model 
}
\description{
Takes a study data - outcome, covariants, genetic data, covariance matrices, and a list defining strata - and perfomrs a stratified 
trait-genotypes association analysis (and test).
}
\usage{
stratLMMTest(Y, W, G, covMatList, IDsList, blockSize = 5000, metaCorBlockSize = 50000, testType = "MetaGLS", verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  An n vector of outcomes/trait values for n perticipants 
}
  \item{W}{
 An n by k matrix of covariates values corresponding to n participants. Row names should be the participants IDs. 
}
  \item{G}{
An n by p matrix of p genotype values for each of n paticipants. 
}
  \item{covMatList}{
 A named list of matrices that model the correaltion structure between the participants (e.g. kiship/GRM matrix).
}
  \item{IDsList}{
 A named list of participants IDd in strata. IDs can overlap between strata. 
}
  \item{blockSize}{
Association analysis of blockSize SNPs is computed together to speed up comoputation.
}
  \item{metaCorBlockSize}{
MetaGLS (MetaCor2 in the paper) test statistics are calculated jointly for metaCorBlockSize SNPs to speed up computation. 
}
  \item{testType}{
 MetaCor test type - MetaCor.fixef (MetaCor1) or MetaGLS (MetaCor2).
}
  \item{verbose}{
If verbose is TRUE, progress is reported via messages. 
}
}

\value{
A matrix providing for each SNP the estimated effect size, SEs, and p-value in each stratum, as well as the meta-analyzed
estimated effect sizes, SEs and p-values, and the results from a test of heterogeneity of effects between strata. 
}

\author{
Tamar Sofer
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



\name{K}
\alias{K}
\docType{data}
\title{
A simple kinship matrix for the simulated genotype matrix.
}
\description{
A kinship matrix to help demonstrate the package usage. Provides the probability that each SNPs is identical by simulation between any two individuals in the simulated data. 
}
\usage{data("K")}

\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(K)
data(W)
K[1:5,1:5] ## uncorrelated individuals since they are from different simulated families
K[which(W$ID_fam == 1),which(W$ID_fam == 1)] ## these individuals are from the same simulate family and have kinship coefficient 0.5.
}
\keyword{datasets}

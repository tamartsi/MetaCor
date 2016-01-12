\name{W}
\alias{W}
\docType{data}
\title{
A simulated design matrix
}
\description{
Simulated design matrix and a few other variables to help demonstrate the package usage.
}
\usage{data("W")}
\format{
  A simulated data frame with covariate information for 6000 individuals. The variables are
  \describe{
    \item{\code{intercept}}{a numeric vector}
    \item{\code{bmi}}{a numeric vector}
    \item{\code{bmi2}}{a numeric vector}
    \item{\code{age}}{a numeric vector}
    \item{\code{ID}}{a numeric vector}
    \item{\code{sex}}{a numeric vector}
    \item{\code{ID_fam}}{Family ID (individuals are assumed correlated within families)}
    \item{\code{sex_stratInd}}{An indicator vector of 1 and 0s. When the indicator is 0, an individual should be removed from a sex-stratified analysis if one wants the strata to be independnet of each other. }
  }
}
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
data(W)
head(W)
}
\keyword{datasets}

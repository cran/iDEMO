\name{idemo.MTTF}
\alias{idemo.MTTF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean of the lifetime distribution
}
\description{
The function idemo.MTTF computes the mean of the lifetime distribution.
}
\usage{
idemo.MTTF(Eta, sEta, W)
}
\arguments{
  \item{Eta}{
the drift rate of the degradation model.
}
  \item{sEta}{
the standard deviation of the random effect.
}
  \item{W}{
the threshold for the degradation path.
}

}

\seealso{
\code{\link{idemo.PDF}},\code{\link{idemo.CDF}}
}

\examples{
library(gsl)
Eta <- 0.00204
sEta <- 0.00042
W <- 10
idemo.MTTF(Eta,sEta,W)
}

\references{

Cheng, Y. S. and Peng, C. Y. (2012), "Integrated Degradation Models in R Using iDEMO." \emph{Journal of Statistical Software}, \bold{49(2)}, 1-22.

Peng, C. Y. and Tseng, S. T. (2009), "Mis-specification Analysis of Linear Degradation Models." \emph{IEEE Transactions on Reliability}, \bold{58}, 444-455.

}

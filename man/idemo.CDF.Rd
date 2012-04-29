\name{idemo.CDF}
\alias{idemo.CDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cumulative density function of the lifetime distribution.
}
\description{
The function idemo.CDF computes the cumulative density function of the lifetime distribution.
}
\usage{
idemo.CDF(t, Eta, sEta, sB, W)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
time.
}
  \item{Eta}{
the drift rate of the degradation model.
}
  \item{sEta}{
the standard deviation of the random effect.
}
  \item{sB}{
the diffusion coefficient of the Brownian Motion.
}
  \item{W}{
the threshold for the degradation path.
}

}

\seealso{
\code{\link{idemo.PDF}},\code{\link{idemo.MTTF}}
}
\examples{
x <- 5000
Eta <- 0.00204
sEta <- 0.00042
sB <- 0.0104
W <- 10
idemo.CDF(x,Eta,sEta,sB,W)
}

\references{

Cheng, Y. S. and Peng, C. Y. (2012), "Integrated Degradation Models in R Using iDEMO." \emph{Journal of Statistical Software}, \bold{49(2)}, 1-22.

Peng, C. Y. and Tseng, S. T. (2009), "Mis-specification Analysis of Linear Degradation Models." \emph{IEEE Transactions on Reliability}, \bold{58}, 444-455.

}

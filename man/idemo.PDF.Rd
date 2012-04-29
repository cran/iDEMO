\name{idemo.PDF}
\alias{idemo.PDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability density function of the lifetime distribution.
}
\description{
The function idemo.PDF computes the probability density function of the lifetime distribution.
}
\usage{
idemo.PDF(t, Eta, sEta, sB, W)
}

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
\code{\link{idemo.CDF}},\code{\link{idemo.MTTF}}
}
\examples{
t <- 5000
Eta <- 0.00204
sEta <- 0.00042
sB <- 0.0104
W <- 10
idemo.PDF(t,Eta,sEta,sB,W)
}

\references{

Cheng, Y. S. and Peng, C. Y. (2012), "Integrated Degradation Models in R Using iDEMO." \emph{Journal of Statistical Software}, \bold{49(2)}, 1-22.

Peng, C. Y. and Tseng, S. T. (2009), "Mis-specification Analysis of Linear Degradation Models." \emph{IEEE Transactions on Reliability}, \bold{58}, 444-455.

}

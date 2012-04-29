\name{OIM.score}
\alias{OIM.score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Observed information matrix calculated by score vector of the degradation model
}
\description{
The function OIM.score computes the observed information matrix calculated by score vector of the degradation model.

}
\usage{
OIM.score(Data, Eta, sEta, sB, sE)
}

\arguments{
  \item{Data}{
the data set.
}
  \item{Eta}{
the drift rate of the degradation model.
}
  \item{sEta}{
the standard deviation of the random effect.
}
  \item{sB}{
the diffusion coefficient of the Brownian motion.
}
  \item{sE}{
the standard deviation of the measurement error.
}
}

\details{
The function OIM.score returns a symmetric 4 by 4 matrix if all the third to fifth 
arguments in the function OIM.score are greater than zero. The order from the first row (column) 
to the last row (column) in the OIM is Eta, sEta^2, sB^2 and sE^2. If any one of the 
third to fifth parameters in the function OIM.score is zero, the corresponding
information matrix will be reduced to a smaller matrix.
}


\seealso{
\code{\link{idemo.FIM}},\code{\link{OIM.hessian}},\code{\link{OIM.robust}}
}

\examples{
library(iDEMO)
data(laser)
Eta <- 0.00204
sEta <- 0.00042
sB <- 0.0104
sE <- 0.0309
OIM.score(laser, Eta, sEta, sB, sE) #full model
OIM.score(laser, Eta, sEta, 0, sE)  #full model without time correlated structure
OIM.score(laser, Eta, 0, sB, 0)     #full model without random effect and measurement error
OIM.score(laser, Eta, sEta, sB, 0)  #full model without measurement error
}


 

\references{

Boldea, O. and Magnus J. R. (2009), "Maximum Likelihood Estimation of the Multivariate Normal Mixture Model." \emph{Journal of the American Statistical Association}, \bold{104}, 1539-1549.

Cheng, Y. S. and Peng, C. Y. (2012), "Integrated Degradation Models in R Using iDEMO." \emph{Journal of Statistical Software}, \bold{49(2)}, 1-22.

Peng, C. Y. and Tseng, S. T. (2009), "Mis-specification Analysis of Linear Degradation Models." \emph{IEEE Transactions on Reliability}, \bold{58}, 444-455.

}


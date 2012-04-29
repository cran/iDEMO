\name{idemo.FIM}
\alias{idemo.FIM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fisher's information matrix of the degradation model
}
\description{
The function idemo.FIM computes the Fisher's information matrix of the degradation model.

}
\usage{
idemo.FIM(sEta, sB, sE, Time)
}
\arguments{
  \item{sEta}{
the standard deviation of the random effect.
}
  \item{sB}{
the diffusion coefficient of the Brownian motion.
}
  \item{sE}{
the standard deviation of the measurement error.
}
  \item{Time}{
the measuring time.
}
}


\details{
The function idemo.FIM returns a symmetric 4 by 4 matrix if all the first 3 
arguments in the function idemo.FIM are greater than zero. Each element is 
the negative of the expectation of the second derivative of the log-likelihood 
function with respect to the parameter. The order from the first row (column) 
to the last row (column) in the FIM is Eta, sEta^2, sB^2 and sE^2. If any one of the 
first three parameters in the function idemo.FIM is zero, the corresponding Fisher's
information matrix will be reduced to a smaller matrix.
}

\seealso{
\code{\link{OIM.hessian}},\code{\link{OIM.score}},\code{\link{OIM.robust}}
}

\examples{
sEta <- 0.00042
sB <- 0.0104
sE <- 0.0309
Time <- seq(250,4000,250)
idemo.FIM(sEta,sB,sE,Time) #full model
idemo.FIM(sEta,0,sE,Time)  #full model without time correlated structure
idemo.FIM(0,sB,0,Time)     #full model without random effect and measurement error
idemo.FIM(sEta,sB,0,Time)  #full model without measurement error
}

\references{

Cheng, Y. S. and Peng, C. Y. (2012), "Integrated Degradation Models in R Using iDEMO." \emph{Journal of Statistical Software}, \bold{49(2)}, 1-22.

Peng, C. Y. and Tseng, S. T. (2009), "Mis-specification Analysis of Linear Degradation Models." \emph{IEEE Transactions on Reliability}, \bold{58}, 444-455.

}

model.selection <-
function(loglik, q, n, ...){

	AIC <- -2 * loglik + 2 * q

	BIC <- -2 * loglik + q * log(n)

	HQC <- -2 * loglik + 2 * q * log(log(n))

	return(c(AIC, BIC, HQC))

}


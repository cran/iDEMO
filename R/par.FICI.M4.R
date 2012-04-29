par.FICI.M4 <-
function(alpha, Par.M4, Cov.Mat.M4, n, ...){

	Eta.M4 <- Par.M4[1]

	sE.M4 <- Par.M4[2]

	CI.lower.eta <- Eta.M4 - qnorm( 1 - alpha / 2 ) * sqrt( Cov.Mat.M4[1,1] / n )

	CI.upper.eta <- Eta.M4 + qnorm( 1 - alpha / 2 ) * sqrt( Cov.Mat.M4[1,1] / n )

	if(Eta.M4>0){
		CI.ln.lower.eta <- exp( log(Eta.M4) - qnorm( 1 - alpha / 2 ) * sqrt( Eta.M4^(-2) * Cov.Mat.M4[1,1] / n ) )
		CI.ln.upper.eta <- exp( log(Eta.M4) + qnorm( 1 - alpha / 2 ) * sqrt( Eta.M4^(-2) * Cov.Mat.M4[1,1] / n ) )
	}else{
		CI.ln.lower.eta <- NA
		CI.ln.upper.eta <- NA
	}

	var.sE <- ( 2 * sE.M4 )^(-2) * Cov.Mat.M4[2,2]

	CI.lower.sE <- sE.M4 - qnorm( 1 - alpha / 2 ) * sqrt( var.sE / n )

	CI.upper.sE <- sE.M4 + qnorm( 1 - alpha / 2 ) * sqrt( var.sE / n )

	CI.ln.lower.sE <- exp( log(sE.M4) - qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sE.M4^2 )^(-2) * Cov.Mat.M4[2,2] / n ) )

	CI.ln.upper.sE <- exp( log(sE.M4) + qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sE.M4^2 )^(-2) * Cov.Mat.M4[2,2] / n ) )

	list( CI.Eta.M4 = c(CI.lower.eta, CI.upper.eta), 

		CI.ln.Eta.M4 = c(CI.ln.lower.eta, CI.ln.upper.eta), 

		CI.sE.M4 = c(CI.lower.sE, CI.upper.sE), 

		CI.ln.sE.M4 = c(CI.ln.lower.sE, CI.ln.upper.sE) )

}


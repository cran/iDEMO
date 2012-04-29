par.FICI.M2 <-
function(alpha, Par.M2, Cov.Mat.M2, n, ...){

	Eta.M2 <- Par.M2[1]

	sB.M2 <- Par.M2[2]

	CI.lower.eta <- Eta.M2 - qnorm( 1 - alpha / 2 ) * sqrt( Cov.Mat.M2[1,1] / n )

	CI.upper.eta <- Eta.M2 + qnorm( 1 - alpha / 2 ) * sqrt( Cov.Mat.M2[1,1] / n )

	if(Eta.M2>0){
		CI.ln.lower.eta<- exp( log(Eta.M2) - qnorm( 1 - alpha / 2 ) * sqrt( Eta.M2^(-2) * Cov.Mat.M2[1,1] / n ) )
		CI.ln.upper.eta<- exp( log(Eta.M2) + qnorm( 1 - alpha / 2 ) * sqrt( Eta.M2^(-2) * Cov.Mat.M2[1,1] / n ) )
	}else{
		CI.ln.lower.eta<- NA
		CI.ln.upper.eta<- NA
	}

	var.sB <- ( 2 * sB.M2 )^(-2) * Cov.Mat.M2[2,2]

	CI.lower.sB <- sB.M2 - qnorm( 1 - alpha / 2 ) * sqrt( var.sB / n )

	CI.upper.sB <- sB.M2 + qnorm( 1 - alpha / 2 ) * sqrt( var.sB / n )

	CI.ln.lower.sB <- exp( log(sB.M2) - qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sB.M2^2 )^(-2) * Cov.Mat.M2[2,2] / n ) )

	CI.ln.upper.sB <- exp( log(sB.M2) + qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sB.M2^2 )^(-2) * Cov.Mat.M2[2,2] / n ) )

	list( CI.Eta.M2 = c(CI.lower.eta, CI.upper.eta), 

		CI.ln.Eta.M2 = c(CI.ln.lower.eta, CI.ln.upper.eta), 

		CI.sB.M2 = c(CI.lower.sB, CI.upper.sB), 

		CI.ln.sB.M2 = c(CI.ln.lower.sB, CI.ln.upper.sB) )

}


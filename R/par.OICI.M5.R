par.OICI.M5 <-
function(alpha, Par.M5, Cov.Mat.M5, ...){

	Eta.M5 <- Par.M5[1]

	sB.M5 <- Par.M5[2]

	sE.M5 <- Par.M5[3]

	CI.lower.eta <- Eta.M5 - qnorm( 1 - alpha / 2 ) * sqrt( Cov.Mat.M5[1,1] )

	CI.upper.eta <- Eta.M5 + qnorm( 1 - alpha / 2 ) * sqrt( Cov.Mat.M5[1,1] )

	if(Eta.M5>0){
		CI.ln.lower.eta<- exp( log(Eta.M5) - qnorm( 1 - alpha / 2 ) * sqrt( Eta.M5^(-2) * Cov.Mat.M5[1,1] ) )
		CI.ln.upper.eta<- exp( log(Eta.M5) + qnorm( 1 - alpha / 2 ) * sqrt( Eta.M5^(-2) * Cov.Mat.M5[1,1] ) )
	}else{
		CI.ln.lower.eta<- NA
		CI.ln.upper.eta<- NA
	}
	var.sB <- ( 2 * sB.M5 )^(-2) * Cov.Mat.M5[2,2]

	CI.lower.sB <- sB.M5 - qnorm( 1 - alpha / 2 ) * sqrt( var.sB )

	CI.upper.sB <- sB.M5 + qnorm( 1 - alpha / 2 ) * sqrt( var.sB )

	CI.ln.lower.sB <- exp( log(sB.M5) - qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sB.M5^2 )^(-2) * Cov.Mat.M5[2,2] ) )

	CI.ln.upper.sB <- exp( log(sB.M5) + qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sB.M5^2 )^(-2) * Cov.Mat.M5[2,2] ) )

	var.sE <- ( 2 * sE.M5 )^(-2) * Cov.Mat.M5[3,3]

	CI.lower.sE <- sE.M5 - qnorm( 1 - alpha / 2 ) * sqrt( var.sE )

	CI.upper.sE <- sE.M5 + qnorm( 1 - alpha / 2 ) * sqrt( var.sE )

	CI.ln.lower.sE <- exp( log(sE.M5) - qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sE.M5^2 )^(-2) * Cov.Mat.M5[3,3] ) )

	CI.ln.upper.sE <- exp( log(sE.M5) + qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sE.M5^2 )^(-2) * Cov.Mat.M5[3,3] ) )

	list( CI.Eta.M5 = c(CI.lower.eta, CI.upper.eta), 

		CI.ln.Eta.M5 = c(CI.ln.lower.eta, CI.ln.upper.eta), 

		CI.sB.M5 = c(CI.lower.sB, CI.upper.sB), 

		CI.ln.sB.M5 = c(CI.ln.lower.sB, CI.ln.upper.sB), 

		CI.sE.M5 = c(CI.lower.sE, CI.upper.sE), 

		CI.ln.sE.M5 = c(CI.ln.lower.sE, CI.ln.upper.sE) )

}


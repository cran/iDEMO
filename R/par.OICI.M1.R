par.OICI.M1 <-
function(alpha, Par.M1, Cov.Mat.M1, ...){

	Eta.M1 <- Par.M1[1]

	sEta.M1 <- Par.M1[2]

	sE.M1 <- Par.M1[3]

	CI.lower.eta <- Eta.M1 - qnorm( 1 - alpha / 2 ) * sqrt( Cov.Mat.M1[1,1] )

	CI.upper.eta <- Eta.M1 + qnorm( 1 - alpha / 2 ) * sqrt( Cov.Mat.M1[1,1] )

	if(Eta.M1>0){
		CI.ln.lower.eta<- exp( log(Eta.M1) - qnorm( 1 - alpha / 2 ) * sqrt( Eta.M1^(-2) * Cov.Mat.M1[1,1] ) )
		CI.ln.upper.eta<- exp( log(Eta.M1) + qnorm( 1 - alpha / 2 ) * sqrt( Eta.M1^(-2) * Cov.Mat.M1[1,1] ) )
	}else{
		CI.ln.lower.eta<- NA
		CI.ln.upper.eta<- NA
	}
	var.s.eta <- ( 2 * sEta.M1 )^(-2) * Cov.Mat.M1[2,2]

	CI.lower.s.eta <- sEta.M1 - qnorm( 1 - alpha / 2 ) * sqrt( var.s.eta )

	CI.upper.s.eta <- sEta.M1 + qnorm( 1 - alpha / 2 ) * sqrt( var.s.eta )

	CI.ln.lower.s.eta <- exp( log(sEta.M1) - qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sEta.M1^2 )^(-2) * Cov.Mat.M1[2,2] ) )

	CI.ln.upper.s.eta <- exp( log(sEta.M1) + qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sEta.M1^2 )^(-2) * Cov.Mat.M1[2,2] ) )

	var.sE <- ( 2 * sE.M1 )^(-2) * Cov.Mat.M1[3,3]

	CI.lower.sE <- sE.M1 - qnorm( 1 - alpha / 2 ) * sqrt( var.sE )

	CI.upper.sE <- sE.M1 + qnorm( 1 - alpha / 2 ) * sqrt( var.sE )

	CI.ln.lower.sE <- exp( log(sE.M1) - qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sE.M1^2 )^(-2) * Cov.Mat.M1[3,3] ) )

	CI.ln.upper.sE <- exp( log(sE.M1) + qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sE.M1^2 )^(-2) * Cov.Mat.M1[3,3] ) )

	list( CI.Eta.M1 = c(CI.lower.eta, CI.upper.eta), 

		CI.ln.Eta.M1 = c(CI.ln.lower.eta, CI.ln.upper.eta), 

		CI.sEta.M1 = c(CI.lower.s.eta, CI.upper.s.eta), 

		CI.ln.sEta.M1 = c(CI.ln.lower.s.eta, CI.ln.upper.s.eta), 

		CI.sE.M1 = c(CI.lower.sE, CI.upper.sE), 

		CI.ln.sE.M1 = c(CI.ln.lower.sE, CI.ln.upper.sE) )

}


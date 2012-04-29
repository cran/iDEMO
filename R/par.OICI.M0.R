par.OICI.M0 <-
function(alpha, Par.M0, Cov.Mat.M0, ...){

	Eta.M0 <- Par.M0[1]

	sEta.M0 <- Par.M0[2]

	sB.M0 <- Par.M0[3]

	sE.M0 <- Par.M0[4]

	CI.lower.eta <- Eta.M0 - qnorm( 1 - alpha / 2 ) * sqrt( Cov.Mat.M0[1,1] )

	CI.upper.eta <- Eta.M0 + qnorm( 1 - alpha / 2 ) * sqrt( Cov.Mat.M0[1,1] )

	if(Eta.M0>0){
		CI.ln.lower.eta<- exp( log(Eta.M0) - qnorm( 1 - alpha / 2 ) * sqrt( Eta.M0^(-2) * Cov.Mat.M0[1,1] ) )
		CI.ln.upper.eta<- exp( log(Eta.M0) + qnorm( 1 - alpha / 2 ) * sqrt( Eta.M0^(-2) * Cov.Mat.M0[1,1] ) )
	}else{
		CI.ln.lower.eta<- NA
		CI.ln.upper.eta<- NA
	}
	var.s.eta <- ( 2 * sEta.M0 )^(-2) * Cov.Mat.M0[2,2]

	CI.lower.s.eta <- sEta.M0 - qnorm( 1 - alpha / 2 ) * sqrt( var.s.eta )

	CI.upper.s.eta <- sEta.M0 + qnorm( 1 - alpha / 2 ) * sqrt( var.s.eta )

	CI.ln.lower.s.eta <- exp( log(sEta.M0) - qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sEta.M0^2 )^(-2) * Cov.Mat.M0[2,2] ) )

	CI.ln.upper.s.eta <- exp( log(sEta.M0) + qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sEta.M0^2 )^(-2) * Cov.Mat.M0[2,2] ) )

	var.sB <- ( 2 * sB.M0 )^(-2) * Cov.Mat.M0[3,3]

	CI.lower.sB <- sB.M0 - qnorm( 1 - alpha / 2 ) * sqrt( var.sB )

	CI.upper.sB <- sB.M0 + qnorm( 1 - alpha / 2 ) * sqrt( var.sB )

	CI.ln.lower.sB <- exp( log(sB.M0) - qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sB.M0^2 )^(-2) * Cov.Mat.M0[3,3] ) )

	CI.ln.upper.sB <- exp( log(sB.M0) + qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sB.M0^2 )^(-2) * Cov.Mat.M0[3,3] ) )

	var.sE <- ( 2 * sE.M0 )^(-2) * Cov.Mat.M0[4,4]

	CI.lower.sE <- sE.M0 - qnorm( 1 - alpha / 2 ) * sqrt( var.sE )

	CI.upper.sE <- sE.M0 + qnorm( 1 - alpha / 2 ) * sqrt( var.sE )

	CI.ln.lower.sE <- exp( log(sE.M0) - qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sE.M0^2 )^(-2) * Cov.Mat.M0[4,4] ) )

	CI.ln.upper.sE <- exp( log(sE.M0) + qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sE.M0^2 )^(-2) * Cov.Mat.M0[4,4] ) )

	list( CI.Eta.M0 = c(CI.lower.eta, CI.upper.eta), 

		CI.ln.Eta.M0 = c(CI.ln.lower.eta, CI.ln.upper.eta), 

		CI.sEta.M0 = c(CI.lower.s.eta, CI.upper.s.eta), 

		CI.ln.sEta.M0 = c(CI.ln.lower.s.eta, CI.ln.upper.s.eta), 

		CI.sB.M0 = c(CI.lower.sB, CI.upper.sB), 

		CI.ln.sB.M0 = c(CI.ln.lower.sB, CI.ln.upper.sB), 

		CI.sE.M0 = c(CI.lower.sE, CI.upper.sE), 

		CI.ln.sE.M0 = c(CI.ln.lower.sE, CI.ln.upper.sE) )


}


MTTF.OICI.M5 <-
function(alpha, W, Eta.M5, Cov.Mat.M5, ...){

	MTTF.approx.M5 <- MTTF.exact.M5 <- W / Eta.M5

	var.MTTF.approx.OI.M5 <- ( -W / Eta.M5^2 )^2 * Cov.Mat.M5[1,1]

	OICI.lower.MTTF.approx.M5 <- OICI.lower.MTTF.exact.M5 <- W / Eta.M5 - qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.OI.M5  )

	OICI.upper.MTTF.approx.M5 <- OICI.upper.MTTF.exact.M5 <- W / Eta.M5 + qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.OI.M5  )


	OICI.ln.lower.MTTF.approx.M5 <- OICI.ln.lower.MTTF.exact.M5 <- exp( log( W / Eta.M5 ) - qnorm( 1 - alpha / 2 ) * sqrt( Eta.M5^(-2) * Cov.Mat.M5[1,1]  ) )

	OICI.ln.upper.MTTF.approx.M5 <- OICI.ln.upper.MTTF.exact.M5 <- exp( log( W / Eta.M5 ) + qnorm( 1 - alpha / 2 ) * sqrt( Eta.M5^(-2) * Cov.Mat.M5[1,1]  ) )


	list( MTTF.approx.M5 = MTTF.approx.M5, 

		CI.MTTF.approx.M5 = c(OICI.lower.MTTF.approx.M5, OICI.upper.MTTF.approx.M5), 

		CI.ln.MTTF.approx.M5 = c(OICI.ln.lower.MTTF.approx.M5, OICI.ln.upper.MTTF.approx.M5),

		MTTF.exact.M5 = MTTF.exact.M5, 

		CI.MTTF.exact.M5 = c(OICI.lower.MTTF.exact.M5, OICI.upper.MTTF.exact.M5), 

		CI.ln.MTTF.exact.M5 = c(OICI.ln.lower.MTTF.exact.M5, OICI.ln.upper.MTTF.exact.M5) )

}


MTTF.OICI.M2 <-
function(alpha, W, Eta.M2, Cov.Mat.M2, ...){

	MTTF.approx.M2 <- MTTF.exact.M2 <- W / Eta.M2

	var.MTTF.approx.OI.M2 <- ( -W / Eta.M2^2 )^2 * Cov.Mat.M2[1,1]

	OICI.lower.MTTF.approx.M2 <- OICI.lower.MTTF.exact.M2 <- W / Eta.M2 - qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.OI.M2  )

	OICI.upper.MTTF.approx.M2 <- OICI.upper.MTTF.exact.M2 <- W / Eta.M2 + qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.OI.M2  )


	OICI.ln.lower.MTTF.approx.M2 <- OICI.ln.lower.MTTF.exact.M2 <- exp( log( W / Eta.M2 ) - qnorm( 1 - alpha / 2 ) * sqrt( Eta.M2^(-2) * Cov.Mat.M2[1,1]  ) )

	OICI.ln.upper.MTTF.approx.M2 <- OICI.ln.upper.MTTF.exact.M2 <- exp( log( W / Eta.M2 ) + qnorm( 1 - alpha / 2 ) * sqrt( Eta.M2^(-2) * Cov.Mat.M2[1,1]  ) )


	list( MTTF.approx.M2 = MTTF.approx.M2, 

		CI.MTTF.approx.M2 = c(OICI.lower.MTTF.approx.M2, OICI.upper.MTTF.approx.M2), 

		CI.ln.MTTF.approx.M2 = c(OICI.ln.lower.MTTF.approx.M2, OICI.ln.upper.MTTF.approx.M2),

		MTTF.exact.M2 = MTTF.exact.M2, 

		CI.MTTF.exact.M2 = c(OICI.lower.MTTF.exact.M2, OICI.upper.MTTF.exact.M2), 

		CI.ln.MTTF.exact.M2 = c(OICI.ln.lower.MTTF.exact.M2, OICI.ln.upper.MTTF.exact.M2) )

}


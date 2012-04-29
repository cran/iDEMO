MTTF.FICI.M5 <-
function(alpha, W, Eta.M5, Cov.Mat.M5, n, ...){

	MTTF.approx.M5 <- MTTF.exact.M5 <- W / Eta.M5

	var.MTTF.approx.FI.M5 <- ( -W / Eta.M5^2 )^2 * Cov.Mat.M5[1,1]

	FICI.lower.MTTF.approx.M5 <- FICI.lower.MTTF.exact.M5 <- W / Eta.M5 - qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.FI.M5 / n )

	FICI.upper.MTTF.approx.M5 <- FICI.upper.MTTF.exact.M5 <- W / Eta.M5 + qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.FI.M5 / n )


	FICI.ln.lower.MTTF.approx.M5 <- FICI.ln.lower.MTTF.exact.M5 <- exp( log( W / Eta.M5 ) - qnorm( 1 - alpha / 2 ) * sqrt( Eta.M5^(-2) * Cov.Mat.M5[1,1] / n ) )

	FICI.ln.upper.MTTF.approx.M5 <- FICI.ln.upper.MTTF.exact.M5 <- exp( log( W / Eta.M5 ) + qnorm( 1 - alpha / 2 ) * sqrt( Eta.M5^(-2) * Cov.Mat.M5[1,1] / n ) )


	list( MTTF.approx.M5 = MTTF.approx.M5, 

		CI.MTTF.approx.M5 = c(FICI.lower.MTTF.approx.M5, FICI.upper.MTTF.approx.M5), 

		CI.ln.MTTF.approx.M5 = c(FICI.ln.lower.MTTF.approx.M5, FICI.ln.upper.MTTF.approx.M5),

		MTTF.exact.M5 = MTTF.exact.M5, 

		CI.MTTF.exact.M5 = c(FICI.lower.MTTF.exact.M5, FICI.upper.MTTF.exact.M5), 

		CI.ln.MTTF.exact.M5 = c(FICI.ln.lower.MTTF.exact.M5, FICI.ln.upper.MTTF.exact.M5) )

}


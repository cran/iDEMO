MTTF.FICI.M2 <-
function(alpha, W, Eta.M2, Cov.Mat.M2, n, ...){

	MTTF.approx.M2 <- MTTF.exact.M2 <- W / Eta.M2

	var.MTTF.approx.FI.M2 <- ( -W / Eta.M2^2 )^2 * Cov.Mat.M2[1,1]

	FICI.lower.MTTF.approx.M2 <- FICI.lower.MTTF.exact.M2 <- W / Eta.M2 - qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.FI.M2 / n )

	FICI.upper.MTTF.approx.M2 <- FICI.upper.MTTF.exact.M2 <- W / Eta.M2 + qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.FI.M2 / n )


	FICI.ln.lower.MTTF.approx.M2 <- FICI.ln.lower.MTTF.exact.M2 <- exp( log( W / Eta.M2 ) - qnorm( 1 - alpha / 2 ) * sqrt( Eta.M2^(-2) * Cov.Mat.M2[1,1] / n ) )

	FICI.ln.upper.MTTF.approx.M2 <- FICI.ln.upper.MTTF.exact.M2 <- exp( log( W / Eta.M2 ) + qnorm( 1 - alpha / 2 ) * sqrt( Eta.M2^(-2) * Cov.Mat.M2[1,1] / n ) )


	list( MTTF.approx.M2 = MTTF.approx.M2, 

		CI.MTTF.approx.M2 = c(FICI.lower.MTTF.approx.M2, FICI.upper.MTTF.approx.M2), 

		CI.ln.MTTF.approx.M2 = c(FICI.ln.lower.MTTF.approx.M2, FICI.ln.upper.MTTF.approx.M2),

		MTTF.exact.M2 = MTTF.exact.M2, 

		CI.MTTF.exact.M2 = c(FICI.lower.MTTF.exact.M2, FICI.upper.MTTF.exact.M2), 

		CI.ln.MTTF.exact.M2 = c(FICI.ln.lower.MTTF.exact.M2, FICI.ln.upper.MTTF.exact.M2) )

}


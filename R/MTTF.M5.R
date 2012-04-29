MTTF.M5 <-
function( Eta.M5, sB.M5, W, ...){

	MTTF.exact.M5 <- MTTF.approx.M5 <- W / Eta.M5

	return( c( MTTF.approx.M5, MTTF.exact.M5 ) )

}


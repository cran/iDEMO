MTTF.M2 <-
function( Eta.M2, sB.M2, W, ...){

	MTTF.exact.M2 <- MTTF.approx.M2 <- W / Eta.M2

	return( c( MTTF.approx.M2, MTTF.exact.M2 ) )

}


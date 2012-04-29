MTTF.M0 <-
function( Eta.M0, sEta.M0, W, ...){

	#approx.

	MTTF.approx.M0 <- W / Eta.M0

	#Exact

	f1 <- function(x)  exp( x^2 / ( 2 * sEta.M0^2 ) )

	if( (Eta.M0/sqrt(2)/sEta.M0) > 10^307 ){

		MTTF.exact.M0 <- W / Eta.M0

	}else{

		MTTF.exact.M0 <- sqrt(2)*W / sEta.M0 * dawson(Eta.M0/sqrt(2)/sEta.M0) #W / sEta.M0^2 * exp( - Eta.M0^2 / ( 2 * sEta.M0^2 ) ) * integrate( f1, 0, Eta.M0 )$value

	}

	return( c( MTTF.approx.M0, MTTF.exact.M0 ) )

}


MTTF.M1 <-
function( Eta.M1, sEta.M1, W, ...){

	#approx.

	MTTF.approx.M1 <- W / Eta.M1


	#Exact

	f1 <- function(x)  exp( x^2 / ( 2 * sEta.M1^2 ) )

	if( (Eta.M1/sqrt(2)/sEta.M1) > 10^307 ){

		MTTF.exact.M1 <- W / Eta.M1

	}else{

		MTTF.exact.M1 <- sqrt(2)*W / sEta.M1 * dawson(Eta.M1/sqrt(2)/sEta.M1) #W / sEta.M1^2 * exp( - Eta.M1^2 / ( 2 * sEta.M1^2 ) ) * integrate( f1, 0, Eta.M1 )$value

	}

	return( c( MTTF.approx.M1, MTTF.exact.M1 ) )

}


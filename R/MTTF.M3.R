MTTF.M3 <-
function( Eta.M3, sEta.M3, W, ...){

    #approx.

    MTTF.approx.M3 <- W / Eta.M3

    #Exact

    f1 <- function(x)  exp( x^2 / ( 2 * sEta.M3^2 ) )

    if( (Eta.M3/sqrt(2)/sEta.M3) > 10^307 ){

	    MTTF.exact.M3 <- W / Eta.M3

    }else{

	    MTTF.exact.M3 <- sqrt(2)*W / sEta.M3 * dawson(Eta.M3/sqrt(2)/sEta.M3) #W / sEta.M3^2 * exp( - Eta.M3^2 / ( 2 * sEta.M3^2 ) ) * integrate( f1, 0, Eta.M3 )$value

    }

    return( c( MTTF.approx.M3, MTTF.exact.M3 ) )

}


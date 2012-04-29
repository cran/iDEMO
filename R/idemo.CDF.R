idemo.CDF <-
function( t, Eta, sEta, sB, W ) {
	if( is.numeric(Eta) & is.numeric(sEta) & is.numeric(sB) & is.numeric(W) & is.numeric(t)){
		if( sEta>=0 & sB>=0 && all(t>=0) ){
			if(sB==0){
				idemo.CDF <- pnorm( ( Eta * t - W ) / ( sEta * t ) )
			}else{
				idemo.CDF <- pnorm( (  Eta * t - W )  / sqrt(  sB^2 * t + sEta^2 * t^2 ) )  + 
      		     	exp(  2 * Eta * W / sB^2 )  * pnorm(   - (  2 * sEta^2 * W * t + sB^2 * (  Eta * t + W ) )  / 
	           		(  sB^2 * sqrt(  sB^2 * t + sEta^2 * t^2 ) ) )  * 
		           	exp(  2 * sEta^2 * W^2 / sB^4 / 8 )  * exp(  2 * sEta^2 * W^2 / sB^4 / 8 )  * 
		           	exp(  2 * sEta^2 * W^2 / sB^4 / 8 )  * exp(  2 * sEta^2 * W^2 / sB^4 / 8 )  * 
		           	exp(  2 * sEta^2 * W^2 / sB^4 / 8 )  * exp(  2 * sEta^2 * W^2 / sB^4 / 8 )  * 
		           	exp(  2 * sEta^2 * W^2 / sB^4 / 8 )  * exp(  2 * sEta^2 * W^2 / sB^4 / 8 ) 
			}
			return(idemo.CDF)
		}else{
			cat("The parameter and time should not be a negative value and vector, respectively!","\n")
		}
	}else{
		cat("The parameter and time should not be a negative value and vector, respectively!","\n")
	}
}


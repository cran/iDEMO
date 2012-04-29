idemo.MTTF <-
function( Eta, sEta, W){
	if( is.numeric(Eta) & is.numeric(sEta) & is.numeric(W)){
		if( sEta>=0 ){
			if(sEta<1e-10){
				idemo.MTTF <- W / Eta
			}else{
				f1 <- function(x)  exp( x^2 / ( 2 * sEta^2 ) )
				idemo.MTTF <- sqrt(2)*W / sEta * dawson(Eta/sqrt(2)/sEta)
			}
			return( idemo.MTTF )
		}else{
			cat("The parameter should not be negative!","\n")
		}
	}else{
		cat("The parameter should not be negative!","\n")
	}
}


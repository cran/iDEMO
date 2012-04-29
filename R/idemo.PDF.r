idemo.PDF <-
function( t, Eta, sEta, sB, W ) {
      if( is.numeric(Eta) & is.numeric(sEta) & is.numeric(sB) & is.numeric(W) & is.numeric(t)){
		if( sEta>=0 & sB>=0 && all(t>=0) ){
			idemo.PDF <- sqrt( W^2 / ( 2 * pi * t^3 * ( sEta^2 * t + sB^2 ) ) ) *
			exp( - ( W - Eta * t )^2 / ( 2 * t * ( sEta^2 * t + sB^2 ) ) )
			return(idemo.PDF)
		}else{
			cat("The parameter and time should not be a negative value and vector, respectively!","\n")
		}
	}else{
		cat("The parameter and time should not be a negative value and vector, respectively!","\n")
      }
}


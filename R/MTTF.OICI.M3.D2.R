MTTF.OICI.M3.D2 <-
function(alpha, W, Eta.M3, sEta.M3, Cov.Mat.M3, h.val, ...){

	if( h.val > min(Eta.M3/2, sEta.M3^2/2) ) h.val <- min(Eta.M3/3, sEta.M3^2/3)

	#approx.

	MTTF.approx.M3 <- W / Eta.M3

	var.MTTF.approx.M3 <- ( -W / Eta.M3^2 )^2 * Cov.Mat.M3[1,1]

	OICI.lower.MTTF.approx.M3 <- MTTF.approx.M3 - qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.M3 )

	OICI.upper.MTTF.approx.M3 <- MTTF.approx.M3 + qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.M3 )


	OICI.ln.lower.MTTF.approx.M3 <- exp( log( MTTF.approx.M3 ) - qnorm( 1 - alpha / 2 ) * sqrt( Eta.M3^(-2) * Cov.Mat.M3[1,1] ) )

	OICI.ln.upper.MTTF.approx.M3 <- exp( log( MTTF.approx.M3 ) + qnorm( 1 - alpha / 2 ) * sqrt( Eta.M3^(-2) * Cov.Mat.M3[1,1] ) )


	#Exact

	if( (Eta.M3/sqrt(2)/sEta.M3) <= 10^307 ){

		f1 <- function(x)  exp( x^2 / ( 2 * sEta.M3^2 ) )

		MTTF.exact.M3 <- sqrt(2)*W / sEta.M3 * dawson(Eta.M3/sqrt(2)/sEta.M3) #W / sEta.M3^2 * exp( - Eta.M3^2 / ( 2 * sEta.M3^2 ) ) * integrate( f1, 0, Eta.M3 )$value

		diff_MTTF.M3 <- numeric( length = 3 )

			sEta2.M3 <- sEta.M3^2

			MTTF.M3.Eta <- function(eta) sqrt(2) * W / sEta.M3 * dawson( eta / sqrt(2) / sEta.M3 )

			MTTF.M3.sEta2 <- function(seta2) sqrt(2) * W / sqrt(seta2) * dawson( Eta.M3 / sqrt(2) / sqrt(seta2) )

			diff_MTTF.M3[1] <- ( MTTF.M3.Eta( Eta.M3 - 2 * h.val ) - 8 * MTTF.M3.Eta( Eta.M3 - h.val ) +

						8 * MTTF.M3.Eta( Eta.M3 + h.val ) - MTTF.M3.Eta( Eta.M3 + 2 * h.val ) ) / 

						( 12 * h.val )

			diff_MTTF.M3[2] <- ( MTTF.M3.sEta2( sEta2.M3 - 2 * h.val ) - 8 * MTTF.M3.sEta2( sEta2.M3 - h.val ) +

						8 * MTTF.M3.sEta2( sEta2.M3 + h.val ) - MTTF.M3.sEta2( sEta2.M3 + 2 * h.val ) ) / 

						( 12 * h.val )

		var.MTTF.exact.M3 <- t(diff_MTTF.M3) %*% Cov.Mat.M3 %*% diff_MTTF.M3

		OICI.lower.MTTF.exact.M3 <- MTTF.exact.M3 - qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.exact.M3 )

		OICI.upper.MTTF.exact.M3 <- MTTF.exact.M3 + qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.exact.M3 )

	
		var2.MTTF.exact <- t( diff_MTTF.M3 / MTTF.exact.M3 ) %*% Cov.Mat.M3 %*% ( diff_MTTF.M3 / MTTF.exact.M3 )

		OICI.ln.lower.MTTF.exact.M3 <- exp( log( MTTF.exact.M3 ) - qnorm( 1 - alpha / 2 ) * sqrt( var2.MTTF.exact ) )

		OICI.ln.upper.MTTF.exact.M3 <- exp( log( MTTF.exact.M3 ) + qnorm( 1 - alpha / 2 ) * sqrt( var2.MTTF.exact ) )

	}else{
		MTTF.exact.M3 <- MTTF.approx.M3
		OICI.lower.MTTF.exact.M3 <- OICI.lower.MTTF.approx.M3
		OICI.upper.MTTF.exact.M3 <- OICI.upper.MTTF.approx.M3
		OICI.ln.lower.MTTF.exact.M3 <- OICI.ln.lower.MTTF.approx.M3
		OICI.ln.upper.MTTF.exact.M3 <- OICI.ln.upper.MTTF.approx.M3
	}
	
	list( MTTF.approx.M3 = MTTF.approx.M3, 

		CI.MTTF.approx.M3 = c(OICI.lower.MTTF.approx.M3, OICI.upper.MTTF.approx.M3), 

		CI.ln.MTTF.approx.M3 = c(OICI.ln.lower.MTTF.approx.M3, OICI.ln.upper.MTTF.approx.M3),

		MTTF.exact.M3 = MTTF.exact.M3, 

		CI.MTTF.exact.M3 = c(OICI.lower.MTTF.exact.M3, OICI.upper.MTTF.exact.M3), 

		CI.ln.MTTF.exact.M3 = c(OICI.ln.lower.MTTF.exact.M3, OICI.ln.upper.MTTF.exact.M3) )

}


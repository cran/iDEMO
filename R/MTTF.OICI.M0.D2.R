MTTF.OICI.M0.D2 <-
function(alpha, W, Eta.M0, sEta.M0, Cov.Mat.M0, h.val, ...){

	if( h.val > min(Eta.M0/2, sEta.M0^2/2) ) h.val <- min(Eta.M0/3, sEta.M0^2/3)

	#approx.

	MTTF.approx.M0 <- W / Eta.M0

	var.MTTF.approx.M0 <- ( -W / Eta.M0^2 )^2 * Cov.Mat.M0[1,1]

	OICI.lower.MTTF.approx.M0 <- MTTF.approx.M0 - qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.M0 )

	OICI.upper.MTTF.approx.M0 <- MTTF.approx.M0 + qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.M0 )


	OICI.ln.lower.MTTF.approx.M0 <- exp( log( MTTF.approx.M0 ) - qnorm( 1 - alpha / 2 ) * sqrt( Eta.M0^(-2) * Cov.Mat.M0[1,1] ) )

	OICI.ln.upper.MTTF.approx.M0 <- exp( log( MTTF.approx.M0 ) + qnorm( 1 - alpha / 2 ) * sqrt( Eta.M0^(-2) * Cov.Mat.M0[1,1] ) )


	#Exact

	if( (Eta.M0/sqrt(2)/sEta.M0) <= 10^307 ){

		f1 <- function(x)  exp( x^2 / ( 2 * sEta.M0^2 ) )

		MTTF.exact.M0 <- sqrt(2)*W / sEta.M0 * dawson(Eta.M0/sqrt(2)/sEta.M0) #sqrt(2)*W / sEta.M0 * dawson(Eta.M0/sqrt(2)/sEta.M0) #W / sEta.M0^2 * exp( - Eta.M0^2 / ( 2 * sEta.M0^2 ) ) * integrate( f1, 0, Eta.M0 )$value

		diff_MTTF.M0 <- numeric( length = 4 )

			sEta2.M0 <- sEta.M0^2

			MTTF.M0.Eta <- function(eta) sqrt(2) * W / sEta.M0 * dawson( eta / sqrt(2) / sEta.M0 )

			MTTF.M0.sEta2 <- function(seta2) sqrt(2) * W / sqrt(seta2) * dawson( Eta.M0 / sqrt(2) / sqrt(seta2) )

			diff_MTTF.M0[1] <- ( MTTF.M0.Eta( Eta.M0 - 2 * h.val ) - 8 * MTTF.M0.Eta( Eta.M0 - h.val ) +

						8 * MTTF.M0.Eta( Eta.M0 + h.val ) - MTTF.M0.Eta( Eta.M0 + 2 * h.val ) ) / 

						( 12 * h.val )

			diff_MTTF.M0[2] <- ( MTTF.M0.sEta2( sEta2.M0 - 2 * h.val ) - 8 * MTTF.M0.sEta2( sEta2.M0 - h.val ) +

						8 * MTTF.M0.sEta2( sEta2.M0 + h.val ) - MTTF.M0.sEta2( sEta2.M0 + 2 * h.val ) ) / 

						( 12 * h.val )

		var.MTTF.exact.M0 <- t(diff_MTTF.M0) %*% Cov.Mat.M0 %*% diff_MTTF.M0

		OICI.lower.MTTF.exact.M0 <- MTTF.exact.M0 - qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.exact.M0 )

		OICI.upper.MTTF.exact.M0 <- MTTF.exact.M0 + qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.exact.M0 )

	

		var2.MTTF.exact <- t( diff_MTTF.M0 / MTTF.exact.M0 ) %*% Cov.Mat.M0 %*% ( diff_MTTF.M0 / MTTF.exact.M0 )

		OICI.ln.lower.MTTF.exact.M0 <- exp( log( MTTF.exact.M0 ) - qnorm( 1 - alpha / 2 ) * sqrt( var2.MTTF.exact ) )

		OICI.ln.upper.MTTF.exact.M0 <- exp( log( MTTF.exact.M0 ) + qnorm( 1 - alpha / 2 ) * sqrt( var2.MTTF.exact ) )

	}else{
		MTTF.exact.M0 <- MTTF.approx.M0
		OICI.lower.MTTF.exact.M0 <- OICI.lower.MTTF.approx.M0
		OICI.upper.MTTF.exact.M0 <- OICI.upper.MTTF.approx.M0
		OICI.ln.lower.MTTF.exact.M0 <- OICI.ln.lower.MTTF.approx.M0
		OICI.ln.upper.MTTF.exact.M0 <- OICI.ln.upper.MTTF.approx.M0
	}

	list( MTTF.approx.M0 = MTTF.approx.M0, 

		CI.MTTF.approx.M0 = c(OICI.lower.MTTF.approx.M0, OICI.upper.MTTF.approx.M0), 

		CI.ln.MTTF.approx.M0 = c(OICI.ln.lower.MTTF.approx.M0, OICI.ln.upper.MTTF.approx.M0),

		MTTF.exact.M0 = MTTF.exact.M0, 

		CI.MTTF.exact.M0 = c(OICI.lower.MTTF.exact.M0, OICI.upper.MTTF.exact.M0), 

		CI.ln.MTTF.exact.M0 = c(OICI.ln.lower.MTTF.exact.M0, OICI.ln.upper.MTTF.exact.M0) )

}


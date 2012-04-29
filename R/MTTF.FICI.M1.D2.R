MTTF.FICI.M1.D2 <-
function(alpha, W, Eta.M1, sEta.M1, Cov.Mat.M1, n, h.val, ...){

	if( h.val > min(Eta.M1/2, sEta.M1^2/2) ) h.val <- min(Eta.M1/3, sEta.M1^2/3)

	#approx.

	MTTF.approx.M1 <- W / Eta.M1

	var.MTTF.approx.M1 <- ( -W / Eta.M1^2 )^2 * Cov.Mat.M1[1,1]

	FICI.lower.MTTF.approx.M1 <- MTTF.approx.M1 - qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.M1 / n )

	FICI.upper.MTTF.approx.M1 <- MTTF.approx.M1 + qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.approx.M1 / n )


	FICI.ln.lower.MTTF.approx.M1 <- exp( log( MTTF.approx.M1 ) - qnorm( 1 - alpha / 2 ) * sqrt( Eta.M1^(-2) * Cov.Mat.M1[1,1] / n ) )

	FICI.ln.upper.MTTF.approx.M1 <- exp( log( MTTF.approx.M1 ) + qnorm( 1 - alpha / 2 ) * sqrt( Eta.M1^(-2) * Cov.Mat.M1[1,1] / n ) )


	#Exact

	if( (Eta.M1/sqrt(2)/sEta.M1) <= 10^307 ){

		f1 <- function(x)  exp( x^2 / ( 2 * sEta.M1^2 ) )

		MTTF.exact.M1 <- sqrt(2)*W / sEta.M1 * dawson(Eta.M1/sqrt(2)/sEta.M1) #W / sEta.M1^2 * exp( - Eta.M1^2 / ( 2 * sEta.M1^2 ) ) * integrate( f1, 0, Eta.M1 )$value

	
		diff_MTTF.M1 <- numeric( length = 3 )

			sEta2.M1 <- sEta.M1^2

			MTTF.M1.Eta <- function(eta) sqrt(2) * W / sEta.M1 * dawson( eta / sqrt(2) / sEta.M1 )

			MTTF.M1.sEta2 <- function(seta2) sqrt(2) * W / sqrt(seta2) * dawson( Eta.M1 / sqrt(2) / sqrt(seta2) )

			diff_MTTF.M1[1] <- ( MTTF.M1.Eta( Eta.M1 - 2 * h.val ) - 8 * MTTF.M1.Eta( Eta.M1 - h.val ) +

						8 * MTTF.M1.Eta( Eta.M1 + h.val ) - MTTF.M1.Eta( Eta.M1 + 2 * h.val ) ) / 

						( 12 * h.val )

			diff_MTTF.M1[2] <- ( MTTF.M1.sEta2( sEta2.M1 - 2 * h.val ) - 8 * MTTF.M1.sEta2( sEta2.M1 - h.val ) +

						8 * MTTF.M1.sEta2( sEta2.M1 + h.val ) - MTTF.M1.sEta2( sEta2.M1 + 2 * h.val ) ) / 

						( 12 * h.val )

		var.MTTF.exact.M1 <- t(diff_MTTF.M1) %*% Cov.Mat.M1 %*% diff_MTTF.M1

		FICI.lower.MTTF.exact.M1 <- MTTF.exact.M1 - qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.exact.M1 / n )

		FICI.upper.MTTF.exact.M1 <- MTTF.exact.M1 + qnorm( 1 - alpha / 2 ) * sqrt( var.MTTF.exact.M1 / n )


		var2.MTTF.exact <- t( diff_MTTF.M1 / MTTF.exact.M1 ) %*% Cov.Mat.M1 %*% ( diff_MTTF.M1 / MTTF.exact.M1 )

		FICI.ln.lower.MTTF.exact.M1 <- exp( log( MTTF.exact.M1 ) - qnorm( 1 - alpha / 2 ) * sqrt( var2.MTTF.exact / n ) )

		FICI.ln.upper.MTTF.exact.M1 <- exp( log( MTTF.exact.M1 ) + qnorm( 1 - alpha / 2 ) * sqrt( var2.MTTF.exact / n ) )

	}else{
		MTTF.exact.M1 <- MTTF.approx.M1
		FICI.lower.MTTF.exact.M1 <- FICI.lower.MTTF.approx.M1
		FICI.upper.MTTF.exact.M1 <- FICI.upper.MTTF.approx.M1
		FICI.ln.lower.MTTF.exact.M1 <- FICI.ln.lower.MTTF.approx.M1
		FICI.ln.upper.MTTF.exact.M1 <- FICI.ln.upper.MTTF.approx.M1
	}

	list( MTTF.approx.M1 = MTTF.approx.M1, 

		CI.MTTF.approx.M1 = c(FICI.lower.MTTF.approx.M1, FICI.upper.MTTF.approx.M1), 

		CI.ln.MTTF.approx.M1 = c(FICI.ln.lower.MTTF.approx.M1, FICI.ln.upper.MTTF.approx.M1),

		MTTF.exact.M1 = MTTF.exact.M1, 

		CI.MTTF.exact.M1 = c(FICI.lower.MTTF.exact.M1, FICI.upper.MTTF.exact.M1), 

		CI.ln.MTTF.exact.M1 = c(FICI.ln.lower.MTTF.exact.M1, FICI.ln.upper.MTTF.exact.M1) )

}


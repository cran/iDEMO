cdf.FICI.logit.M2.D <-
function( Par.M2, t, Cov.Matrix.FI.M2, alpha, W, n, ... ) {

	Eta.M2 <- Par.M2[1]

	sB.M2 <- Par.M2[2]

	cdf.FICI.logit.M2.cal <- function(t){

		X0 <- 2 * Eta.M2 * W / sB.M2^2
		X1 <- ( Eta.M2 * t - W ) / sqrt( sB.M2^2 * t )
		X2 <- - ( Eta.M2 * t + W ) / sqrt( sB.M2^2 * t )

		diff.cdf.M2 <- matrix( 0, 2, 1 )
		sB2.M2 <- sB.M2^2

		if( X2< -38){
			cdf.M2 <- pnorm(X1)
			Ft.M2 <- expression( pnorm( ( Eta.M2 * t - W ) / sqrt( sB2.M2 * t ) ) )
			diff.cdf.M2[1,1] <- eval(D(Ft.M2, "Eta.M2"))
			diff.cdf.M2[2,1] <- eval(D(Ft.M2, "sB2.M2"))
		}else{
			U <- dnorm(X2)
			V <- pnorm(X2)
			if(floor(X0/709)==0){
				U <- U * exp(X0)
				V <- V * exp(X0)
			}else{
				for( i in 1 : floor(X0/709) ){
					U <- U * exp(709)
					V <- V * exp(709)
				}
				U <- U * exp( X0 %% 709 )
				V <- V * exp( X0 %% 709 )
			}

      	      cdf.M2 <- pnorm(X1) + V

			Ft.M2 <- expression( pnorm( ( Eta.M2 * t - W ) / sqrt( sB2.M2 * t ) ) +

		            exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) * pnorm( - ( Eta.M2 * t + W ) / 

				sqrt( sB2.M2 * t ) ) * exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) *

				exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) * exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) *

				exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) * exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) *

				exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) * exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) )

			diff.cdf.M2 <- matrix( 0, 2, 1 )

			diff.cdf.M2[1,1] <- eval(D(Ft.M2, "Eta.M2"))

			diff.cdf.M2[2,1] <- eval(D(Ft.M2, "sB2.M2"))

		}


		var.cdf.FI.M2 <-   c( t(  diff.cdf.M2 )  %*% Cov.Matrix.FI.M2 %*% diff.cdf.M2 ) / n 

		W.FI.M2 <-  exp(   qnorm( 1 - alpha / 2 ) *  sqrt(  var.cdf.FI.M2 ) /  ( cdf.M2 *  ( 1 - cdf.M2 ) ) ) 
  
		FICI.loWer.cdf.M2 <-  cdf.M2 /  ( cdf.M2 +  ( 1 - cdf.M2 ) *  W.FI.M2  ) 

		FICI.upper.cdf.M2 <-  cdf.M2 /  ( cdf.M2 +  ( 1 - cdf.M2 ) /  W.FI.M2  ) 


		return( c( cdf.M2, FICI.loWer.cdf.M2, FICI.upper.cdf.M2 ) )

	}


	M <- matrix( apply( as.matrix(t), 1, cdf.FICI.logit.M2.cal ), 3, length(t), byrow=FALSE )

	cdf.M2.out <- M[1,]

	FICI.loWer.cdf.M2.out <- M[2,]

	FICI.upper.cdf.M2.out <- M[3,]

	list( cdf.M2 = cdf.M2.out, lcl.cdf.M2 = FICI.loWer.cdf.M2.out, ucl.cdf.M2 = FICI.upper.cdf.M2.out )

}


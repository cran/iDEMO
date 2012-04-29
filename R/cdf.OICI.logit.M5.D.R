cdf.OICI.logit.M5.D <-
function( Par.M5, t, Cov.Matrix.OI.M5, alpha, W, ... ) {

	Eta.M5 <- Par.M5[1]

	sB.M5 <- Par.M5[2]

	sE.M5 <- Par.M5[3]

	cdf.OICI.logit.M5.cal <- function(t){

		X0 <- 2 * Eta.M5 * W / sB.M5^2
		X1 <- ( Eta.M5 * t - W ) / sqrt( sB.M5^2 * t )
		X2 <- - ( Eta.M5 * t + W ) / sqrt( sB.M5^2 * t )

		diff.cdf.M5 <- matrix( 0, 3, 1 )
		sB2.M5 <- sB.M5^2

		if( X2< -38){
			cdf.M5 <- pnorm(X1)
			Ft.M5 <- expression( pnorm( ( Eta.M5 * t - W ) / sqrt( sB2.M5 * t ) ) )
			diff.cdf.M5[1,1] <- eval(D(Ft.M5, "Eta.M5"))
			diff.cdf.M5[2,1] <- eval(D(Ft.M5, "sB2.M5"))
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

      	      cdf.M5 <- pnorm(X1) + V

			Ft.M5 <- expression( pnorm( ( Eta.M5 * t - W ) / sqrt( sB2.M5 * t ) ) +

		            exp( 2 * Eta.M5 * W / sB2.M5 ) * pnorm( - ( Eta.M5 * t + W ) / 

				sqrt( sB2.M5 * t ) ) )

			diff.cdf.M5 <- matrix( 0, 3, 1 )

			diff.cdf.M5[1,1] <- eval(D(Ft.M5, "Eta.M5"))

			diff.cdf.M5[2,1] <- eval(D(Ft.M5, "sB2.M5"))

		}


		var.cdf.OI.M5 <-   c( t(  diff.cdf.M5 )  %*% Cov.Matrix.OI.M5 %*% diff.cdf.M5 ) 

		W.OI.M5 <-  exp(   qnorm( 1 - alpha / 2 ) *  sqrt(  var.cdf.OI.M5 ) /  ( cdf.M5 *  ( 1 - cdf.M5 ) ) ) 
  
		OICI.loWer.cdf.M5 <-  cdf.M5 /  ( cdf.M5 +  ( 1 - cdf.M5 ) *  W.OI.M5  ) 

		OICI.upper.cdf.M5 <-  cdf.M5 /  ( cdf.M5 +  ( 1 - cdf.M5 ) /  W.OI.M5  ) 


		return( c( cdf.M5, OICI.loWer.cdf.M5, OICI.upper.cdf.M5 ) )

	}


	M <- matrix( apply( as.matrix(t), 1, cdf.OICI.logit.M5.cal ), 3, length(t), byrow=FALSE )

	cdf.M5.out <- M[1,]

	OICI.loWer.cdf.M5.out <- M[2,]

	OICI.upper.cdf.M5.out <- M[3,]

	list( cdf.M5 = cdf.M5.out, lcl.cdf.M5 = OICI.loWer.cdf.M5.out, ucl.cdf.M5 = OICI.upper.cdf.M5.out )

}


cdf.OICI.logit.M3.D <-
function( Par.M3, t, Cov.Matrix.OI.M3, alpha, W, ... ) {

	Eta.M3 <- Par.M3[1]

	sEta.M3 <- Par.M3[2]

	sB.M3 <- Par.M3[3]

	cdf.OICI.logit.M3.cal <- function(t){

		X0 <- 2 * Eta.M3 * W / sB.M3^2 + 2 * sEta.M3^2 * W^2 / sB.M3^4
		sB2t_sEta2t2 <-  sB.M3^2 * t + sEta.M3^2 * t^2
		X1 <- ( Eta.M3 * t - W ) / sqrt( sB2t_sEta2t2 ) 
		X2 <- - (  2 * sEta.M3^2 * W * t + sB.M3^2 * (  Eta.M3 * t + W ) )  / 
			(  sB.M3^2 * sqrt(  sB2t_sEta2t2 ) ) 

		diff.cdf.M3 <-  matrix( 0, 3, 1 ) 
		sEta2.M3 <- sEta.M3^2
		sB2.M3 <- sB.M3^2

		if( X2< -38){
			cdf.M3 <- pnorm( (  Eta.M3 * t - W )  / sqrt(  sB.M3^2 * t + sEta.M3^2 * t^2 ) )

			Ft.M3 <- expression( pnorm( ( Eta.M3 * t - W ) / sqrt( sB2.M3 * t + sEta2.M3 * t^2 ) ) )

			diff.cdf.M3[1,1] <- eval(D(Ft.M3, "Eta.M3"))

			diff.cdf.M3[2,1] <- eval(D(Ft.M3, "sEta2.M3"))

			diff.cdf.M3[3,1] <- eval(D(Ft.M3, "sB2.M3"))

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

      	      cdf.M3 <- pnorm(X1) + V

			Ft.M3 <- expression( pnorm( ( Eta.M3 * t - W ) / sqrt( sB2.M3 * t + sEta2.M3 * t^2 ) ) + 

				exp( 2 * Eta.M3 * W / sB2.M3 ) * pnorm( -( 2 * sEta2.M3 * W * t + sB2.M3 * ( Eta.M3 * t + W ) ) /

				( sB2.M3 * sqrt( sB2.M3 * t + sEta2.M3 * t^2 ) ) ) * exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) *

		            exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) * exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) * 

		            exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) * exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) * 

		            exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) * exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) * 

		            exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) )

			diff.cdf.M3 <- matrix( 0, 3, 1 )

			diff.cdf.M3[1,1] <- eval(D(Ft.M3, "Eta.M3"))

			diff.cdf.M3[2,1] <- eval(D(Ft.M3, "sEta2.M3"))

			diff.cdf.M3[3,1] <- eval(D(Ft.M3, "sB2.M3"))

		}

		var.cdf.OI.M3 <-   c( t(  diff.cdf.M3 )  %*% Cov.Matrix.OI.M3 %*% diff.cdf.M3 ) 

		W.OI.M3 <-  exp(   qnorm( 1 - alpha / 2 ) *  sqrt(  var.cdf.OI.M3 ) /  ( cdf.M3 *  ( 1 - cdf.M3 ) ) ) 
  
		OICI.lower.cdf.M3 <-  cdf.M3 /  ( cdf.M3 +  ( 1 - cdf.M3 ) *  W.OI.M3  ) 

		OICI.upper.cdf.M3 <-  cdf.M3 /  ( cdf.M3 +  ( 1 - cdf.M3 ) /  W.OI.M3  ) 


		return( c( cdf.M3, OICI.lower.cdf.M3, OICI.upper.cdf.M3 ) )

	}


	M <- matrix( apply( as.matrix(t), 1, cdf.OICI.logit.M3.cal ), 3, length(t), byrow=FALSE )

	cdf.M3.out <- M[1,]

	OICI.lower.cdf.M3.out <- M[2,]

	OICI.upper.cdf.M3.out <- M[3,]

	list( cdf.M3 = cdf.M3.out, lcl.cdf.M3 = OICI.lower.cdf.M3.out, ucl.cdf.M3 = OICI.upper.cdf.M3.out )


}


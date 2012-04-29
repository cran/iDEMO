cdf.OICI.logit.M0.D <-
function( Par.M0, t, Cov.Matrix.OI.M0, alpha, W, ... ) {

	Eta.M0 <- Par.M0[1]

	sEta.M0 <- Par.M0[2]

	sB.M0 <- Par.M0[3]

	sE.M0 <- Par.M0[4]

	cdf.OICI.logit.M0.cal <- function(t){

		X0 <- 2 * Eta.M0 * W / sB.M0^2 + 2 * sEta.M0^2 * W^2 / sB.M0^4
		sB2t_sEta2t2 <-  sB.M0^2 * t + sEta.M0^2 * t^2
		X1 <- ( Eta.M0 * t - W ) / sqrt( sB2t_sEta2t2 ) 
		X2 <- - (  2 * sEta.M0^2 * W * t + sB.M0^2 * (  Eta.M0 * t + W ) )  / 
			(  sB.M0^2 * sqrt(  sB2t_sEta2t2 ) ) 

		diff.cdf.M0 <-  matrix( 0, 4, 1 ) 
		sEta2.M0 <- sEta.M0^2
		sB2.M0 <- sB.M0^2

		if( X2< -38){
			cdf.M0 <- pnorm( (  Eta.M0 * t - W )  / sqrt(  sB.M0^2 * t + sEta.M0^2 * t^2 ) )

			Ft.M0 <- expression( pnorm( ( Eta.M0 * t - W ) / sqrt( sB2.M0 * t + sEta2.M0 * t^2 ) ) )

			diff.cdf.M0[1,1] <- eval(D(Ft.M0, "Eta.M0"))

			diff.cdf.M0[2,1] <- eval(D(Ft.M0, "sEta2.M0"))

			diff.cdf.M0[3,1] <- eval(D(Ft.M0, "sB2.M0"))

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

      	      cdf.M0 <- pnorm(X1) + V

			Ft.M0 <- expression( pnorm( ( Eta.M0 * t - W ) / sqrt( sB2.M0 * t + sEta2.M0 * t^2 ) ) + 

				exp( 2 * Eta.M0 * W / sB2.M0 ) * pnorm( -( 2 * sEta2.M0 * W * t + sB2.M0 * ( Eta.M0 * t + W ) ) /

				( sB2.M0 * sqrt( sB2.M0 * t + sEta2.M0 * t^2 ) ) ) * exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) *

		            exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) * exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) * 

		            exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) * exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) * 

		            exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) * exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) * 

		            exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) )

			diff.cdf.M0[1,1] <- eval(D(Ft.M0, "Eta.M0"))

			diff.cdf.M0[2,1] <- eval(D(Ft.M0, "sEta2.M0"))

			diff.cdf.M0[3,1] <- eval(D(Ft.M0, "sB2.M0"))

		}

		var.cdf.OI.M0 <-   c( t(  diff.cdf.M0 )  %*% Cov.Matrix.OI.M0 %*% diff.cdf.M0 ) 

		W.OI.M0 <-  exp(   qnorm( 1 - alpha / 2 ) *  sqrt(  var.cdf.OI.M0 ) /  ( cdf.M0 *  ( 1 - cdf.M0 ) ) ) 
  
		OICI.lower.cdf.M0 <-  cdf.M0 /  ( cdf.M0 +  ( 1 - cdf.M0 ) *  W.OI.M0  ) 

		OICI.upper.cdf.M0 <-  cdf.M0 /  ( cdf.M0 +  ( 1 - cdf.M0 ) /  W.OI.M0  ) 


		return( c( cdf.M0, OICI.lower.cdf.M0, OICI.upper.cdf.M0 ) )

	}


	M <- matrix( apply( as.matrix(t), 1, cdf.OICI.logit.M0.cal ), 3, length(t), byrow=FALSE )

	cdf.M0.out <- M[1,]

	OICI.lower.cdf.M0.out <- M[2,]

	OICI.upper.cdf.M0.out <- M[3,]

	list( cdf.M0 = cdf.M0.out, lcl.cdf.M0 = OICI.lower.cdf.M0.out, ucl.cdf.M0 = OICI.upper.cdf.M0.out )


}


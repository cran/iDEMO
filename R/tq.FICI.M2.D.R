tq.FICI.M2.D <-
function( P2, alpha, W, Eta.M2, sB.M2, Cov.Mat.M2, n, range.t, ... ){

	tp.FICI.M2.cal <- function(p2){

		t.percM2 <- function(t){

			x0 <- 2 * Eta.M2 * W / sB.M2^2

			x2 <- - ( Eta.M2 * t + W ) / sqrt( sB.M2^2 * t )

			if( x2< -38){
	      	      pnorm( ( Eta.M2 * t - W ) / sqrt( sB.M2^2 * t ) ) - p2
			}else{
				V <- pnorm(x2)
				if(floor(x0/709)==0){
					V <- V * exp(x0)
				}else{
					for( i in 1 : floor(x0/709) ){
						V <- V * exp(709)
					}
					V <- V * exp( x0 %% 709 )
				}
	      	      pnorm( ( Eta.M2 * t - W ) / sqrt( sB.M2^2 * t ) ) + V - p2
			}
		}
  
		tp <- uniroot( t.percM2, range.t )$root

		sB2.M2 <- sB.M2^2

		X0 <- 2 * Eta.M2 * W / sB.M2^2

		X1 <- ( Eta.M2 * tp - W ) / sqrt( sB.M2^2 * tp )

		X2 <- - ( Eta.M2 * tp + W ) / sqrt( sB.M2^2 * tp )

		diff.cdf.M2 <- matrix( 0, 2, 1 )

		if( X2< -38){

			Ft.M2 <- expression( pnorm( ( Eta.M2 * tp - W ) / sqrt( sB2.M2 * tp ) ) )

			diff.cdf.M2[1,1] <- eval(D(Ft.M2, "Eta.M2"))

			diff.cdf.M2[2,1] <- eval(D(Ft.M2, "sB2.M2"))

		}else{

			Ft.M2 <- expression( pnorm( ( Eta.M2 * tp - W ) / sqrt( sB2.M2 * tp ) ) +

		            exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) * pnorm( - ( Eta.M2 * tp + W ) / 

				sqrt( sB2.M2 * tp ) ) * exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) *

				exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) * exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) *

				exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) * exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) *

				exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) * exp( 2 * Eta.M2 * W / sB2.M2 / 8 ) )

			diff.cdf.M2 <- matrix( 0, 2, 1 )

			diff.cdf.M2[1,1] <- eval(D(Ft.M2, "Eta.M2"))

			diff.cdf.M2[2,1] <- eval(D(Ft.M2, "sB2.M2"))

		}


		f_T.M2 <- function(t) sqrt( W^2 / ( 2 * pi * t^3 * sB.M2^2 ) ) *

	                    exp( - ( W - Eta.M2 * t )^2 / ( 2 * t * sB.M2^2 ) )


		var.tp.FI.M2 <- ( t(diff.cdf.M2) %*% Cov.Mat.M2 %*% diff.cdf.M2 ) / ( f_T.M2(tp)^2 ) / n

		CI.loWer.tp.M2 <- tp - qnorm( 1 - alpha / 2 ) * sqrt(var.tp.FI.M2)

		CI.upper.tp.M2 <- tp + qnorm( 1 - alpha / 2 ) * sqrt(var.tp.FI.M2)




		var.tp.FI.ln.M2 <- ( t(diff.cdf.M2) %*% Cov.Mat.M2 %*% diff.cdf.M2 ) / ( tp * f_T.M2(tp) )^2 / n

		CI.ln.loWer.tp.M2 <- exp( log(tp) - qnorm( 1 - alpha / 2 )*sqrt(var.tp.FI.ln.M2) )

		CI.ln.upper.tp.M2 <- exp( log(tp) + qnorm( 1 - alpha / 2 )*sqrt(var.tp.FI.ln.M2) )

		c( tp, CI.loWer.tp.M2, CI.upper.tp.M2, CI.ln.loWer.tp.M2, CI.ln.upper.tp.M2 )

	}
	matrix( apply( as.matrix(P2), 1, tp.FICI.M2.cal ), length(P2), 5, byrow=TRUE )

}


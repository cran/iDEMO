tq.OICI.M5.D <-
function( P5, alpha, W, Eta.M5, sB.M5, Cov.Mat.M5, range.t, ... ){

	tp.OICI.M5.cal <- function(p5){

		t.percM5 <- function(t){

			x0 <- 2 * Eta.M5 * W / sB.M5^2

			x2 <- - ( Eta.M5 * t + W ) / sqrt( sB.M5^2 * t )

			if( x2< -38){
	      	      pnorm( ( Eta.M5 * t - W ) / sqrt( sB.M5^2 * t ) ) - p5
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
	      	      pnorm( ( Eta.M5 * t - W ) / sqrt( sB.M5^2 * t ) ) + V - p5
			}
		}
  
		tp <- uniroot( t.percM5, range.t )$root

		sB2.M5 <- sB.M5^2

		X0 <- 2 * Eta.M5 * W / sB.M5^2

		X1 <- ( Eta.M5 * tp - W ) / sqrt( sB.M5^2 * tp )

		X2 <- - ( Eta.M5 * tp + W ) / sqrt( sB.M5^2 * tp )

		diff.cdf.M5 <- matrix( 0, 3, 1 )

		if( X2< -38){

			Ft.M5 <- expression( pnorm( ( Eta.M5 * tp - W ) / sqrt( sB2.M5 * tp ) ) )

			diff.cdf.M5[1,1] <- eval(D(Ft.M5, "Eta.M5"))

			diff.cdf.M5[2,1] <- eval(D(Ft.M5, "sB2.M5"))

		}else{

			Ft.M5 <- expression( pnorm( ( Eta.M5 * tp - W ) / sqrt( sB2.M5 * tp ) ) +

		            exp( 2 * Eta.M5 * W / sB2.M5 ) * pnorm( - ( Eta.M5 * tp + W ) / 

				sqrt( sB2.M5 * tp ) ) )

			diff.cdf.M5 <- matrix( 0, 3, 1 )

			diff.cdf.M5[1,1] <- eval(D(Ft.M5, "Eta.M5"))

			diff.cdf.M5[2,1] <- eval(D(Ft.M5, "sB2.M5"))

		}


		f_T.M5 <- function(t) sqrt( W^2 / ( 2 * pi * t^3 * sB.M5^2 ) ) *

	                    exp( - ( W - Eta.M5 * t )^2 / ( 2 * t * sB.M5^2 ) )


		var.tp.OI.M5 <- ( t(diff.cdf.M5) %*% Cov.Mat.M5 %*% diff.cdf.M5 ) / ( f_T.M5(tp)^2 )

		CI.loWer.tp.M5 <- tp - qnorm( 1 - alpha / 2 ) * sqrt(var.tp.OI.M5)

		CI.upper.tp.M5 <- tp + qnorm( 1 - alpha / 2 ) * sqrt(var.tp.OI.M5)




		var.tp.OI.ln.M5 <- ( t(diff.cdf.M5) %*% Cov.Mat.M5 %*% diff.cdf.M5 ) / ( tp * f_T.M5(tp) )^2

		CI.ln.loWer.tp.M5 <- exp( log(tp) - qnorm( 1 - alpha / 2 )*sqrt(var.tp.OI.ln.M5) )

		CI.ln.upper.tp.M5 <- exp( log(tp) + qnorm( 1 - alpha / 2 )*sqrt(var.tp.OI.ln.M5) )

		c( tp, CI.loWer.tp.M5, CI.upper.tp.M5, CI.ln.loWer.tp.M5, CI.ln.upper.tp.M5 )

	}
	matrix( apply( as.matrix(P5), 1, tp.OICI.M5.cal ), length(P5), 5, byrow=TRUE )

}


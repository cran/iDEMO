tq.FICI.M0.D <-
function( P0, alpha, W, Eta.M0, sEta.M0, sB.M0, Cov.Mat.M0, n, range.t, ... ){

	tp.FICI.M0.cal <- function(p0){

		t.percM0 <- function(t){

			sB2t_sEta2t2 <- sB.M0^2 * t + sEta.M0^2 * t^2

			x0 <- 2 * Eta.M0 * W / sB.M0^2 + 2 * sEta.M0^2 * W^2 / sB.M0^4

      	      x1 <- ( Eta.M0 * t - W ) / sqrt( sB2t_sEta2t2 )

			x12 <- 2 * Eta.M0 * W / sB.M0^2

	            x2 <- - ( 2 * sEta.M0^2 * W * t + sB.M0^2 * ( Eta.M0 * t + W ) ) /

				( sB.M0^2 * sqrt( sB2t_sEta2t2 ) )

            	x3 <- 2 * sEta.M0^2 * W^2 / sB.M0^4 

			if( x2< -38){
	      	      pnorm(x1) - p0
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
	      	      pnorm(x1) + V - p0
			}

		}
  
		tp <- uniroot( t.percM0, range.t )$root


		diff.cdf.M0 <- matrix( 0, 4, 1 )

		sB2t_sEta2t2 <- sB.M0^2 * tp + sEta.M0^2 * tp^2

		X1 <- ( Eta.M0 * tp - W ) / sqrt( sB2t_sEta2t2 )

		X2 <- -( 2 * sEta.M0^2 * W * tp + sB.M0^2 * ( Eta.M0 * tp + W ) ) /

			( sB.M0^2 * sqrt( sB2t_sEta2t2 ) )

		X3 <- exp( ( sEta.M0^2 * W^2 / sB.M0^4 ) / 4 )

		X0 <- 2 * Eta.M0 * W / sB.M0^2 + 2 * sEta.M0^2 * W^2 / sB.M0^4

		sEta2.M0 <- sEta.M0^2

		sB2.M0 <- sB.M0^2

		if( X2< -38){

			Ft.M0 <- expression( pnorm( ( Eta.M0 * tp - W ) / sqrt( sB2.M0 * tp + sEta2.M0 * tp^2 ) ) )

			diff.cdf.M0[1,1] <- eval(D(Ft.M0, "Eta.M0"))

			diff.cdf.M0[2,1] <- eval(D(Ft.M0, "sEta2.M0"))

			diff.cdf.M0[3,1] <- eval(D(Ft.M0, "sB2.M0"))

		}else{

			Ft.M0 <- expression( pnorm( ( Eta.M0 * tp - W ) / sqrt( sB2.M0 * tp + sEta2.M0 * tp^2 ) ) + 

				exp( 2 * Eta.M0 * W / sB2.M0 ) * pnorm( -( 2 * sEta2.M0 * W * tp + sB2.M0 * ( Eta.M0 * tp + W ) ) /

				( sB2.M0 * sqrt( sB2.M0 * tp + sEta2.M0 * tp^2 ) ) ) * exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) *

		            exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) * exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) * 

		            exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) * exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) * 

		            exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) * exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) * 

		            exp( sEta2.M0 * W^2 / ( 4 * sB2.M0^2 ) ) )

			diff.cdf.M0[1,1] <- eval(D(Ft.M0, "Eta.M0"))

			diff.cdf.M0[2,1] <- eval(D(Ft.M0, "sEta2.M0"))

			diff.cdf.M0[3,1] <- eval(D(Ft.M0, "sB2.M0"))

		}

		sEta2t_sB2 <- sEta.M0^2 * tp + sB.M0^2

		f_T.M0 <- function(t) sqrt( W^2 / ( 2 * pi * t^3 * sEta2t_sB2 ) ) *

				exp( - ( W - Eta.M0 * t )^2 / ( 2 * t * sEta2t_sB2 ) )


		var.tp.FI.M0 <- ( t(diff.cdf.M0) %*% Cov.Mat.M0 %*% diff.cdf.M0 ) / ( f_T.M0(tp)^2 ) / n

		CI.lower.tp.M0 <- tp - qnorm( 1 - alpha / 2 ) * sqrt(var.tp.FI.M0)

		CI.upper.tp.M0 <- tp + qnorm( 1 - alpha / 2 ) * sqrt(var.tp.FI.M0)




		var.tp.FI.ln.M0 <- ( t(diff.cdf.M0) %*% Cov.Mat.M0 %*% diff.cdf.M0 ) / ( tp * f_T.M0(tp) )^2 / n

		CI.ln.lower.tp.M0 <- exp( log(tp) - qnorm( 1 - alpha / 2 ) * sqrt(var.tp.FI.ln.M0) )

		CI.ln.upper.tp.M0 <- exp( log(tp) + qnorm( 1 - alpha / 2 ) * sqrt(var.tp.FI.ln.M0) )

		c( tp, CI.lower.tp.M0, CI.upper.tp.M0, CI.ln.lower.tp.M0, CI.ln.upper.tp.M0, diff.cdf.M0 )

	}

	matrix( apply( as.matrix(P0), 1, tp.FICI.M0.cal ), length(P0), 9, byrow=TRUE )

}


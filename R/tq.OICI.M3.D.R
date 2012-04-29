tq.OICI.M3.D <-
function( P3, alpha, W, Eta.M3, sEta.M3, sB.M3, Cov.Mat.M3, range.t, ... ){

	tp.OICI.M3.cal <- function(p3){

		t.percM3 <- function(t){

			sB2t_sEta2t2 <- sB.M3^2 * t + sEta.M3^2 * t^2

			x0 <- 2 * Eta.M3 * W / sB.M3^2 + 2 * sEta.M3^2 * W^2 / sB.M3^4

      	      x1 <- ( Eta.M3 * t - W ) / sqrt( sB2t_sEta2t2 )

			x12 <- 2 * Eta.M3 * W / sB.M3^2

	            x2 <- - ( 2 * sEta.M3^2 * W * t + sB.M3^2 * ( Eta.M3 * t + W ) ) /

				( sB.M3^2 * sqrt( sB2t_sEta2t2 ) )

            	x3 <- 2 * sEta.M3^2 * W^2 / sB.M3^4 

			if( x2< -38){
	      	      pnorm(x1) - p3
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
	      	      pnorm(x1) + V - p3
			}

		}
  
		tp <- uniroot( t.percM3, range.t )$root


		diff.cdf.M3 <- matrix( 0, 3, 1 )

		sB2t_sEta2t2 <- sB.M3^2 * tp + sEta.M3^2 * tp^2

		X1 <- ( Eta.M3 * tp - W ) / sqrt( sB2t_sEta2t2 )

		X2 <- -( 2 * sEta.M3^2 * W * tp + sB.M3^2 * ( Eta.M3 * tp + W ) ) /

			( sB.M3^2 * sqrt( sB2t_sEta2t2 ) )

		X3 <- exp( ( sEta.M3^2 * W^2 / sB.M3^4 ) / 4 )

		X0 <- 2 * Eta.M3 * W / sB.M3^2 + 2 * sEta.M3^2 * W^2 / sB.M3^4

		sEta2.M3 <- sEta.M3^2

		sB2.M3 <- sB.M3^2

		if( X2< -38){

			Ft.M3 <- expression( pnorm( ( Eta.M3 * tp - W ) / sqrt( sB2.M3 * tp + sEta2.M3 * tp^2 ) ) )

			diff.cdf.M3[1,1] <- eval(D(Ft.M3, "Eta.M3"))

			diff.cdf.M3[2,1] <- eval(D(Ft.M3, "sEta2.M3"))

			diff.cdf.M3[3,1] <- eval(D(Ft.M3, "sB2.M3"))

		}else{

			Ft.M3 <- expression( pnorm( ( Eta.M3 * tp - W ) / sqrt( sB2.M3 * tp + sEta2.M3 * tp^2 ) ) + 

				exp( 2 * Eta.M3 * W / sB2.M3 ) * pnorm( -( 2 * sEta2.M3 * W * tp + sB2.M3 * ( Eta.M3 * tp + W ) ) /

				( sB2.M3 * sqrt( sB2.M3 * tp + sEta2.M3 * tp^2 ) ) ) * exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) *

		            exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) * exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) * 

		            exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) * exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) * 

		            exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) * exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) * 

		            exp( sEta2.M3 * W^2 / ( 4 * sB2.M3^2 ) ) )

			diff.cdf.M3 <- matrix( 0, 3, 1 )

			diff.cdf.M3[1,1] <- eval(D(Ft.M3, "Eta.M3"))

			diff.cdf.M3[2,1] <- eval(D(Ft.M3, "sEta2.M3"))

			diff.cdf.M3[3,1] <- eval(D(Ft.M3, "sB2.M3"))

		}

		sEta2t_sB2 <- sEta.M3^2 * tp + sB.M3^2

		f_T.M3 <- function(t) sqrt( W^2 / ( 2 * pi * t^3 * sEta2t_sB2 ) ) *

				exp( - ( W - Eta.M3 * t )^2 / ( 2 * t * sEta2t_sB2 ) )


		var.tp.OI.M3 <- ( t(diff.cdf.M3) %*% Cov.Mat.M3 %*% diff.cdf.M3 ) / ( f_T.M3(tp)^2 )

		CI.lower.tp.M3 <- tp - qnorm( 1 - alpha / 2 ) * sqrt(var.tp.OI.M3)

		CI.upper.tp.M3 <- tp + qnorm( 1 - alpha / 2 ) * sqrt(var.tp.OI.M3)


		var.tp.OI.ln.M3 <- ( t(diff.cdf.M3) %*% Cov.Mat.M3 %*% diff.cdf.M3 ) / ( tp * f_T.M3(tp) )^2

		CI.ln.lower.tp.M3 <- exp( log(tp) - qnorm( 1 - alpha / 2 ) * sqrt(var.tp.OI.ln.M3) )

		CI.ln.upper.tp.M3 <- exp( log(tp) + qnorm( 1 - alpha / 2 ) * sqrt(var.tp.OI.ln.M3) )

		c( tp, CI.lower.tp.M3, CI.upper.tp.M3, CI.ln.lower.tp.M3, CI.ln.upper.tp.M3, diff.cdf.M3 )

	}

	matrix( apply( as.matrix(P3), 1, tp.OICI.M3.cal ), length(P3), 8, byrow=TRUE )

}


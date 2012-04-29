tq.FICI.M1.D <-
function( P1, alpha, W, Eta.M1, sEta.M1, Cov.Mat.M1, n, range.t, ... ){

	tp.FICI.M1.cal <- function(p1){

		t.percM1 <- function(t){

			pnorm( ( Eta.M1 * t - W ) / ( sEta.M1 * t ) ) - p1

		}
  
		tp <- uniroot( t.percM1, range.t )$root

		sEta2.M1 <- sEta.M1^2

		Ft.M1 <- expression( pnorm( ( Eta.M1 * tp - W ) / ( sqrt(sEta2.M1) * tp ) ) )

		diff.cdf.M1 <- matrix( 0, 3, 1 )

		diff.cdf.M1[1,1] <- eval(D(Ft.M1, "Eta.M1"))

		diff.cdf.M1[2,1] <- eval(D(Ft.M1, "sEta2.M1"))


		f_T.M1 <- function(t) sqrt( W^2 / ( 2 * pi * t^4 * sEta.M1^2 ) ) *

				exp( - ( W - Eta.M1 * t )^2 / ( 2 * t^2 * sEta.M1^2 ) )


		var.tp.FI.M1 <- ( t(diff.cdf.M1) %*% Cov.Mat.M1 %*% diff.cdf.M1 ) / ( f_T.M1(tp)^2 ) / n

		CI.lower.tp.M1 <- tp - qnorm( 1 - alpha / 2 ) * sqrt(var.tp.FI.M1)

		CI.upper.tp.M1 <- tp + qnorm( 1 - alpha / 2 ) * sqrt(var.tp.FI.M1)




		var.tp.FI.ln.M1 <- ( t(diff.cdf.M1) %*% Cov.Mat.M1 %*% diff.cdf.M1 ) / ( tp * f_T.M1(tp) )^2 / n

		CI.ln.lower.tp.M1 <- exp( log(tp) - qnorm( 1 - alpha / 2 )*sqrt(var.tp.FI.ln.M1) )

		CI.ln.upper.tp.M1 <- exp( log(tp) + qnorm( 1 - alpha / 2 )*sqrt(var.tp.FI.ln.M1) )

		c( tp, CI.lower.tp.M1, CI.upper.tp.M1, CI.ln.lower.tp.M1, CI.ln.upper.tp.M1 )

	}

	matrix( apply( as.matrix(P1), 1, tp.FICI.M1.cal ), length(P1), 5, byrow=TRUE )

}


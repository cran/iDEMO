cdf.FICI.logit.M1.D <-
function( Par.M1, t, Cov.Matrix.FI.M1, alpha, W, n, ... ) {

	Eta.M1 <- Par.M1[1]

	sEta.M1 <- Par.M1[2]

	sE.M1 <- Par.M1[3]

	cdf.FICI.logit.M1.cal <- function(t){

		cdf.M1 <- pnorm( ( Eta.M1 * t - W ) / ( sEta.M1 * t ) )

		sEta2.M1 <- sEta.M1^2

		Ft.M1 <- expression( pnorm( ( Eta.M1 * t - W ) / ( sqrt(sEta2.M1) * t ) ) )

		diff.cdf.M1 <- matrix( 0, 3, 1 )

		diff.cdf.M1[1,1] <- eval(D(Ft.M1, "Eta.M1"))

		diff.cdf.M1[2,1] <- eval(D(Ft.M1, "sEta2.M1"))


		var.cdf.FI.M1 <-   c( t(  diff.cdf.M1 )  %*% Cov.Matrix.FI.M1 %*% diff.cdf.M1 ) / n 

		W.FI.M1 <-  exp(   qnorm( 1 - alpha / 2 ) *  sqrt(  var.cdf.FI.M1 ) /  ( cdf.M1 *  ( 1 - cdf.M1 ) ) ) 
  
		FICI.lower.cdf.M1 <-  cdf.M1 /  ( cdf.M1 +  ( 1 - cdf.M1 ) *  W.FI.M1  ) 

		FICI.upper.cdf.M1 <-  cdf.M1 /  ( cdf.M1 +  ( 1 - cdf.M1 ) /  W.FI.M1  ) 


		return( c( cdf.M1, FICI.lower.cdf.M1, FICI.upper.cdf.M1 ) )

	}


	M <- matrix( apply( as.matrix(t), 1, cdf.FICI.logit.M1.cal ), 3, length(t), byrow=FALSE )

	cdf.M1.out <- M[1,]

	FICI.lower.cdf.M1.out <- M[2,]

	FICI.upper.cdf.M1.out <- M[3,]

	list( cdf.M1 = cdf.M1.out, lcl.cdf.M1 = FICI.lower.cdf.M1.out, ucl.cdf.M1 = FICI.upper.cdf.M1.out )

}


cdf.OICI.logit.M1.D <-
function( Par.M1, t, Cov.Matrix.OI.M1, alpha, W, ... ) {

	Eta.M1 <- Par.M1[1]

	sEta.M1 <- Par.M1[2]

	sE.M1 <- Par.M1[3]

	cdf.OICI.logit.M1.cal <- function(t){

		cdf.M1 <- pnorm( ( Eta.M1 * t - W ) / ( sEta.M1 * t ) )

		sEta2.M1 <- sEta.M1^2

		Ft.M1 <- expression( pnorm( ( Eta.M1 * t - W ) / ( sqrt(sEta2.M1) * t ) ) )

		diff.cdf.M1 <- matrix( 0, 3, 1 )

		diff.cdf.M1[1,1] <- eval(D(Ft.M1, "Eta.M1"))

		diff.cdf.M1[2,1] <- eval(D(Ft.M1, "sEta2.M1"))


		var.cdf.OI.M1 <-   c( t(  diff.cdf.M1 )  %*% Cov.Matrix.OI.M1 %*% diff.cdf.M1 ) 

		W.OI.M1 <-  exp(   qnorm( 1 - alpha / 2 ) *  sqrt(  var.cdf.OI.M1 ) /  ( cdf.M1 *  ( 1 - cdf.M1 ) ) ) 
  
		OICI.lower.cdf.M1 <-  cdf.M1 /  ( cdf.M1 +  ( 1 - cdf.M1 ) *  W.OI.M1  ) 

		OICI.upper.cdf.M1 <-  cdf.M1 /  ( cdf.M1 +  ( 1 - cdf.M1 ) /  W.OI.M1  ) 


		return( c( cdf.M1, OICI.lower.cdf.M1, OICI.upper.cdf.M1 ) )

	}


	M <- matrix( apply( as.matrix(t), 1, cdf.OICI.logit.M1.cal ), 3, length(t), byrow=FALSE )

	cdf.M1.out <- M[1,]

	OICI.lower.cdf.M1.out <- M[2,]

	OICI.upper.cdf.M1.out <- M[3,]

	list( cdf.M1 = cdf.M1.out, lcl.cdf.M1 = OICI.lower.cdf.M1.out, ucl.cdf.M1 = OICI.upper.cdf.M1.out )

}


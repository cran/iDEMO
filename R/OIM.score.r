OIM.score <- function(Data, Eta, sEta, sB, sE){

	Time <- Data[,1]

	if( sEta>=0 & sB>=0 & sE>=0 && all(Time>0) ){

		m <- length(Time)

		Q <- matrix(rep(Time, length(Time)), ncol = length(Time))

		Q[row(Q) > col(Q)] <- sort(Q[row(Q) < col(Q)])

		Qinv <- matrix( 0, m, m )
		a <- numeric( length = m )
		a[1] <- 1 / Time[1]
		for(i in 2:m) a[i] <- 1 / ( Time[i] - Time[i-1] )
		for(i in 2:(m-1)) Qinv[ i, c( i-1, i, i+1 ) ] <- c( -a[i], a[i]+a[i+1], -a[i+1] )
		Qinv[ 1, c( 1, 2 ) ] <- c( a[1]+a[2], -a[2] )
		Qinv[ m, c( m-1, m ) ] <- c( -a[m], a[m] )

		if(sEta!=0 && sB!=0 && sE!=0){

			O <- sB^2 * Q + sE^2 * diag( x=1, m, m )

			S <- sEta^2 * Time %*% t(Time) + O

			tSt <- c( t(Time) %*% solve(S) %*% Time )

			tOt <- c( t(Time) %*% solve(O) %*% Time )

			one_sEta2tOt <- c( 1 + sEta^2 * tOt )

			OttO <- solve(O) %*% Time %*% t(Time) %*% solve(O)

			OQ <- solve(O) %*% Q

			tOQOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Time )

			OQO <- solve(O) %*% Q %*% solve(O)

			OQOttO <- solve(O) %*% Q %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

			OttOQO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% Q %*% solve(O)

			tOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Time )

			OO <- solve(O) %*% solve(O)

			OOttO <- solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

			OttOO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O)


			uu <- function( vec.y, ... ){

				u <- matrix(0,4,1)

				tSy <- t(Time) %*% solve(S) %*% vec.y

				yEtat <- ( vec.y - Eta * Time )

				u[1,1] <- tSy - Eta * tSt

				u[2,1] <- - tOt / ( 2 * one_sEta2tOt ) +  t(yEtat) %*% OttO %*% yEtat / 

				        ( 2 * one_sEta2tOt^2 )

				u[3,1] <- - sum( diag(OQ) ) / 2  +  sEta^2 * tOQOt / ( 2 * one_sEta2tOt )  +

				        t(yEtat) %*% OQO %*% yEtat / 2  -  sEta^2 / ( 2 * one_sEta2tOt ) * 

				        t(yEtat) %*% ( OQOttO + OttOQO ) %*% yEtat  +

				        sEta^4 * tOQOt / ( 2 * one_sEta2tOt^2 ) * 

				        t(yEtat) %*% OttO %*% yEtat

				u[4,1] <- - sum( diag( solve(O) ) ) / 2  +  sEta^2 * tOOt / ( 2 * one_sEta2tOt )  +

				        t(yEtat) %*% OO %*% yEtat / 2  -  sEta^2 / ( 2 * one_sEta2tOt ) * 

				        t(yEtat) %*% ( OOttO + OttOO ) %*% yEtat  +

				        sEta^4 * tOOt / ( 2 * one_sEta2tOt^2 ) * 

				        t(yEtat) %*% OttO %*% yEtat

				u %*% t(u)
        
			}

			uu.cal <- apply( Data[,-1], 2, uu )

			OI.B <- matrix( apply( uu.cal, 1, sum ), 4, 4, byrow=TRUE )

			return(OI.B)

		}

		if(sEta!=0 && sB==0 && sE!=0){

			O <- sE^2 * diag( x=1, m, m )

			S <- sEta^2 * Time %*% t(Time) + O

			tOt <- t(Time) %*% solve(O) %*% Time 

			one_sEta2tOt <- c( 1 + sEta^2 * tOt )

			tSt <- c( t(Time) %*% solve(S) %*% Time )

			OttO <- solve(O) %*% Time %*% t(Time) %*% solve(O)			

			tOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Time )

			OO <- solve(O) %*% solve(O)

			trace.Oinv <- sum( diag( solve(O) ) )

			OOttO <- solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

			OttOO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O)

			uu <- function( vec.y, ... ){

				u <- matrix(0,3,1)

				tSy <- t(Time) %*% solve(S) %*% vec.y

				yEtat <- ( vec.y - Eta * Time )

				u[1,1] <- tSy - Eta * tSt

				u[2,1] <- - tOt / ( 2 * one_sEta2tOt ) +  t(yEtat) %*% OttO %*% yEtat / 

				        ( 2 * one_sEta2tOt^2 )

				u[3,1] <- - trace.Oinv / 2  +  sEta^2 * tOOt / ( 2 * one_sEta2tOt )  +

				        t(yEtat) %*% OO %*% yEtat / 2  -  sEta^2 / ( 2 * one_sEta2tOt ) * 

				        t(yEtat) %*% ( OOttO + OttOO ) %*% yEtat  +

				        sEta^4 * tOOt / ( 2 * one_sEta2tOt^2 ) * 

				        t(yEtat) %*% OttO %*% yEtat

				u %*% t(u)
        
			}

			uu.cal <- apply( Data[,-1], 2, uu )

			OI.B <- matrix( apply( uu.cal, 1, sum ), 3, 3, byrow=TRUE )

			return(OI.B)

		}

		if(sEta==0 && sB!=0 && sE==0){

			tSt <- t(Time) %*% ( Qinv / sB^2 ) %*% Time

			OQ <- diag(m) / sB^2

			OQO <- Qinv / sB^4

			uu <- function( vec.y, ... ){

				u <- matrix(0,2,1)

				tSy <- t(Time) %*% ( Qinv / sB^2 ) %*% vec.y

				yEtat <- ( vec.y - Eta * Time )

				u[1,1] <- tSy - Eta * tSt

				u[2,1] <- - sum( diag(OQ) ) / 2  +  t(yEtat) %*% OQO %*% yEtat / 2

				u %*% t(u)
				
			}

			uu.cal <- apply( Data[,-1], 2, uu )

			OI.B <- matrix( apply( uu.cal, 1, sum ), 2, 2, byrow=TRUE )

			return(OI.B)

		}

		if(sEta!=0 && sB!=0 && sE==0){

			O <- sB^2 * Q

			Oinv <- Qinv / sB^2

			S <- sEta^2 * Time %*% t(Time) + O

			tSt <- c( t(Time) %*% solve(S) %*% Time )
			
			tOt <- t(Time) %*% Oinv %*% Time

			one_sEta2tOt <- c( 1 + sEta^2 * tOt )

			OttO <- Oinv %*% Time %*% t(Time) %*% Oinv
			
			OQ <- Oinv %*% Q

			tOQOt <- c( t(Time) %*% Oinv %*% Q %*% Oinv %*% Time )

			OQO <- Oinv %*% Q %*% Oinv

			OQOttO <- Oinv %*% Q %*% Oinv %*% Time %*% t(Time) %*% Oinv

			OttOQO <- Oinv %*% Time %*% t(Time) %*% Oinv %*% Q %*% Oinv

			uu <- function( vec.y, ... ){

				u <- matrix(0,3,1)

				tSy <- t(Time) %*% solve(S) %*% vec.y

				yEtat <- ( vec.y - Eta * Time )

				u[1,1] <- tSy - Eta * tSt

				u[2,1] <- - tOt / ( 2 * one_sEta2tOt ) +  t(yEtat) %*% OttO %*% yEtat / 

								( 2 * one_sEta2tOt^2 )

				u[3,1] <- - sum( diag(OQ) ) / 2  +  sEta^2 * tOQOt / ( 2 * one_sEta2tOt )  +

								t(yEtat) %*% OQO %*% yEtat / 2  -  sEta^2 / ( 2 * one_sEta2tOt ) * 

								t(yEtat) %*% ( OQOttO + OttOQO ) %*% yEtat  +

								sEta^4 * tOQOt / ( 2 * one_sEta2tOt^2 ) * 

								t(yEtat) %*% OttO %*% yEtat

				u %*% t(u)
        
			}

			uu.cal <- apply( Data[,-1], 2, uu )

			OI.B <- matrix( apply( uu.cal, 1, sum ), 3, 3, byrow=TRUE )

			return(OI.B)

		}

		if(sEta==0 && sB==0 && sE!=0){

			tOt <- c( t(Time) %*% diag(m) %*% Time / sE^2 )

			uu <- function( vec.y, ... ){

				u <- matrix(0,2,1)

				tOy <- c( t(Time) %*% diag(m) %*% vec.y / sE^2 )

				yEtat <- ( vec.y - Eta * Time )

				u[1,1] <- tOy - Eta * tOt

				u[2,1] <- - m / ( 2 * sE^2 ) + c( t(yEtat) %*% diag(m) %*% yEtat /

					( 2 * sE^4 ) )

				u %*% t(u)
				
			}

			uu.cal <- apply( Data[,-1], 2, uu )

			OI.B <- matrix( apply( uu.cal, 1, sum ), 2, 2, byrow=TRUE )

			return(OI.B)

		}

		if(sEta==0 && sB!=0 && sE!=0){

			S <- O <- sB^2 * Q + sE^2 * diag( x=1, m, m )

			tOt <- c( t(Time) %*% solve(O) %*% Time )

			OQ <- solve(O) %*% Q

			OQO <- solve(O) %*% Q %*% solve(O)

			OO <- solve(O) %*% solve(O)


			uu <- function( vec.y, ... ){

				u <- matrix(0,3,1)

				tOy <- t(Time) %*% solve(O) %*% vec.y

				yEtat <- ( vec.y - Eta * Time )

				u[1,1] <- tOy - Eta * tOt

				u[2,1] <- - sum( diag(OQ) ) / 2 + t(yEtat) %*% OQO %*% yEtat / 2

				u[3,1] <- - sum( diag( solve(O) ) ) / 2 + t(yEtat) %*% OO %*% yEtat / 2

				u %*% t(u)
				
			}

			uu.cal <- apply( Data[,-1], 2, uu )

			OI.B <- matrix( apply( uu.cal, 1, sum ), 3, 3, byrow=TRUE )

			return(OI.B)

		}

	}else{
		cat("The parameter and measuring time should not be negative and should be a positive vector, respectively!","\n")
	}

}
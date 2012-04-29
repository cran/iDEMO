OIM.hessian <- function(Data, Eta, sEta, sB, sE){

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

			tOQOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Time )

			tOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Time )

			OttO <- solve(O) %*% Time %*% t(Time) %*% solve(O)

			OQOttO <- solve(O) %*% Q %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

			OttOQO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% Q %*% solve(O)

			OOttO <- solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

			OttOO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O)

			OQOQ <- solve(O) %*% Q %*% solve(O) %*% Q

			tOQOQOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Q %*% solve(O) %*% Time )

			OQOQO <- solve(O) %*% Q %*% solve(O) %*% Q %*% solve(O)

			OQOQOttO <- solve(O) %*% Q %*% solve(O) %*% Q %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

			OQOttOQO <- solve(O) %*% Q %*% solve(O) %*% Time %*% t(Time) %*% solve(O) %*% Q %*% solve(O)

			OttOQOQO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Q %*% solve(O)

			QOO <- Q %*% solve(O) %*% solve(O)

			tOQOOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% solve(O) %*% Time )

			tOOQOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Q %*% solve(O) %*% Time )

			OOQO <- solve(O) %*% solve(O) %*% Q %*% solve(O)

			OQOO <- solve(O) %*% Q %*% solve(O) %*% solve(O)

			OOQOttO <- solve(O) %*% solve(O) %*% Q %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

			OQOOttO <- solve(O) %*% Q %*% solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

			OQOttOO <- solve(O) %*% Q %*% solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O)

			OOttOQO <- solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O) %*% Q %*% solve(O)

			OttOOQO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O) %*% Q %*% solve(O)

			OttOQOO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% Q %*% solve(O) %*% solve(O)

			OO <- solve(O) %*% solve(O)

			tOOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% solve(O) %*% Time )

			OOO <- solve(O) %*% solve(O) %*% solve(O)

			OOOttO <- solve(O) %*% solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

			OOttOO <- solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O)

			OttOOO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O) %*% solve(O)

			A.cal <- function( vec.y, ... ){

				tOy <- c( t(Time) %*% solve(O) %*% vec.y )

				tOQOy <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% vec.y )

				tOOy <- c( t(Time) %*% solve(O) %*% solve(O) %*% vec.y )

				yEtat <- c( vec.y - Eta * Time )

				M.A.temp <- matrix( 0, 4, 4 )

				M.A.temp[1,1] <- - tSt

				M.A.temp[1,2] <- M.A.temp[2,1] <- (- tOt * tOy + Eta * tOt^2 ) / one_sEta2tOt^2

				M.A.temp[1,3] <- M.A.temp[3,1] <- - tOQOy + Eta * tOQOt - sEta^4 * tOQOt / one_sEta2tOt^2 * ( tOt * tOy - Eta * tOt^2 ) +

                                            sEta^2 / one_sEta2tOt * ( tOQOt * tOy + tOt * tOQOy - 2 * Eta * tOt * tOQOt )

				M.A.temp[1,4] <- M.A.temp[4,1] <- - tOOy + Eta * tOOt - sEta^4 * tOOt / one_sEta2tOt^2 * ( tOt * tOy - Eta * tOt^2 ) +

                                            sEta^2 / one_sEta2tOt * ( tOOt * tOy + tOt * tOOy - 2 * Eta * tOt * tOOt )

				M.A.temp[2,2] <- tOt^2 / ( 2 * one_sEta2tOt^2 ) - tOt / one_sEta2tOt^3 * t(yEtat) %*% OttO %*% yEtat

				M.A.temp[2,3] <- M.A.temp[3,2] <- tOQOt / ( 2 * one_sEta2tOt^2 ) - t(yEtat) %*% ( OQOttO + OttOQO ) %*% yEtat / ( 2 * one_sEta2tOt^2 ) +

                                            sEta^2 * tOQOt / one_sEta2tOt^3 * t(yEtat) %*% OttO %*% yEtat

				M.A.temp[2,4] <- M.A.temp[4,2] <- tOOt / ( 2 * one_sEta2tOt^2 ) - t(yEtat) %*% ( OOttO + OttOO ) %*% yEtat / ( 2 * one_sEta2tOt^2 ) +

                                            sEta^2 * tOOt / one_sEta2tOt^3 * t(yEtat) %*% OttO %*% yEtat

				M.A.temp[3,3] <- sum( diag( OQOQ ) ) / 2 - sEta^2 * tOQOQOt / one_sEta2tOt + sEta^4 * tOQOt^2 / ( 2 * one_sEta2tOt^2 ) -

                                            t(yEtat) %*% OQOQO %*% yEtat + sEta^2 / one_sEta2tOt * t(yEtat) %*% ( OQOQOttO + OQOttOQO + OttOQOQO ) %*% yEtat -

                                            sEta^4 * tOQOt / one_sEta2tOt^2 * t(yEtat) %*% ( OQOttO + OttOQO ) %*% yEtat -
    
                                            sEta^4 * tOQOQOt / one_sEta2tOt^2 * t(yEtat) %*% OttO %*% yEtat + sEta^6 * tOQOt^2 / one_sEta2tOt^3 *

                                            t(yEtat) %*% OttO %*% yEtat

				M.A.temp[3,4] <- M.A.temp[4,3] <- sum( diag( QOO ) ) / 2 + sEta^4 * tOOt * tOQOt / ( 2 * one_sEta2tOt^2 ) - sEta^2 / ( 2 * one_sEta2tOt ) *

                                            ( tOQOOt + tOOQOt ) - 0.5 * t(yEtat) %*% ( OOQO + OQOO ) %*% yEtat -

                                            sEta^4 * tOOt / ( 2 * one_sEta2tOt^2 ) * t(yEtat) %*% ( OQOttO + OttOQO ) %*% yEtat +

                                            sEta^2 / ( 2 *one_sEta2tOt ) * t(yEtat) %*% ( OOQOttO + OQOOttO + OQOttOO + OOttOQO + OttOOQO + OttOQOO ) %*% yEtat +

                                            ( sEta^6 * tOOt * tOQOt / one_sEta2tOt^3 - sEta^4 * ( tOOQOt + tOQOOt ) / ( 2 *one_sEta2tOt^2 ) ) *

                                            t(yEtat) %*% OttO %*% yEtat -

                                            sEta^4 * tOQOt / ( 2 * one_sEta2tOt^2 ) * t(yEtat) %*% ( OOttO + OttOO ) %*% yEtat

				M.A.temp[4,4] <- sum( diag( OO ) ) / 2  + sEta^4 * tOOt^2 / ( 2 * one_sEta2tOt^2 ) - sEta^2 * tOOOt / one_sEta2tOt -

                                            t(yEtat) %*% OOO %*% yEtat - sEta^4 * tOOt / one_sEta2tOt^2 * t(yEtat) %*% ( OOttO + OttOO ) %*% yEtat +   

                                            sEta^2 / one_sEta2tOt * t(yEtat) %*% ( OOOttO + OOttOO + OttOOO ) %*% yEtat +

                                            sEta^6 * tOOt^2 / one_sEta2tOt^3 * t(yEtat) %*% OttO %*% yEtat -

                                            sEta^4 * tOOOt / one_sEta2tOt^2 * t(yEtat) %*% OttO %*% yEtat
 
				M.A.temp

			}

			A.cal.2 <- apply( Data[,-1], 2, A.cal )

			OI.A <- (-1) * matrix( apply( A.cal.2, 1, sum ), 4, 4 )

			return(OI.A)

		}

		if(sEta!=0 && sB==0 && sE!=0){

			O <- sE^2 * diag( x=1, m, m )

			S <- sEta^2 * Time %*% t(Time) + O

			tSt <- c( t(Time) %*% solve(S) %*% Time )

			tSSt <- c( t(Time) %*% solve(S) %*% solve(S) %*% Time )

			SS <- solve(S) %*% solve(S)

			SttS <- solve(S) %*% Time %*% t(Time) %*% solve(S)

			SSttS <- solve(S) %*% solve(S) %*% Time %*% t(Time) %*% solve(S)

			SttSS <- solve(S) %*% Time %*% t(Time) %*% solve(S) %*% solve(S)

			SSS <- solve(S) %*% solve(S) %*% solve(S)
    
			A.cal <- function( vec.y, ... ){

				tSy <- c( t(Time) %*% solve(S) %*% vec.y )

				tSSy <- c( t(Time) %*% solve(S) %*% solve(S) %*% vec.y )

				yEtat <- c( vec.y - Eta * Time )

				M.A.temp <- matrix( 0, 3, 3 )

				M.A.temp[1,1] <- - tSt

				M.A.temp[1,2] <- M.A.temp[2,1] <- - tSt * tSy + Eta * tSt^2

				M.A.temp[1,3] <- M.A.temp[3,1] <- - tSSy + Eta * tSSt

				M.A.temp[2,2] <- tSt^2 / 2 - tSt * t(yEtat) %*% SttS %*% yEtat

				M.A.temp[2,3] <- M.A.temp[3,2] <- tSSt / 2 - t(yEtat) %*% ( SSttS + SttSS ) %*% yEtat / 2

				M.A.temp[3,3] <- sum( diag( SS ) ) / 2 - t(yEtat) %*% SSS %*% yEtat
 
				M.A.temp

			}

			A.cal.2 <- apply( Data[,-1], 2, A.cal )

			OI.A <- (-1) * matrix( apply( A.cal.2, 1, sum ), 3, 3 )

			return(OI.A)

		}
		if(sEta==0 && sB!=0 && sE==0){

			S <- sB^2 * Q

			tQt <- c( t(Time) %*% Qinv %*% Time )

			SttS <- solve(S) %*% Time %*% t(Time) %*% solve(S)

			SSttS <- solve(S) %*% solve(S) %*% Time %*% t(Time) %*% solve(S)

			SttSS <- solve(S) %*% Time %*% t(Time) %*% solve(S) %*% solve(S)

			SSS <- solve(S) %*% solve(S) %*% solve(S)
    
			A.cal <- function( vec.y, ... ){

    				tQy <- c( t(Time) %*% Qinv %*% vec.y )

				yEtat <- c( vec.y - Eta * Time )

    				M.A.temp <- matrix( 0, 2, 2 )

    				M.A.temp[1,1] <- - tQt / sB^2

    				M.A.temp[1,2] <- M.A.temp[2,1] <- ( - tQy + Eta * tQt ) / sB^4

    				M.A.temp[2,2] <- m / ( 2 * sB^4 ) - t(yEtat) %*% Qinv %*% yEtat / sB^6
 
    				M.A.temp

			}

			A.cal.2 <- apply( Data[,-1], 2, A.cal )

			OI.A <- (-1) * matrix( apply( A.cal.2, 1, sum ), 2, 2 )

			return(OI.A)

		}

		if(sEta!=0 && sB!=0 && sE==0){

			O <- sB^2 * Q

			S <- sEta^2 * Time %*% t(Time) + O

			tSt <- c( t(Time) %*% solve(S) %*% Time )

			tSQSt <- c( t(Time) %*% solve(S) %*% Q %*% solve(S) %*% Time )

			SttS <- solve(S) %*% Time %*% t(Time) %*% solve(S)

			SQSttS <- solve(S) %*% Q %*% solve(S) %*% Time %*% t(Time) %*% solve(S)

			SttSQS <- solve(S) %*% Time %*% t(Time) %*% solve(S) %*% Q %*% solve(S)

			SQSQ <- solve(S) %*% Q %*% solve(S) %*% Q

			SQSQS <- solve(S) %*% Q %*% solve(S) %*% Q %*% solve(S)
    
			A.cal <- function( vec.y, ... ){

				tSy <- c( t(Time) %*% solve(S) %*% vec.y )

				tSQSy <- c( t(Time) %*% solve(S) %*% Q %*% solve(S) %*% vec.y )

				yEtat <- c( vec.y - Eta * Time )

				M.A.temp <- matrix( 0, 3, 3 )

				M.A.temp[1,1] <- - tSt

				M.A.temp[1,2] <- M.A.temp[2,1] <- - tSt * tSy + Eta * tSt^2

				M.A.temp[1,3] <- M.A.temp[3,1] <- - tSQSy + Eta * tSQSt

				M.A.temp[2,2] <- tSt^2 / 2 - tSt * t(yEtat) %*% SttS %*% yEtat

				M.A.temp[2,3] <- M.A.temp[3,2] <- tSQSt / 2 - t(yEtat) %*% ( SQSttS + SttSQS ) %*% yEtat / 2

				M.A.temp[3,3] <- sum( diag( SQSQ ) ) / 2 - t(yEtat) %*% SQSQS %*% yEtat
 
				M.A.temp

			}

			A.cal.2 <- apply( Data[,-1], 2, A.cal )

			OI.A <- (-1) * matrix( apply( A.cal.2, 1, sum ), 3, 3 )

			return(OI.A)

		}

		if(sEta==0 && sB==0 && sE!=0){

			tOt <- c( t(Time) %*% diag(m) %*% Time / sE^2 )

			tOOt <- c( t(Time) %*% diag(m) %*% Time / sE^4 )

			OOO <- diag(m) / sE^6
			
			A.cal <- function( vec.y, ... ){

				tOOy <- c( t(Time) %*% diag(m) %*% vec.y / sE^4 )

				yEtat <- c( vec.y - Eta * Time )

				M.A.temp <- matrix( 0, 2, 2 )

				M.A.temp[1,1] <- - tOt

				M.A.temp[1,2] <- M.A.temp[2,1] <- - tOOy + Eta * tOOt

				M.A.temp[2,2] <- m / ( 2 * sE^4 ) - t(yEtat) %*% OOO %*% yEtat
 
				M.A.temp

			}

			A.cal.2 <- apply( Data[,-1], 2, A.cal )

			OI.A <- (-1) * matrix( apply( A.cal.2, 1, sum ), 2, 2 )

			return(OI.A)

		}

		if(sEta==0 && sB!=0 && sE!=0){

			S <- O <- sB^2 * Q + sE^2 * diag( x=1, m, m )

			tOt <- c( t(Time) %*% solve(O) %*% Time )

			tOQOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Time )

			tOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Time )

			OQOQ <- solve(O) %*% Q %*% solve(O) %*% Q

			OQOQO <- solve(O) %*% Q %*% solve(O) %*% Q %*% solve(O)

			QOO <- Q %*% solve(O) %*% solve(O)

			OOQO <- solve(O) %*% solve(O) %*% Q %*% solve(O)

			OQOO <- solve(O) %*% Q %*% solve(O) %*% solve(O)

			OO <- solve(O) %*% solve(O)

			OOO <- solve(O) %*% solve(O) %*% solve(O)

			A.cal <- function( vec.y, ... ){

				tOQOy <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% vec.y )

				tOOy <- c( t(Time) %*% solve(O) %*% solve(O) %*% vec.y )

				yEtat <- c( vec.y - Eta * Time )

				M.A.temp <- matrix( 0, 3, 3 )

				M.A.temp[1,1] <- - tOt

				M.A.temp[1,2] <- M.A.temp[2,1] <- - tOQOy + Eta * tOQOt

				M.A.temp[1,3] <- M.A.temp[3,1] <- - tOOy + Eta * tOOt

				M.A.temp[2,2] <- sum( diag( OQOQ ) ) / 2 - t(yEtat) %*% OQOQO %*% yEtat

				M.A.temp[2,3] <- M.A.temp[3,2] <- sum( diag( QOO ) ) / 2 - t(yEtat) %*% ( OOQO + OQOO ) %*% yEtat / 2

				M.A.temp[3,3] <- sum( diag( OO ) ) / 2 - t(yEtat) %*% OOO %*% yEtat
 
				M.A.temp

			}

			A.cal.2 <- apply( Data[,-1], 2, A.cal )

			OI.A <- (-1) * matrix( apply( A.cal.2, 1, sum ), 3, 3 )

			return(OI.A)

		}

	}else{
		cat("The parameter and measuring time should not be negative and should be a positive vector, respectively!","\n")
	}

}
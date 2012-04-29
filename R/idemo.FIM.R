idemo.FIM <-
function( sEta, sB, sE, Time){
	if( is.numeric(sE) & is.numeric(sEta) & is.numeric(sB) & is.numeric(Time)){
		if( sEta>=0 & sB>=0 & sE>=0 && all(Time>0) ){
			m <- length(Time)
			Q <- matrix(rep(Time, length(Time)), ncol = length(Time))
			Q[row(Q) > col(Q)] <- sort(Q[row(Q) < col(Q)])
			O <- sB^2 * Q + sE^2 * diag( x=1, m, m )
			S <- sEta^2 * Time %*% t(Time) + O
			tSt <- c( t(Time) %*% solve(S) %*% Time )
			tOt <- c( t(Time) %*% solve(O) %*% Time )
			one_sEta2tOt <- c( 1 + sEta^2 * tOt )
			tOQOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Time )
			tOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Time )
			QOQO <- Q %*% solve(O) %*% Q %*% solve(O)
			tOQOQOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Q %*% solve(O) %*% Time )
			QOO <- Q %*% solve(O) %*% solve(O)
			tOOQOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Q %*% solve(O) %*% Time )
			OO <- solve(O) %*% solve(O)
			tOOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% solve(O) %*% Time )

			if(sEta!=0 && sB!=0 && sE!=0){
				FI <- matrix( 0, 4, 4 )
				FI[1,1] <- tSt
				FI[2,2] <- (tOt)^2 / ( 2 * one_sEta2tOt^2 )
				FI[2,3] <- FI[3,2] <- tOQOt / ( 2 * one_sEta2tOt^2 )
				FI[2,4] <- FI[4,2] <- tOOt / ( 2 * one_sEta2tOt^2 )
				FI[3,3] <- sum( diag(QOQO) ) / 2 + sEta^4 * tOQOt^2 / ( 2 * one_sEta2tOt^2 ) - 
						sEta^2 * tOQOQOt / one_sEta2tOt
				FI[3,4] <- FI[4,3] <- sum( diag(QOO) ) / 2 + sEta^4 * tOQOt * tOOt / ( 2 * one_sEta2tOt^2 ) -
						sEta^2 * tOOQOt / one_sEta2tOt
				FI[4,4] <- sum( diag(OO) ) / 2 + sEta^4 * tOOt^2 / ( 2 * one_sEta2tOt^2 ) - 
						sEta^2 * tOOOt / one_sEta2tOt
				return(FI)
			}

			if(sEta!=0 && sB==0 && sE!=0){
				FI <- matrix( 0, 3, 3 )
				tSSt <- c( t(Time) %*% solve(S) %*% solve(S) %*% Time )
				SS <- solve(S) %*% solve(S)
				FI[1,1] <- tSt
				FI[2,2] <- tSt^2 / 2
				FI[2,3] <- FI[3,2] <- tSSt / 2
				FI[3,3] <- sum( diag(SS) ) / 2
				return(FI)
			}

			if(sEta==0 && sB!=0 && sE==0){
				Q.inv <- matrix( 0, m, m )
				a <- numeric( length = m )
				a[1] <- 1 / Time[1]
				for(i in 2:m) a[i] <- 1 / ( Time[i] - Time[i-1] )
				for(i in 2:(m-1)) Q.inv[ i, c( i-1, i, i+1 ) ] <- c( -a[i], a[i]+a[i+1], -a[i+1] )
				Q.inv[ 1, c( 1, 2 ) ] <- c( a[1]+a[2], -a[2] )
				Q.inv[ m, c( m-1, m ) ] <- c( -a[m], a[m] )
				tQt <- c( t(Time) %*% Q.inv %*% Time )
				FI <- matrix( 0, 2, 2 )
				FI[1,1] <- tQt / sB^2
				FI[2,2] <- m / ( 2 * sB^4 )
				return(FI)
			}

			if(sEta!=0 && sB!=0 && sE==0){
				tSQSt <- c( t(Time) %*% solve(S) %*% Q %*% solve(S) %*% Time )
				QSQS <- Q %*% solve(S) %*% Q %*% solve(S)
				FI <- matrix( 0, 3, 3 )
				FI[1,1] <- tSt
				FI[2,2] <- (tSt)^2 / 2
				FI[2,3] <- FI[3,2] <- tSQSt / 2
				FI[3,3] <- sum( diag(QSQS) ) / 2
				return(FI)
			}

			if(sEta==0 && sB==0 && sE!=0){
				tOt <- c( t(Time) %*% diag(m) %*% Time / sE^2 )
				FI <- matrix( 0, 2, 2 )
				FI[1,1] <- tOt
				FI[2,2] <- m / ( 2 * sE^4 )
				return(FI)
			}

			if(sEta==0 && sB!=0 && sE!=0){
				FI <- matrix( 0, 3, 3 )
				FI[1,1] <- tOt
				FI[2,2] <- sum( diag(QOQO) ) / 2
				FI[2,3] <- FI[3,2] <- sum( diag(QOO) ) / 2
				FI[3,3] <- sum( diag(OO) ) / 2
				return(FI)
			}
		}else{
			cat("The parameter and measuring time should not be negative and should be a positive vector, respectively!","\n")
		}
	}else{
		cat("The parameter and measuring time should not be negative and should be a positive vector, respectively!","\n")
	}
}


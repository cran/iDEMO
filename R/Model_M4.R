Model_M4 <-
function(Data, int.eta, int.se, W, q, alpha, FIM, OIMA, OIMB, OIMC, env.in){

	Data <- as.matrix(Data)

	time <- Data[,1]

	m <- length(time); n <- length(Data[1,-1])

	Eta.M4 <- sum( t(time) %*% Data[,-1] ) / ( n * t(time) %*% time )

	sE.M4 <- sqrt( sum( diag( t( Data[,-1] - Eta.M4 * time ) %*% ( Data[,-1] - Eta.M4 * time ) ) ) / ( n * m ) )

	S <- sE.M4^2 * diag(m)

	par.M4 <- c(Eta.M4, sE.M4)

	assign('par.M4', par.M4, env.in)

	loglik.M4 <- - n * m * log( 2 * pi ) / 2 - n * log( det(S) ) / 2 - 

            sum( diag( as.matrix( t( Data[,-1] - Eta.M4 * time ) ) %*% ( diag(m) / sE.M4^2 ) %*% as.matrix( Data[,-1] - Eta.M4 * time ) ) ) / 2

	if(FIM == 1){

		Cov.Matrix.FI.M4 <- solve( FI.M4( sE.M4, time, m ) )

		par.FI.M4 <- par.FICI.M4(alpha, par.M4, Cov.Matrix.FI.M4, n )

	}

	if(OIMA == 1){

		A <- OI.A.M4( Data, Eta.M4, sE.M4, time, m )

		Cov.Matrix.OI.M4.A <- solve(A)

		par.OI.M4.A <- par.OICI.M4( alpha, par.M4, Cov.Matrix.OI.M4.A )

	}

	if(OIMB == 1){

		B <- OI.B.M4( Data, Eta.M4, sE.M4, time, m )

		Cov.Matrix.OI.M4.B <- solve(B)

		par.OI.M4.B <- par.OICI.M4( alpha, par.M4, Cov.Matrix.OI.M4.B )

	}

	if(OIMC == 1){

		A <- OI.A.M4( Data, Eta.M4, sE.M4, time, m )

		B <- OI.B.M4( Data, Eta.M4, sE.M4, time, m )

		Cov.Matrix.OI.M4.C <- solve(A) %*% B %*% solve(A)

		par.OI.M4.C <- par.OICI.M4( alpha, par.M4, Cov.Matrix.OI.M4.C )

	}


	eta.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	se.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)

	ind <- 1

	if(FIM==1){
		eta.CI.M[ind,] <- t(c(par.FI.M4$CI.Eta.M4,par.FI.M4$CI.ln.Eta.M4))
		se.CI.M[ind,] <- t(c(par.FI.M4$CI.sE.M4,par.FI.M4$CI.ln.sE.M4))
		ind <- ind + 1
	}

	if(OIMA==1){
		eta.CI.M[ind,] <- t(c(par.OI.M4.A$CI.Eta.M4,par.OI.M4.A$CI.ln.Eta.M4))
		se.CI.M[ind,] <- t(c(par.OI.M4.A$CI.sE.M4,par.OI.M4.A$CI.ln.sE.M4))
		ind <- ind + 1
	}

	if(OIMB==1){
		eta.CI.M[ind,] <- t(c(par.OI.M4.B$CI.Eta.M4,par.OI.M4.B$CI.ln.Eta.M4))
		se.CI.M[ind,] <- t(c(par.OI.M4.B$CI.sE.M4,par.OI.M4.B$CI.ln.sE.M4))
		ind <- ind + 1
	}

	if(OIMC==1){
		eta.CI.M[ind,] <- t(c(par.OI.M4.C$CI.Eta.M4,par.OI.M4.C$CI.ln.Eta.M4))
		se.CI.M[ind,] <- t(c(par.OI.M4.C$CI.sE.M4,par.OI.M4.C$CI.ln.sE.M4))
		ind <- ind + 1
	}

	eta.CI.M <- data.frame(eta.CI.M)
	se.CI.M <- data.frame(se.CI.M)

	row.name <- c("FIM", "OIM(hessian matrix)", "OIM(score vector)",
			"OIM(robust estimator)")

	name.ind <- c(FIM,OIMA,OIMB,OIMC)
	row.names(eta.CI.M) <- row.name[which(name.ind==1)]
	row.names(se.CI.M) <- row.name[which(name.ind==1)]

	names(eta.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	names(se.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")


	cat(paste(rep("#",70),sep="",collapse=""),"\n")
	cat("      Model: Measurement error ( traditional regression model )    ","\n")
	cat(paste(rep("#",70),sep="",collapse=""),"\n\n")

	if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
		cat(paste(rep("-",60),sep="",collapse=""),"\n")
		cat(paste("      MLE and ", 100*(1-alpha),"% confidence interval of parameters",sep=""),"\n")
		cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
	}else{
		cat(paste(rep("-",60),sep="",collapse=""),"\n")
		cat("      MLE of parameters","\n")
		cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
	}

	cat("eta =",format(Eta.M4,digits=22),"\n")
	if(ind>1) print(eta.CI.M)
	cat("\n")

	cat("sigma_epsilon =",format(sE.M4,digits=22),"\n")
	if(ind>1) print(se.CI.M)
	cat("\n")

	cat("log-likelihood =",format(loglik.M4,digits=22),"\n")
	cat("\n")

}


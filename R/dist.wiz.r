dist.wiz <-
function(Data, int.eta, int.seta, int.sb, int.se, optim.alg='Nelder-Mead'){

	Data <- as.matrix(Data)

	time <- Data[,1]

	m <- length(time); n <- length(Data[1,-1])

	Q <- matrix(rep(time, length(time)), ncol = length(time))

	Q[row(Q) > col(Q)] <- sort(Q[row(Q) < col(Q)])


	##---M0---

	fit_M0 <- function(){

		llk_M0 <- function(par){

			Eta <- par[1]

			sEta <- par[2]

			sB <- par[3]

			sE <- par[4]

			Omega <- sB^2 * Q + sE^2 * diag(m)

			Sigma <- sEta^2 * time %*% t(time) + Omega

			(-1) * ( - n * m * log( 2 * pi ) / 2 - n * log( det(Sigma) ) / 2 - 

	            sum( diag( as.matrix( t( Data[,-1] - Eta * time ) ) %*% solve(Sigma) %*% as.matrix( Data[,-1] - Eta * time ) ) ) / 2 )

		}

		initial.value.M0 <- c(int.eta, int.seta,int.sb, int.se)

		iter.M0 <- optim( initial.value.M0, llk_M0, control = list(factr = 1e-15), method=optim.alg )

		for(i in 1:5){

			initial.value.M0 <- iter.M0$par

			iter.M0 <- optim( initial.value.M0, llk_M0, control = list(factr = 1e-15), method=optim.alg )

		}

		cnt_M0 <- 0

		while( iter.M0 $ convergence != 0 && cnt_M0 < 10){

			initial.value.M0 <- iter.M0 $ par

			iter.M0 <- optim( initial.value.M0, llk_M0, control = list(factr = 1e-15), method=optim.alg )

			cnt_M0 <- cnt_M0 + 1

			if(cnt_M0 >= 10) {

				Eta.M0 <- NA

				sEta.M0 <- NA

				sB.M0 <- NA

				sE.M0 <- NA

			}

		}

		if(iter.M0 $ convergence == 0){

			Eta.M0 <- iter.M0 $ par[1]

			sEta.M0 <- abs( iter.M0 $ par[2] )

			sB.M0 <- abs( iter.M0 $ par[3] )

			sE.M0 <- abs( iter.M0 $ par[4] )

		}

		return(c(Eta.M0,sEta.M0,sB.M0,sE.M0))

	}

	par.M0 <- try(fit_M0(),T)

	Eta.M0 <- as.numeric(par.M0[1])

	sEta.M0 <- as.numeric(par.M0[2])

	sB.M0 <- as.numeric(par.M0[3])

	sE.M0 <- as.numeric(par.M0[4])

	if( all(!is.na(par.M0)) ){

		O <- sB.M0^2 * Q + sE.M0^2 * diag(m)

		S <- sEta.M0^2 * time %*% t(time) + O

		loglik.M0 <- - n * m * log( 2 * pi ) / 2 - n * log( det(S) ) / 2 - 

	            sum( diag( as.matrix( t( Data[,-1] - Eta.M0 * time ) ) %*% solve(S) %*% as.matrix( Data[,-1] - Eta.M0 * time ) ) ) / 2

		AIC.M0 <- model.selection( loglik.M0, length(par.M0), n )[1]

		BIC.M0 <- model.selection( loglik.M0, length(par.M0), n )[2]

		HQC.M0 <- model.selection( loglik.M0, length(par.M0), n )[3]

	}else{

		loglik.M0 <- NA

		AIC.M0 <- NA

		BIC.M0 <- NA

		HQC.M0 <- NA

	}

	##---M1---

	fit_M1 <- function(){

		llk_M1 <- function(par){

			Eta <- par[1]

			sEta <- par[2]

			sE <- par[3]

			Omega <- sE^2 * diag(m)

			Sigma <- sEta^2 * time %*% t(time) + Omega

			(-1) * ( - n * m * log( 2 * pi ) / 2 - n * log( det(Sigma) ) / 2 - 

			sum( diag( as.matrix( t( Data[,-1] - Eta * time ) ) %*% solve(Sigma) %*% as.matrix( Data[,-1] - Eta * time ) ) ) / 2 )

		}

		initial.value.M1 <- c(int.eta, int.seta, int.se)

		iter.M1 <- optim( initial.value.M1, llk_M1, control = list(factr = 1e-15), method=optim.alg )

		for(i in 1:5){

			initial.value.M1 <- iter.M1$par

			iter.M1 <- optim( initial.value.M1, llk_M1, control = list(factr = 1e-15), method=optim.alg )

		}

		cnt_M1 <- 0

		while( iter.M1 $ convergence != 0 && cnt_M1 < 10){

			initial.value.M1 <- iter.M1 $ par

			iter.M1 <- optim( initial.value.M1, llk_M1, control = list(factr = 1e-15), method=optim.alg )

			cnt_M1 <- cnt_M1 + 1

			if(cnt_M1 >= 10) {

				Eta.M1 <- NA

				sEta.M1 <- NA

				sE.M1 <- NA

			}

		}

		if(iter.M1 $ convergence == 0){

			Eta.M1 <- iter.M1 $ par[1]

			sEta.M1 <- abs( iter.M1 $ par[2] )

			sE.M1 <- abs( iter.M1 $ par[3] )

		}

		return(c(Eta.M1,sEta.M1,sE.M1))

	}

	par.M1 <- try(fit_M1(),T)

	Eta.M1 <- as.numeric(par.M1[1])

	sEta.M1 <- as.numeric(par.M1[2])

	sE.M1 <- as.numeric(par.M1[3])

	if( all(!is.na(par.M1)) ){

		O <- sE.M1^2 * diag(m)

		S <- sEta.M1^2 * time %*% t(time) + O

		loglik.M1 <- - n * m * log( 2 * pi ) / 2 - n * log( det(S) ) / 2 - 

	            sum( diag( as.matrix( t( Data[,-1] - Eta.M1 * time ) ) %*% solve(S) %*% as.matrix( Data[,-1] - Eta.M1 * time ) ) ) / 2

		AIC.M1 <- model.selection( loglik.M1, length(par.M1), n )[1]

		BIC.M1 <- model.selection( loglik.M1, length(par.M1), n )[2]

		HQC.M1 <- model.selection( loglik.M1, length(par.M1), n )[3]

	}else{

		loglik.M1 <- NA

		AIC.M1 <- NA

		BIC.M1 <- NA

		HQC.M1 <- NA

	}

	##---M2---

  	Q.inv <- matrix( 0, m, m )

	a <- numeric( length = m )

	a[1] <- 1 / time[1]

	for(i in 2:m) a[i] <- 1 / ( time[i] - time[i-1] )

  	for(i in 2:(m-1)) Q.inv[ i, c( i-1, i, i+1 ) ] <- c( -a[i], a[i]+a[i+1], -a[i+1] )

	Q.inv[ 1, c( 1, 2 ) ] <- c( a[1]+a[2], -a[2] )

 	Q.inv[ m, c( m-1, m ) ] <- c( -a[m], a[m] )


	Eta.M2 <- sum( t(time) %*% Q.inv %*% Data[,-1] ) / ( n * t(time) %*% Q.inv %*% time )

	sB.M2 <- sqrt( sum( diag( t( Data[,-1] - Eta.M2 * time ) %*% Q.inv %*% ( Data[,-1] - Eta.M2 * time ) ) ) / ( n * m ) )


	loglik.M2 <- - n * m * log( 2 * pi ) / 2 - n * log( det( sB.M2^2 * Q ) ) / 2 - 

            sum( diag( as.matrix( t( Data[,-1] - Eta.M2 * time ) ) %*% ( Q.inv / sB.M2^2 ) %*% as.matrix( Data[,-1] - Eta.M2 * time ) ) ) / 2


	par.M2 <- c(Eta.M2, sB.M2)

	AIC.M2 <- model.selection( loglik.M2, length(par.M2), n )[1]

	BIC.M2 <- model.selection( loglik.M2, length(par.M2), n )[2]

	HQC.M2 <- model.selection( loglik.M2, length(par.M2), n )[3]

	##---M3---

	fit_M3 <- function(){

		llk_M3 <- function(par){

			Eta <- par[1]

			sEta <- par[2]

			sB <- par[3]

			Omega <- sB^2 * Q

			Sigma <- sEta^2 * time %*% t(time) + Omega

			(-1) * ( - n * m * log( 2 * pi ) / 2 - n * log( det(Sigma) ) / 2 - 

			sum( diag( as.matrix( t( Data[,-1] - Eta * time ) ) %*% solve(Sigma) %*% as.matrix( Data[,-1] - Eta * time ) ) ) / 2 )

		}

		initial.value.M3 <- c(int.eta, int.seta,int.sb)

		iter.M3 <- optim( initial.value.M3, llk_M3, control = list(factr = 1e-15), method=optim.alg )

		for(i in 1:5){

			initial.value.M3 <- iter.M3$par

			iter.M3 <- optim( initial.value.M3, llk_M3, control = list(factr = 1e-15), method=optim.alg )

		}

		cnt_M3 <- 0

		while( iter.M3 $ convergence != 0 && cnt_M3 < 10){

			initial.value.M3 <- iter.M3 $ par

			iter.M3 <- optim( initial.value.M3, llk_M3, control = list(factr = 1e-15), method=optim.alg )

			cnt_M3 <- cnt_M3 + 1

			if(cnt_M3 >= 10) {

				Eta.M3 <- NA

				sEta.M3 <- NA

				sB.M3 <- NA

			}

		}

		if(iter.M3 $ convergence == 0){

			Eta.M3 <- iter.M3 $ par[1]

			sEta.M3 <- abs( iter.M3 $ par[2] )

			sB.M3 <- abs( iter.M3 $ par[3] )

		}

		return(c(Eta.M3,sEta.M3,sB.M3))

	}
	par.M3 <- try(fit_M3(),T)

	Eta.M3 <- as.numeric(par.M3[1])

	sEta.M3 <- as.numeric(par.M3[2])

	sB.M3 <- as.numeric(par.M3[3])

	if( all(!is.na(par.M3)) ){

		O <- sB.M3^2 * Q

		O.inv <- Q.inv / sB.M3^2 

		S <- sEta.M3^2 * time %*% t(time) + O

		loglik.M3 <- - n * m * log( 2 * pi ) / 2 - n * log( det(S) ) / 2 - 

	            sum( diag( as.matrix( t( Data[,-1] - Eta.M3 * time ) ) %*% solve(S) %*% as.matrix( Data[,-1] - Eta.M3 * time ) ) ) / 2

		AIC.M3 <- model.selection( loglik.M3, length(par.M3), n )[1]

		BIC.M3 <- model.selection( loglik.M3, length(par.M3), n )[2]

		HQC.M3 <- model.selection( loglik.M3, length(par.M3), n )[3]

	}else{

		loglik.M3 <- NA

		AIC.M3 <- NA

		BIC.M3 <- NA

		HQC.M3 <- NA

	}

	#---M4---

	Eta.M4 <- sum( t(time) %*% Data[,-1] ) / ( n * t(time) %*% time )

	sE.M4 <- sqrt( sum( diag( t( Data[,-1] - Eta.M4 * time ) %*% ( Data[,-1] - Eta.M4 * time ) ) ) / ( n * m ) )

	S <- sE.M4^2 * diag(m)


	loglik.M4 <- - n * m * log( 2 * pi ) / 2 - n * log( det(S) ) / 2 - 

            sum( diag( as.matrix( t( Data[,-1] - Eta.M4 * time ) ) %*% ( diag(m) / sE.M4^2 ) %*% as.matrix( Data[,-1] - Eta.M4 * time ) ) ) / 2


	par.M4 <- c(Eta.M4, sE.M4)

	AIC.M4 <- model.selection( loglik.M4, length(par.M4), n )[1]

	BIC.M4 <- model.selection( loglik.M4, length(par.M4), n )[2]

	HQC.M4 <- model.selection( loglik.M4, length(par.M4), n )[3]


	#---M5---

	fit_M5 <- function(){

		llk_M5 <- function(par){

			Eta <- par[1]

			sB <- par[2]

			sE <- par[3]

			Sigma <- sB^2 * Q + sE^2 * diag(m)

			(-1) * ( - n * m * log( 2 * pi ) / 2 - n * log( det(Sigma) ) / 2 - 

			sum( diag( as.matrix( t( Data[,-1] - Eta * time ) ) %*% solve(Sigma) %*% as.matrix( Data[,-1] - Eta * time ) ) ) / 2 )

		}

		initial.value.M5 <- c(int.eta,int.sb, int.se)

		iter.M5 <- optim( initial.value.M5, llk_M5, control = list(factr = 1e-15), method=optim.alg )

		for(i in 1:5){

			initial.value.M5 <- iter.M5$par

			iter.M5 <- optim( initial.value.M5, llk_M5, control = list(factr = 1e-15), method=optim.alg )

		}

		cnt_M5 <- 0

		while( iter.M5 $ convergence != 0 && cnt_M5 < 10){

			initial.value.M5 <- iter.M5 $ par

			iter.M5 <- optim( initial.value.M5, llk_M5, control = list(factr = 1e-15), method=optim.alg )

			cnt_M5 <- cnt_M5 + 1

			if(cnt_M5 >= 10) {

				Eta.M5 <- NA

				sB.M5 <- NA

				sE.M5 <- NA

			}

		}

		if(iter.M5 $ convergence == 0){

			Eta.M5 <- iter.M5 $ par[1]

			sB.M5 <- abs( iter.M5 $ par[2] )

			sE.M5 <- abs( iter.M5 $ par[3] )

		}

		return(c(Eta.M5,sB.M5,sE.M5))

	}
	par.M5 <- try(fit_M5(),T)

	Eta.M5 <- as.numeric(par.M5[1])

	sB.M5 <- as.numeric(par.M5[2])

	sE.M5 <- as.numeric(par.M5[3])

	if( all(!is.na(par.M5)) ){

		S <- sB.M5^2 * Q + sE.M5^2 * diag(m)

		loglik.M5 <- - n * m * log( 2 * pi ) / 2 - n * log( det(S) ) / 2 - 

	            sum( diag( as.matrix( t( Data[,-1] - Eta.M5 * time ) ) %*% solve(S) %*% as.matrix( Data[,-1] - Eta.M5 * time ) ) ) / 2

		AIC.M5 <- model.selection( loglik.M5, length(par.M5), n )[1]

		BIC.M5 <- model.selection( loglik.M5, length(par.M5), n )[2]

		HQC.M5 <- model.selection( loglik.M5, length(par.M5), n )[3]

	}else{

		loglik.M5 <- NA

		AIC.M5 <- NA

		BIC.M5 <- NA

		HQC.M5 <- NA

	}

	c.eta  <- c( Eta.M0,  Eta.M1,  Eta.M2, Eta.M3,  Eta.M4, Eta.M5 )
	c.seta <- c( sEta.M0, sEta.M1, NA,     sEta.M3, NA,     NA     )
	c.sb   <- c( sB.M0,   NA,      sB.M2,  sB.M3,   NA,     sB.M5  )
	c.se   <- c( sE.M0,   sE.M1,   NA,     NA,      sE.M4,  sE.M5  )
	c.lik <- c( loglik.M0, loglik.M1, loglik.M2, loglik.M3, loglik.M4, loglik.M5 )

	c.aic_val <- c( AIC.M0, AIC.M1, AIC.M2, AIC.M3, AIC.M4, AIC.M5 )
	c.aic_rank <- numeric(6)
	c.aic_rank[which(is.na(c.aic_val))] <- NA
	c.aic_rank[which(!is.na(c.aic_val))] <- rank(c.aic_val[which(!is.na(c.aic_val))])

	c.bic_val <- c( BIC.M0, BIC.M1, BIC.M2, BIC.M3, BIC.M4, BIC.M5 )
	c.bic_rank <- numeric(6)
	c.bic_rank[which(is.na(c.bic_val))] <- NA
	c.bic_rank[which(!is.na(c.bic_val))] <- rank(c.bic_val[which(!is.na(c.bic_val))])

	c.hqc_val <- c( HQC.M0, HQC.M1, HQC.M2, HQC.M3, HQC.M4, HQC.M5 )
	c.hqc_rank <- numeric(6)
	c.hqc_rank[which(is.na(c.hqc_val))] <- NA
	c.hqc_rank[which(!is.na(c.hqc_val))] <- rank(c.hqc_val[which(!is.na(c.hqc_val))])

	c.aic <- paste(round(c.aic_val,digits=4),'(',c.aic_rank,')')
	c.bic <- paste(round(c.bic_val,digits=4),'(',c.bic_rank,')')
	c.hqc <- paste(round(c.hqc_val,digits=4),'(',c.hqc_rank,')')

	model.name <- c('M0','M1','M2','M3','M4','M5')

	T <- data.frame(c.eta, c.seta, c.sb, c.se, c.lik, c.aic, c.bic, c.hqc)
	names(T) <- c('eta', 'sigma_eta', 'sigma_B', 'sigma_epsilon', 'log-likelihood', 'AIC(rank)', 'BIC(rank)', 'HQC(rank)')
	print(T)
	cat('\n')
	if(sum(is.na(c.lik))>=1){
		cat('Please try different intial values and/or the optimization algorithm.','\n')
	}
}

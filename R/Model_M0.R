Model_M0 <-
function(Data, int.eta, int.seta, int.sb, int.se, W, q, alpha, FIM, OIMA, OIMB, OIMC, t2, h.val, pseudoVal, MTTFval, quanval, pdfplotval, cdfplotval, ppplotval, qqplotval, kstestval, adtestval, dataname, optim.alg, env.in ){

	Data <- as.matrix(Data)

	time <- Data[,1]

	m <- length(time); n <- length(Data[1,-1])

	Q <- matrix(rep(time, length(time)), ncol = length(time))

	Q[row(Q) > col(Q)] <- sort(Q[row(Q) < col(Q)])

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

			stop('Please try different initial values (and/or the optimization algorithm)\n  to obtain more stable and more accurate estimates of the unknown parameters.')

		}
	}

	if(iter.M0 $ convergence == 0){

		Eta.M0 <- iter.M0 $ par[1]

		sEta.M0 <- abs( iter.M0 $ par[2] )

		sB.M0 <- abs( iter.M0 $ par[3] )

		sE.M0 <- abs( iter.M0 $ par[4] )

	}

	O <- sB.M0^2 * Q + sE.M0^2 * diag(m)

	S <- sEta.M0^2 * time %*% t(time) + O

	par.M0 <- c(Eta.M0, sEta.M0, sB.M0, sE.M0)

	assign('par.M0', par.M0, env.in)

	loglik.M0 <- - n * m * log( 2 * pi ) / 2 - n * log( det(S) ) / 2 - 

	            sum( diag( as.matrix( t( Data[,-1] - Eta.M0 * time ) ) %*% solve(S) %*% as.matrix( Data[,-1] - Eta.M0 * time ) ) ) / 2

	pdf_M0 <- function(t){

			sqrt( W^2 / ( 2 * pi * t^3 * ( sEta.M0^2 * t + sB.M0^2 ) ) ) *

			exp( - ( W - Eta.M0 * t )^2 / ( 2 * t * ( sEta.M0^2 * t + sB.M0^2 ) ) )

			}

	cdf_M0 <- function(t){

		sB2t_sEta2t2 <- sB.M0^2 * t + sEta.M0^2 * t^2

		x0 <- 2 * Eta.M0 * W / sB.M0^2 + 2 * sEta.M0^2 * W^2 / sB.M0^4

     	      x1 <- ( Eta.M0 * t - W ) / sqrt( sB2t_sEta2t2 )

		x12 <- 2 * Eta.M0 * W / sB.M0^2

            x2 <- - ( 2 * sEta.M0^2 * W * t + sB.M0^2 * ( Eta.M0 * t + W ) ) /

			( sB.M0^2 * sqrt( sB2t_sEta2t2 ) )

           	x3 <- 2 * sEta.M0^2 * W^2 / sB.M0^4 

		if( x2< -38){
      	      pnorm(x1)
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
      	      pnorm(x1) + V
		}
	}

	cdf_M0_inf <- function(){

		x0 <- 2 * Eta.M0 * W / sB.M0^2 + 2 * sEta.M0^2 * W^2 / sB.M0^4

     	      x1 <- Eta.M0 / sEta.M0

            x2 <- - 2 * sEta.M0 * W / sB.M0^2 - Eta.M0 / sEta.M0

		if( x2< -38){
      	      pnorm(x1)
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
      	      pnorm(x1) + V
		}
	}

	cdf_M0_INF <- cdf_M0_inf()
  
	tq_M0 <- function(p0){

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
  
		uniroot( t.percM0, range(t2) )$root

	}

	if(FIM == 1){

		Cov.Matrix.FI.M0 <- solve( FI.M0( sEta.M0, Q, O, S, time ) )

		par.FI.M0 <- par.FICI.M0(alpha, par.M0, Cov.Matrix.FI.M0, n )

		if(MTTFval==1 && Eta.M0>0){
			MTTF.FI.M0 <- MTTF.FICI.M0.D2( alpha, W, Eta.M0, sEta.M0, Cov.Matrix.FI.M0, n, h.val )
		}
		if(quanval==1 && cdf_M0_INF>max(q)) tq.FI.M0 <- try(tq.FICI.M0.D( q, alpha, W, Eta.M0, sEta.M0, sB.M0, Cov.Matrix.FI.M0, n, range(t2) ),T)

		if(cdfplotval==1) CDF.FI.M0 <- cdf.FICI.logit.M0.D( Par.M0 = par.M0, t=t2, Cov.Matrix.FI.M0, alpha=alpha, W=W, n )

	}

	if(OIMA == 1){

		A <- OI.A.M0( Data, Eta.M0, sEta.M0, Q, O, S, time )

		Cov.Matrix.OI.M0.A <- solve(A)

		par.OI.M0.A <- par.OICI.M0( alpha, par.M0, Cov.Matrix.OI.M0.A )

		if(MTTFval==1 && Eta.M0>0){
			MTTF.OIA.M0 <- MTTF.OICI.M0.D2( alpha, W, Eta.M0, sEta.M0, Cov.Matrix.OI.M0.A, h.val )
		}
		if(quanval==1 && cdf_M0_INF>max(q)) tq.OIA.M0 <- try(tq.OICI.M0.D( q, alpha, W, Eta.M0, sEta.M0, sB.M0, Cov.Matrix.OI.M0.A, range(t2) ),T)

		if(cdfplotval==1) CDF.OIA.M0 <- cdf.OICI.logit.M0.D( Par.M0 = par.M0, t=t2, Cov.Matrix.OI.M0.A, alpha=alpha, W=W )

	}

	if(OIMB == 1){

		B <- OI.B.M0( Data, Eta.M0, sEta.M0, Q, O, S, time )

		Cov.Matrix.OI.M0.B <- solve(B)

		par.OI.M0.B <- par.OICI.M0( alpha, par.M0, Cov.Matrix.OI.M0.B )

		if(MTTFval==1 && Eta.M0>0){
			MTTF.OIB.M0 <- MTTF.OICI.M0.D2( alpha, W, Eta.M0, sEta.M0, Cov.Matrix.OI.M0.B, h.val )
		}
		if(quanval==1 && cdf_M0_INF>max(q)) tq.OIB.M0 <- try(tq.OICI.M0.D( q, alpha, W, Eta.M0, sEta.M0, sB.M0, Cov.Matrix.OI.M0.B, range(t2) ),T)

		if(cdfplotval==1) CDF.OIB.M0 <- cdf.OICI.logit.M0.D( Par.M0 = par.M0, t=t2, Cov.Matrix.OI.M0.B, alpha=alpha, W=W )

	}

	if(OIMC == 1){

		A <- OI.A.M0( Data, Eta.M0, sEta.M0, Q, O, S, time )

		B <- OI.B.M0( Data, Eta.M0, sEta.M0, Q, O, S, time )

		Cov.Matrix.OI.M0.C <- solve(A) %*% B %*% solve(A)

		par.OI.M0.C <- par.OICI.M0( alpha, par.M0, Cov.Matrix.OI.M0.C )

		if(MTTFval==1 && Eta.M0>0){
			MTTF.OIC.M0 <- MTTF.OICI.M0.D2( alpha, W, Eta.M0, sEta.M0, Cov.Matrix.OI.M0.C, h.val )
		}
		if(quanval==1 && cdf_M0_INF>max(q)) tq.OIC.M0 <- try(tq.OICI.M0.D( q, alpha, W, Eta.M0, sEta.M0, sB.M0, Cov.Matrix.OI.M0.C, range(t2) ),T)

		if(cdfplotval==1) CDF.OIC.M0 <- cdf.OICI.logit.M0.D( Par.M0 = par.M0, t=t2, Cov.Matrix.OI.M0.C, alpha=alpha, W=W )

	}

	if(MTTFval==1 && Eta.M0>0){
		MTTF.approx.M0 <- MTTF.M0( par.M0[1], par.M0[2], W )[1]
		MTTF.exact.M0 <- MTTF.M0( par.M0[1], par.M0[2], W )[2]
	}

	if(quanval==1 && cdf_M0_INF>max(q)){

		fun_tp <- function(p0){

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
  
			tp <- uniroot( t.percM0, range(t2) )$root

		}

		tt <- try(apply( t(q), 2, fun_tp),T)

	}

	eta.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	seta.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	sb.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	se.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	MTTF.approx.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	MTTF.exact.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	tq.CI.M <- array(0,c(sum(FIM,OIMA,OIMB,OIMC),4,length(q)))

	ind <- 1

	if(FIM==1){
		eta.CI.M[ind,] <- t(c(par.FI.M0$CI.Eta.M0,par.FI.M0$CI.ln.Eta.M0))
		seta.CI.M[ind,] <- t(c(par.FI.M0$CI.sEta.M0,par.FI.M0$CI.ln.sEta.M0))
		sb.CI.M[ind,] <- t(c(par.FI.M0$CI.sB.M0,par.FI.M0$CI.ln.sB.M0))
		se.CI.M[ind,] <- t(c(par.FI.M0$CI.sE.M0,par.FI.M0$CI.ln.sE.M0))
		if(MTTFval==1 && Eta.M0>0){
			MTTF.approx.CI.M[ind,] <- t(c(MTTF.FI.M0$CI.MTTF.approx.M0, MTTF.FI.M0$CI.ln.MTTF.approx.M0))
			MTTF.exact.CI.M[ind,] <- t(c(MTTF.FI.M0$CI.MTTF.exact.M0, MTTF.FI.M0$CI.ln.MTTF.exact.M0))
		}
		if(quanval==1 && cdf_M0_INF>max(q) && length(tq.FI.M0)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.FI.M0[i,2:5])
			}
		}
		ind <- ind + 1
	}

	if(OIMA==1){
		eta.CI.M[ind,] <- t(c(par.OI.M0.A$CI.Eta.M0,par.OI.M0.A$CI.ln.Eta.M0))
		seta.CI.M[ind,] <- t(c(par.OI.M0.A$CI.sEta.M0,par.OI.M0.A$CI.ln.sEta.M0))
		sb.CI.M[ind,] <- t(c(par.OI.M0.A$CI.sB.M0,par.OI.M0.A$CI.ln.sB.M0))
		se.CI.M[ind,] <- t(c(par.OI.M0.A$CI.sE.M0,par.OI.M0.A$CI.ln.sE.M0))
		if(MTTFval==1 && Eta.M0>0){
			MTTF.approx.CI.M[ind,] <- t(c(MTTF.OIA.M0$CI.MTTF.approx.M0, MTTF.OIA.M0$CI.ln.MTTF.approx.M0))
			MTTF.exact.CI.M[ind,] <- t(c(MTTF.OIA.M0$CI.MTTF.exact.M0, MTTF.OIA.M0$CI.ln.MTTF.exact.M0))
		}
		if(quanval==1 && cdf_M0_INF>max(q) && length(tq.OIA.M0)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.OIA.M0[i,2:5])
			}
		}
		ind <- ind + 1
	}

	if(OIMB==1){
		eta.CI.M[ind,] <- t(c(par.OI.M0.B$CI.Eta.M0,par.OI.M0.B$CI.ln.Eta.M0))
		seta.CI.M[ind,] <- t(c(par.OI.M0.B$CI.sEta.M0,par.OI.M0.B$CI.ln.sEta.M0))
		sb.CI.M[ind,] <- t(c(par.OI.M0.B$CI.sB.M0,par.OI.M0.B$CI.ln.sB.M0))
		se.CI.M[ind,] <- t(c(par.OI.M0.B$CI.sE.M0,par.OI.M0.B$CI.ln.sE.M0))
		if(MTTFval==1 && Eta.M0>0){
			MTTF.approx.CI.M[ind,] <- t(c(MTTF.OIB.M0$CI.MTTF.approx.M0, MTTF.OIB.M0$CI.ln.MTTF.approx.M0))
			MTTF.exact.CI.M[ind,] <- t(c(MTTF.OIB.M0$CI.MTTF.exact.M0, MTTF.OIB.M0$CI.ln.MTTF.exact.M0))
		}
		if(quanval==1 && cdf_M0_INF>max(q) && length(tq.OIB.M0)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.OIB.M0[i,2:5])
			}
		}
		ind <- ind + 1
	}

	if(OIMC==1){
		eta.CI.M[ind,] <- t(c(par.OI.M0.C$CI.Eta.M0,par.OI.M0.C$CI.ln.Eta.M0))
		seta.CI.M[ind,] <- t(c(par.OI.M0.C$CI.sEta.M0,par.OI.M0.C$CI.ln.sEta.M0))
		sb.CI.M[ind,] <- t(c(par.OI.M0.C$CI.sB.M0,par.OI.M0.C$CI.ln.sB.M0))
		se.CI.M[ind,] <- t(c(par.OI.M0.C$CI.sE.M0,par.OI.M0.C$CI.ln.sE.M0))
		if(MTTFval==1 && Eta.M0>0){
			MTTF.approx.CI.M[ind,] <- t(c(MTTF.OIC.M0$CI.MTTF.approx.M0, MTTF.OIC.M0$CI.ln.MTTF.approx.M0))
			MTTF.exact.CI.M[ind,] <- t(c(MTTF.OIC.M0$CI.MTTF.exact.M0, MTTF.OIC.M0$CI.ln.MTTF.exact.M0))
		}
		if(quanval==1 && cdf_M0_INF>max(q) && length(tq.OIC.M0)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.OIC.M0[i,2:5])
			}
		}
		ind <- ind + 1
	}

	eta.CI.M <- data.frame(eta.CI.M)
	seta.CI.M <- data.frame(seta.CI.M)
	sb.CI.M <- data.frame(sb.CI.M)
	se.CI.M <- data.frame(se.CI.M)
	if(MTTFval==1 && Eta.M0>0){
		MTTF.approx.CI.M <- data.frame(MTTF.approx.CI.M)
		MTTF.exact.CI.M <- data.frame(MTTF.exact.CI.M)
	}
	if(quanval==1 && cdf_M0_INF>max(q)) tq.CI.M <- data.frame(tq.CI.M)

	row.name <- c("FIM", "OIM(Hessian matrix)", "OIM(score vector)",
			"OIM(robust matrix)")

	name.ind <- c(FIM,OIMA,OIMB,OIMC)
	row.names(eta.CI.M) <- row.name[which(name.ind==1)]
	row.names(seta.CI.M) <- row.name[which(name.ind==1)]
	row.names(sb.CI.M) <- row.name[which(name.ind==1)]
	row.names(se.CI.M) <- row.name[which(name.ind==1)]
	if(MTTFval==1 && Eta.M0>0){
		row.names(MTTF.approx.CI.M) <- row.name[which(name.ind==1)]
		row.names(MTTF.exact.CI.M) <- row.name[which(name.ind==1)]
	}
	if(quanval==1 && cdf_M0_INF>max(q)) row.names(tq.CI.M) <- row.name[which(name.ind==1)]

	names(eta.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	names(seta.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	names(sb.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	names(se.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	if(MTTFval==1 && Eta.M0>0){
		names(MTTF.approx.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
		names(MTTF.exact.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	}
	if(quanval==1 && cdf_M0_INF>max(q)) names(tq.CI.M) <- rep(c("LCI","UCI","LCI.ln","UCI.ln"),length(q))


	cat(paste(rep("#",70),sep="",collapse=""),"\n")
	cat("      Model: Random effect + Brownian motion + Measurement error    ","\n")
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
	cat("eta =",format(Eta.M0,digits=22),"\n")
	if(ind>1) print(eta.CI.M)
	cat("\n")

	cat("sigma_eta =",format(sEta.M0,digits=22),"\n")
	if(ind>1) print(seta.CI.M)
	cat("\n")

	cat("sigma_B =",format(sB.M0,digits=22),"\n")
	if(ind>1) print(sb.CI.M)
	cat("\n")

	cat("sigma_epsilon =",format(sE.M0,digits=22),"\n")
	if(ind>1) print(se.CI.M)
	cat("\n")

	cat("log-likelihood =",format(loglik.M0,digits=22),"\n")
      cat("Optimization Algorithm:",optim.alg,"\n")
	cat("\n")

	if(MTTFval==1){
		if(Eta.M0>0){
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				cat(paste(rep("-",60),sep="",collapse=""),"\n")
				cat(paste("        MTTF and ", 100*(1-alpha),"% confidence interval",sep=""),"\n")
				cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
			}else{
				cat(paste(rep("-",60),sep="",collapse=""),"\n")
				cat("        MTTF","\n")
				cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
			}
			cat("MTTF(exact)=",MTTF.exact.M0,"\n")
			if(ind>1) print(MTTF.exact.CI.M)
			cat("\n")

			cat("MTTF(approx.)=",MTTF.approx.M0,"\n")
			if(ind>1) print(MTTF.approx.CI.M)
			cat("\n")
		}else{
			cat("\n\n\n")
			cat('Warning message: the estimated drift rate is negative, so it is meaningless to calculate the MTTF.','\n\n\n\n')
		}
	}

	if(quanval==1 && cdf_M0_INF>max(q)){
		if(mode(tt)=="numeric"){
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				cat(paste(rep("-",60),sep="",collapse=""),"\n")
				cat(paste("        qth-quantile and ", 100*(1-alpha),"% confidence interval",sep=""),"\n")
				cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
			}else{
				cat(paste(rep("-",60),sep="",collapse=""),"\n")
				cat("        qth-quantile","\n")
				cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
			}
			for(i in 1: length(q)){
				cat("t(",q[i],") =",tt[i],"\n")
				if(ind>1) print(tq.CI.M[,(4*(i-1)+1):(4*i)])
				cat("\n")
			}
		}else{
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				cat(paste(rep("-",60),sep="",collapse=""),"\n")
				cat(paste("        qth-quantile and ", 100*(1-alpha),"% confidence interval",sep=""),"\n")
				cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
			}else{
				cat(paste(rep("-",60),sep="",collapse=""),"\n")
				cat("        qth-quantile","\n")
				cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
			}
			cat("The estimated quantile has failed. Please reset the domain\n to search for the target quantile (or theoretical quantile)\n or check the correctness of the threshold value.","\n\n")
		}
		cat("limit_{ t \ to infty } F_T(t) = ",cdf_M0_INF,"\n\n")
	}
	if(quanval==1 && cdf_M0_INF<max(q)){
		cat("limit_{ t \ to infty } F_T(t) = ",cdf_M0_INF,"\n\n")
	}


	if( (pdfplotval==1 || cdfplotval==1) && pseudoVal=="0" ){
		if(pdfplotval==1){
			x11()
			plot( t2, pdf_M0(t2), type="l", main="PDF estimation",
				 xlab="t", ylab=expression(f[T](t)), lwd=2, col=1 )
		}
		if(cdfplotval==1){
			x11()
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				plot( t2, apply(as.matrix(t2),1,cdf_M0), type="l", ylim=c(0,1), main='CDF and confidence interval with logit transformation',
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}else{
				plot( t2, apply(as.matrix(t2),1,cdf_M0), type="l", ylim=c(0,1), main="CDF estimation",
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}
			ind2 <- 2
			if(FIM==1){
				lines( t2, CDF.FI.M0$lcl.cdf.M0, lty = ind2, col = ind2 )
				lines( t2, CDF.FI.M0$ucl.cdf.M0, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMA==1){
				lines( t2, CDF.OIA.M0$lcl.cdf.M0, lty = ind2, col = ind2 )
				lines( t2, CDF.OIA.M0$ucl.cdf.M0, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMB==1){
				lines( t2, CDF.OIB.M0$lcl.cdf.M0, lty = ind2, col = ind2 )
				lines( t2, CDF.OIB.M0$ucl.cdf.M0, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMC==1){
				lines( t2, CDF.OIC.M0$lcl.cdf.M0, lty = ind2, col = ind2 )
				lines( t2, CDF.OIC.M0$ucl.cdf.M0, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			text1 <- c(
			paste("FIM ", 100*(1-alpha),"% CI",sep=""),paste("OIM (Hessian matrix) ", 100*(1-alpha),"% CI",sep=""),
			paste("OIM (score vector) ", 100*(1-alpha),"% CI",sep=""),
			paste("OIM (robust matrix) ", 100*(1-alpha),"% CI",sep=""))

			ind3 <- which(c(FIM,OIMA,OIMB,OIMC)==1)

			if(length(ind3)!=0){
				legend('bottomright', c("CDF estimation",text1[ind3]),
				col=c(1,1+(1:length(ind3))),lty=c(1,1+(1:length(ind3))),lwd=c(2,rep(1,length(ind3))),box.col="white",bg="white",inset = .01)
			}
		}
	}

	if( (pdfplotval==1 || cdfplotval==1) && pseudoVal=="1" ){
		fit <- function(y,...)	lm(y~0+time,data=data.frame(Data))$coef
		PFT <- W / apply( Data[,-1], 2, fit )
		if(pdfplotval==1){
			x11()
			plot( t2, pdf_M0(t2), type="l", main="PDF estimation",
				 xlab="t", ylab=expression(f[T](t)), lwd=2, col=1 )
			if( all(PFT>0) && Eta.M0>0 ){
				points(sort(PFT), rep(0,length(PFT)), pch=1 )
				legend('topright', c("PDF estimation","Pseudo failure time"),
				col=c(1,1),lty=c(1,-1),lwd=c(2,-1),pch=c(-1,1),merge=TRUE,box.col="white",bg="white",inset = .01)
			}
		}
		if(cdfplotval==1){
			x11()
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				plot( t2, apply(as.matrix(t2),1,cdf_M0), type="l", ylim=c(0,1), main='CDF and confidence interval with logit transformation',
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}else{
				plot( t2, apply(as.matrix(t2),1,cdf_M0), type="l", ylim=c(0,1), main="CDF estimation",
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}
			ind2 <- 2
			if(FIM==1){
				lines( t2, CDF.FI.M0$lcl.cdf.M0, lty = ind2, col = ind2 )
				lines( t2, CDF.FI.M0$ucl.cdf.M0, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMA==1){
				lines( t2, CDF.OIA.M0$lcl.cdf.M0, lty = ind2, col = ind2 )
				lines( t2, CDF.OIA.M0$ucl.cdf.M0, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMB==1){
				lines( t2, CDF.OIB.M0$lcl.cdf.M0, lty = ind2, col = ind2 )
				lines( t2, CDF.OIB.M0$ucl.cdf.M0, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMC==1){
				lines( t2, CDF.OIC.M0$lcl.cdf.M0, lty = ind2, col = ind2 )
				lines( t2, CDF.OIC.M0$ucl.cdf.M0, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			text1 <- c(
			paste("FIM ", 100*(1-alpha),"% CI",sep=""),paste("OIM (Hessian matrix) ", 100*(1-alpha),"% CI",sep=""),
			paste("OIM (score vector) ", 100*(1-alpha),"% CI",sep=""),
			paste("OIM (robust matrix) ", 100*(1-alpha),"% CI",sep=""))

			ind3 <- which(c(FIM,OIMA,OIMB,OIMC)==1)

			if( all(PFT>0) && Eta.M0>0 ){
				y <- (1:length(PFT)-0.5) / length(PFT)
				points(sort(PFT), y, pch=1 )

				if(length(ind3)!=0){
					legend('bottomright', c("CDF estimation",text1[ind3],"Pseudo failure time"),
					col=c(1,1+(1:length(ind3)),1),lty=c(1,1+(1:length(ind3)),-1),lwd=c(2,rep(1,length(ind3)),-1),pch=c(rep(-1,(length(ind3)+1)),1),merge=TRUE,box.col="white",bg="white",inset = .01)
				}else{
					legend('bottomright', c("CDF estimation","Pseudo failure time"),
					col=c(1,1),lty=c(1,-1),lwd=c(2,-1),pch=c(-1,1),merge=TRUE,box.col="white",bg="white",inset = .01)
				}
			}
		}
		if(all(PFT>0) && Eta.M0>0){
			if(ppplotval==1){
				x11()
				plot((1:length(PFT)-0.5)/length(PFT),apply(as.matrix(sort(PFT)),1,cdf_M0),main="Probability-probability plot",xlab="Theoretical cumulative probabilities",ylab="Sample cumulative probabilities",xlim=c(0,1),ylim=c(0,1))
				abline(0,1)
			}
			if(qqplotval==1 && cdf_M0_INF>max(q)){
				qx <- try(apply(as.matrix((1:length(PFT)-0.5)/length(PFT)),1,tq_M0),T)
				qy <- sort(PFT)
				if(mode(qx)=="numeric"){
					x11()
					plot(qx,qy,main="Quantile-quantile plot",xlab="Theoretical quantiles",ylab="Sample quantiles",xlim=c(min(qx,qy),max(qx,qy)),ylim=c(min(qx,qy),max(qx,qy)))
					abline(0,1)
				}else{
					cat(paste(rep("-",60),sep="",collapse=""),"\n")
					cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
					cat("Warning Message: Quantile-Quantile plot is unable to appear\n because the estimated quantile has failed.","\n\n")
				}
			}
			if(kstestval==1 || adtestval==1){
				cat(paste(rep("-",60),sep="",collapse=""),"\n")
				cat("        Goodness-of-fit tests","\n")
				cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
			}
			if(kstestval==1){
				n <- length(PFT)
				value.M0 <- NULL
				for(i in 1:n) value.M0 <- c(value.M0, cdf_M0(sort(PFT)[i]))
				cat("Kolmogorov-Smirnov test :","\n")
				cat("statistic =",(STATISTIC <- max(abs(value.M0 - (1:n - 1)/(n)))),"\n")
				cat("p-value =",(1 - .C("pkolmogorov2x", p = as.double(STATISTIC), as.integer(n), PACKAGE = "stats")$p),"\n\n")
			}
			if(adtestval==1){
				value.M0 <- NULL
				for(i in 1:n) value.M0 <- c(value.M0, cdf_M0(sort(PFT)[i]))
				cat("Anderson-Darling test","\n")
				cat("statistic =",(STATISTIC.M0 <- ad.test(value.M0)$statistic),"\n")
				cat("p-value =",(ad.test(value.M0)$p.value),"\n\n")
			}
		}else{
			cat("Warning message: because of the negative PFT or the negative drift rate estimation,","\n",
			"the functions in the Goodness-of-Fit will not work because it is meaningless to do so.","\n")
		}
	}
}


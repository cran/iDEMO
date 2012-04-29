Model_M5 <-
function(Data, int.eta, int.sb, int.se, W, q, alpha, FIM, OIMA, OIMB, OIMC, t2, pseudoVal, MTTFval, quanval, pdfplotval, cdfplotval, ppplotval, qqplotval, kstestval, adtestval, dataname, optim.alg, env.in){

	Data <- as.matrix(Data)

	time <- Data[,1]

	m <- length(time); n <- length(Data[1,-1])

	Q <- matrix(rep(time, length(time)), ncol = length(time))

	Q[row(Q) > col(Q)] <- sort(Q[row(Q) < col(Q)])

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

			stop('Please try different initial values (and/or the optimization algorithm)\n  to obtain more stable and more accurate estimates of the unknown parameters.')

		}
	}

	if(iter.M5 $ convergence == 0){

		Eta.M5 <- iter.M5 $ par[1]

		sB.M5 <- abs( iter.M5 $ par[2] )

		sE.M5 <- abs( iter.M5 $ par[3] )

	}

	S <- O <- sB.M5^2 * Q + sE.M5^2 * diag(m)

	par.M5 <- c(Eta.M5, sB.M5, sE.M5)

	assign('par.M5', par.M5, env.in)

	loglik.M5 <- - n * m * log( 2 * pi ) / 2 - n * log( det(S) ) / 2 - 

	            sum( diag( as.matrix( t( Data[,-1] - Eta.M5 * time ) ) %*% solve(S) %*% as.matrix( Data[,-1] - Eta.M5 * time ) ) ) / 2

	if(FIM == 1){

		Cov.Matrix.FI.M5 <- solve( FI.M5( Q, O, time ) )

		par.FI.M5 <- par.FICI.M5(alpha, par.M5, Cov.Matrix.FI.M5, n )

		if(MTTFval==1 && Eta.M5>0) MTTF.FI.M5 <- MTTF.FICI.M5( alpha, W, Eta.M5, Cov.Matrix.FI.M5, n )

		if(quanval==1) tq.FI.M5 <- try(tq.FICI.M5.D( q, alpha, W, Eta.M5, sB.M5, Cov.Matrix.FI.M5, n, range(t2) ),T)

		if(cdfplotval==1) CDF.FI.M5 <- cdf.FICI.logit.M5.D( Par.M5 = par.M5, t=t2, Cov.Matrix.FI.M5, alpha=alpha, W=W, n )

	}

	if(OIMA == 1){

		A <- OI.A.M5( Data, Eta.M5, Q, O, time )

		Cov.Matrix.OI.M5.A <- solve(A)

		par.OI.M5.A <- par.OICI.M5( alpha, par.M5, Cov.Matrix.OI.M5.A )

		if(MTTFval==1 && Eta.M5>0) MTTF.OIA.M5 <- MTTF.OICI.M5( alpha, W, Eta.M5, Cov.Matrix.OI.M5.A )

		if(quanval==1) tq.OIA.M5 <- try(tq.OICI.M5.D( q, alpha, W, Eta.M5, sB.M5, Cov.Matrix.OI.M5.A, range(t2) ),T)

		if(cdfplotval==1) CDF.OIA.M5 <- cdf.OICI.logit.M5.D( Par.M5 = par.M5, t=t2, Cov.Matrix.OI.M5.A, alpha=alpha, W=W )

	}

	if(OIMB == 1){

		B <- OI.B.M5( Data, Eta.M5, Q, O, time )

		Cov.Matrix.OI.M5.B <- solve(B)

		par.OI.M5.B <- par.OICI.M5( alpha, par.M5, Cov.Matrix.OI.M5.B )

		if(MTTFval==1 && Eta.M5>0) MTTF.OIB.M5 <- MTTF.OICI.M5( alpha, W, Eta.M5, Cov.Matrix.OI.M5.B )

		if(quanval==1) tq.OIB.M5 <- try(tq.OICI.M5.D( q, alpha, W, Eta.M5, sB.M5, Cov.Matrix.OI.M5.B, range(t2) ),T)

		if(cdfplotval==1) CDF.OIB.M5 <- cdf.OICI.logit.M5.D( Par.M5 = par.M5, t=t2, Cov.Matrix.OI.M5.B, alpha=alpha, W=W )

	}

	if(OIMC == 1){

		A <- OI.A.M5( Data, Eta.M5, Q, O, time )

		B <- OI.B.M5( Data, Eta.M5, Q, O, time )

		Cov.Matrix.OI.M5.C <- solve(A) %*% B %*% solve(A)

		par.OI.M5.C <- par.OICI.M5( alpha, par.M5, Cov.Matrix.OI.M5.C )

		if(MTTFval==1 && Eta.M5>0) MTTF.OIC.M5 <- MTTF.OICI.M5( alpha, W, Eta.M5, Cov.Matrix.OI.M5.C )

		if(quanval==1) tq.OIC.M5 <- try(tq.OICI.M5.D( q, alpha, W, Eta.M5, sB.M5, Cov.Matrix.OI.M5.C, range(t2) ),T)

		if(cdfplotval==1) CDF.OIC.M5 <- cdf.OICI.logit.M5.D( Par.M5 = par.M5, t=t2, Cov.Matrix.OI.M5.C, alpha=alpha, W=W )

	}

	if(MTTFval==1 && Eta.M5>0){

		MTTF.approx.M5 <- MTTF.M5( par.M5[1], par.M5[2], W )[1]

		MTTF.exact.M5 <- MTTF.M5( par.M5[1], par.M5[2], W )[2]

	}

	if(quanval==1){

		fun_tp <- function(p5){

			t.percM5 <- function(t){

				x0 <- 2 * Eta.M5 * W / sB.M5^2

				x2 <- - ( Eta.M5 * t + W ) / sqrt( sB.M5^2 * t )

				if( x2< -38){
		      	      pnorm( ( Eta.M5 * t - W ) / sqrt( sB.M5^2 * t ) ) - p5
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
		      	      pnorm( ( Eta.M5 * t - W ) / sqrt( sB.M5^2 * t ) ) + V - p5
				}

			}
  
			tp <- uniroot( t.percM5, range(t2) )$root

		}

		tt <- try(apply( t(q), 2, fun_tp),T)

	}

	pdf_M5 <- function(t){

			sqrt( W^2 / ( 2 * pi * t^3 * sB.M5^2 ) ) *

			exp( - ( W - Eta.M5 * t )^2 / ( 2 * t * sB.M5^2 ) )

			}

	cdf_M5 <- function(t){
		x0 <- 2 * Eta.M5 * W / sB.M5^2

		x2 <- - ( Eta.M5 * t + W ) / sqrt( sB.M5^2 * t )

		if( x2< -38){
      	      pnorm( ( Eta.M5 * t - W ) / sqrt( sB.M5^2 * t ) )
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
      	      pnorm( ( Eta.M5 * t - W ) / sqrt( sB.M5^2 * t ) ) + V
		}
	}

	tq_M5 <- function(p5){

		t.percM5 <- function(t){

			x0 <- 2 * Eta.M5 * W / sB.M5^2

			x2 <- - ( Eta.M5 * t + W ) / sqrt( sB.M5^2 * t )

			if( x2< -38){
	      	      pnorm( ( Eta.M5 * t - W ) / sqrt( sB.M5^2 * t ) ) - p5
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
	      	      pnorm( ( Eta.M5 * t - W ) / sqrt( sB.M5^2 * t ) ) + V - p5
			}
		}
  
		uniroot( t.percM5, range(t2) )$root
	}

	eta.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	sb.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	se.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	MTTF.approx.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	MTTF.exact.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	tq.CI.M <- array(0,c(sum(FIM,OIMA,OIMB,OIMC),4,length(q)))

	ind <- 1

	if(FIM==1){
		eta.CI.M[ind,] <- t(c(par.FI.M5$CI.Eta.M5,par.FI.M5$CI.ln.Eta.M5))
		sb.CI.M[ind,] <- t(c(par.FI.M5$CI.sB.M5,par.FI.M5$CI.ln.sB.M5))
		se.CI.M[ind,] <- t(c(par.FI.M5$CI.sE.M5,par.FI.M5$CI.ln.sE.M5))
		if(MTTFval==1 && Eta.M5>0) MTTF.approx.CI.M[ind,] <- t(c(MTTF.FI.M5$CI.MTTF.approx.M5, MTTF.FI.M5$CI.ln.MTTF.approx.M5))
		if(MTTFval==1 && Eta.M5>0) MTTF.exact.CI.M[ind,] <- t(c(MTTF.FI.M5$CI.MTTF.exact.M5, MTTF.FI.M5$CI.ln.MTTF.exact.M5))
		if(quanval==1 && length(tq.FI.M5)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.FI.M5[i,2:5])
			}
		}
		ind <- ind + 1
	}

	if(OIMA==1){
		eta.CI.M[ind,] <- t(c(par.OI.M5.A$CI.Eta.M5,par.OI.M5.A$CI.ln.Eta.M5))
		sb.CI.M[ind,] <- t(c(par.OI.M5.A$CI.sB.M5,par.OI.M5.A$CI.ln.sB.M5))
		se.CI.M[ind,] <- t(c(par.OI.M5.A$CI.sE.M5,par.OI.M5.A$CI.ln.sE.M5))
		if(MTTFval==1 && Eta.M5>0) MTTF.approx.CI.M[ind,] <- t(c(MTTF.OIA.M5$CI.MTTF.approx.M5, MTTF.OIA.M5$CI.ln.MTTF.approx.M5))
		if(MTTFval==1 && Eta.M5>0) MTTF.exact.CI.M[ind,] <- t(c(MTTF.OIA.M5$CI.MTTF.exact.M5, MTTF.OIA.M5$CI.ln.MTTF.exact.M5))
		if(quanval==1 && length(tq.OIA.M5)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.OIA.M5[i,2:5])
			}
		}
		ind <- ind + 1
	}

	if(OIMB==1){
		eta.CI.M[ind,] <- t(c(par.OI.M5.B$CI.Eta.M5,par.OI.M5.B$CI.ln.Eta.M5))
		sb.CI.M[ind,] <- t(c(par.OI.M5.B$CI.sB.M5,par.OI.M5.B$CI.ln.sB.M5))
		se.CI.M[ind,] <- t(c(par.OI.M5.B$CI.sE.M5,par.OI.M5.B$CI.ln.sE.M5))
		if(MTTFval==1 && Eta.M5>0) MTTF.approx.CI.M[ind,] <- t(c(MTTF.OIB.M5$CI.MTTF.approx.M5, MTTF.OIB.M5$CI.ln.MTTF.approx.M5))
		if(MTTFval==1 && Eta.M5>0) MTTF.exact.CI.M[ind,] <- t(c(MTTF.OIB.M5$CI.MTTF.exact.M5, MTTF.OIB.M5$CI.ln.MTTF.exact.M5))
		if(quanval==1 && length(tq.OIB.M5)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.OIB.M5[i,2:5])
			}
		}
		ind <- ind + 1
	}

	if(OIMC==1){
		eta.CI.M[ind,] <- t(c(par.OI.M5.C$CI.Eta.M5,par.OI.M5.C$CI.ln.Eta.M5))
		sb.CI.M[ind,] <- t(c(par.OI.M5.C$CI.sB.M5,par.OI.M5.C$CI.ln.sB.M5))
		se.CI.M[ind,] <- t(c(par.OI.M5.C$CI.sE.M5,par.OI.M5.C$CI.ln.sE.M5))
		if(MTTFval==1 && Eta.M5>0) MTTF.approx.CI.M[ind,] <- t(c(MTTF.OIC.M5$CI.MTTF.approx.M5, MTTF.OIC.M5$CI.ln.MTTF.approx.M5))
		if(MTTFval==1 && Eta.M5>0) MTTF.exact.CI.M[ind,] <- t(c(MTTF.OIC.M5$CI.MTTF.exact.M5, MTTF.OIC.M5$CI.ln.MTTF.exact.M5))
		if(quanval==1 && length(tq.OIC.M5)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.OIC.M5[i,2:5])
			}
		}
		ind <- ind + 1
	}


	eta.CI.M <- data.frame(eta.CI.M)
	sb.CI.M <- data.frame(sb.CI.M)
	se.CI.M <- data.frame(se.CI.M)
	if(MTTFval==1 && Eta.M5>0) MTTF.approx.CI.M <- data.frame(MTTF.approx.CI.M)
	if(MTTFval==1 && Eta.M5>0) MTTF.exact.CI.M <- data.frame(MTTF.exact.CI.M)
	if(quanval==1) tq.CI.M <- data.frame(tq.CI.M)

	row.name <- c("FIM", "OIM(Hessian matrix)", "OIM(score vector)",
			"OIM(robust matrix)")

	name.ind <- c(FIM,OIMA,OIMB,OIMC)
	row.names(eta.CI.M) <- row.name[which(name.ind==1)]
	row.names(sb.CI.M) <- row.name[which(name.ind==1)]
	row.names(se.CI.M) <- row.name[which(name.ind==1)]
	if(MTTFval==1 && Eta.M5>0) row.names(MTTF.approx.CI.M) <- row.name[which(name.ind==1)]
	if(MTTFval==1 && Eta.M5>0) row.names(MTTF.exact.CI.M) <- row.name[which(name.ind==1)]
	if(quanval==1) row.names(tq.CI.M) <- row.name[which(name.ind==1)]

	names(eta.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	names(sb.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	names(se.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	if(MTTFval==1 && Eta.M5>0) names(MTTF.approx.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	if(MTTFval==1 && Eta.M5>0) names(MTTF.exact.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	if(quanval==1) names(tq.CI.M) <- rep(c("LCI","UCI","LCI.ln","UCI.ln"),length(q))

	cat(paste(rep("#",70),sep="",collapse=""),"\n")
	cat("      Model: Brownian motion + Measurement error    ","\n")
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

	cat("eta =",format(Eta.M5,digits=22),"\n")
	if(ind>1) print(eta.CI.M)
	cat("\n")

	cat("sigma_B =",format(sB.M5,digits=22),"\n")
	if(ind>1) print(sb.CI.M)
	cat("\n")

	cat("sigma_epsilon =",format(sE.M5,digits=22),"\n")
	if(ind>1) print(se.CI.M)
	cat("\n")

	cat("log-likelihood =",format(loglik.M5,digits=22),"\n")
      cat("Optimization Algorithm:",optim.alg,"\n")
	cat("\n")

	if(MTTFval==1){
		if(Eta.M5>0){
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				cat(paste(rep("-",60),sep="",collapse=""),"\n")
				cat(paste("        MTTF and ", 100*(1-alpha),"% confidence interval",sep=""),"\n")
				cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
			}else{
				cat(paste(rep("-",60),sep="",collapse=""),"\n")
				cat("        MTTF","\n")
				cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
			}
			cat("MTTF(exact)=",MTTF.exact.M5,"\n")
			if(ind>1) print(MTTF.exact.CI.M)
			cat("\n")

			cat("MTTF(approx.)=",MTTF.approx.M5,"\n")
			if(ind>1) print(MTTF.approx.CI.M)
			cat("\n")
		}else{
			cat("\n\n\n")
			cat('Warning message: the estimated drift rate is negative, so it is meaningless to calculate the MTTF.','\n\n\n\n')
		}
	}

	if(quanval==1){
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
		cat("limit_{ t \ to infty } F_T(t) = ",1,"\n\n")
	}

	if( (pdfplotval==1 || cdfplotval==1) && pseudoVal=="0" ){
		if(pdfplotval==1){
			x11()
			plot( t2, pdf_M5(t2), type="l", main="PDF estimation",
				 xlab="t", ylab=expression(f[T](t)), lwd=2, col=1 )
		}
		if(cdfplotval==1){
			x11()
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				plot( t2, apply(as.matrix(t2),1,cdf_M5), type="l", ylim=c(0,1), main='CDF and confidence interval with logit transformation',
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}else{
				plot( t2, apply(as.matrix(t2),1,cdf_M5), type="l", ylim=c(0,1), main="CDF estimation",
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}
			ind2 <- 2
			if(FIM==1){
				lines( t2, CDF.FI.M5$lcl.cdf.M5, lty = ind2, col = ind2 )
				lines( t2, CDF.FI.M5$ucl.cdf.M5, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMA==1){
				lines( t2, CDF.OIA.M5$lcl.cdf.M5, lty = ind2, col = ind2 )
				lines( t2, CDF.OIA.M5$ucl.cdf.M5, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMB==1){
				lines( t2, CDF.OIB.M5$lcl.cdf.M5, lty = ind2, col = ind2 )
				lines( t2, CDF.OIB.M5$ucl.cdf.M5, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMC==1){
				lines( t2, CDF.OIC.M5$lcl.cdf.M5, lty = ind2, col = ind2 )
				lines( t2, CDF.OIC.M5$ucl.cdf.M5, lty = ind2, col = ind2 )
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
			plot( t2, pdf_M5(t2), type="l", main="PDF estimation",
				 xlab="t", ylab=expression(f[T](t)), lwd=2, col=1 )
			if( all(PFT>0) && Eta.M5>0 ){
				points(sort(PFT), rep(0,length(PFT)), pch=1 )
				legend('topright', c("PDF estimation","Pseudo failure time"),
				col=c(1,1),lty=c(1,-1),lwd=c(2,-1),pch=c(-1,1),merge=TRUE,box.col="white",bg="white",inset = .01)
			}
		}
		if(cdfplotval==1){
			x11()
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				plot( t2, apply(as.matrix(t2),1,cdf_M5), type="l", ylim=c(0,1), main='CDF and confidence interval with logit transformation',
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}else{
				plot( t2, apply(as.matrix(t2),1,cdf_M5), type="l", ylim=c(0,1), main="CDF estimation",
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}
			ind2 <- 2
			if(FIM==1){
				lines( t2, CDF.FI.M5$lcl.cdf.M5, lty = ind2, col = ind2 )
				lines( t2, CDF.FI.M5$ucl.cdf.M5, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMA==1){
				lines( t2, CDF.OIA.M5$lcl.cdf.M5, lty = ind2, col = ind2 )
				lines( t2, CDF.OIA.M5$ucl.cdf.M5, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMB==1){
				lines( t2, CDF.OIB.M5$lcl.cdf.M5, lty = ind2, col = ind2 )
				lines( t2, CDF.OIB.M5$ucl.cdf.M5, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMC==1){
				lines( t2, CDF.OIC.M5$lcl.cdf.M5, lty = ind2, col = ind2 )
				lines( t2, CDF.OIC.M5$ucl.cdf.M5, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			text1 <- c(
			paste("FIM ", 100*(1-alpha),"% CI",sep=""),paste("OIM (Hessian matrix) ", 100*(1-alpha),"% CI",sep=""),
			paste("OIM (score vector) ", 100*(1-alpha),"% CI",sep=""),
			paste("OIM (robust matrix) ", 100*(1-alpha),"% CI",sep=""))

			ind3 <- which(c(FIM,OIMA,OIMB,OIMC)==1)

			if( all(PFT>0) && Eta.M5>0 ){
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
		if(all(PFT>0) && Eta.M5>0){
			if(ppplotval==1){
				x11()
				plot((1:length(PFT)-0.5)/length(PFT),apply(as.matrix(sort(PFT)),1,cdf_M5),main="Probability-probability plot",xlab="Theoretical cumulative probabilities",ylab="Sample cumulative probabilities",xlim=c(0,1),ylim=c(0,1))
				abline(0,1)
			}
			if(qqplotval==1){
				qx <- try(apply(as.matrix((1:length(PFT)-0.5)/length(PFT)),1,tq_M5),T)
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
				value.M5 <- NULL
				for(i in 1:n) value.M5 <- c(value.M5, cdf_M5(sort(PFT)[i]))
				cat("Kolmogorov-Smirnov test :","\n")
				cat("statistic =",(STATISTIC <- max(abs(value.M5 - (1:n - 1)/(n)))),"\n")
				cat("p-value =",(1 - .C("pkolmogorov2x", p = as.double(STATISTIC), as.integer(n), PACKAGE = "stats")$p),"\n\n")
			}
			if(adtestval==1){
				value.M5 <- NULL
				for(i in 1:n) value.M5 <- c(value.M5, cdf_M5(sort(PFT)[i]))
				cat("Anderson-Darling test","\n")
				cat("statistic =",(STATISTIC.M5 <- ad.test(value.M5)$statistic),"\n")
				cat("p-value =",(ad.test(value.M5)$p.value),"\n\n")
			}
		}else{
			cat("Warning message: because of the negative PFT or the negative drift rate estimation,","\n",
			"the functions in the Goodness-of-Fit will not work because it is meaningless to do so.","\n")
		}
	}
}


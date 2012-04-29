Model_M2 <-
function(Data, int.eta, int.sb, W, q, alpha, FIM, OIMA, OIMB, OIMC, t2, pseudoVal, MTTFval, quanval, pdfplotval, cdfplotval, ppplotval, qqplotval, kstestval, adtestval, dataname, env.in){

	Data <- as.matrix(Data)

	time <- Data[,1]

	m <- length(time); n <- length(Data[1,-1])

	Q <- matrix(rep(time, length(time)), ncol = length(time))

	Q[row(Q) > col(Q)] <- sort(Q[row(Q) < col(Q)])

  	Q.inv <- matrix( 0, m, m )

	a <- numeric( length = m )

	a[1] <- 1 / time[1]

	for(i in 2:m) a[i] <- 1 / ( time[i] - time[i-1] )

  	for(i in 2:(m-1)) Q.inv[ i, c( i-1, i, i+1 ) ] <- c( -a[i], a[i]+a[i+1], -a[i+1] )

	Q.inv[ 1, c( 1, 2 ) ] <- c( a[1]+a[2], -a[2] )

 	Q.inv[ m, c( m-1, m ) ] <- c( -a[m], a[m] )


	Eta.M2 <- sum( t(time) %*% Q.inv %*% Data[,-1] ) / ( n * t(time) %*% Q.inv %*% time )

	sB.M2 <- sqrt( sum( diag( t( Data[,-1] - Eta.M2 * time ) %*% Q.inv %*% ( Data[,-1] - Eta.M2 * time ) ) ) / ( n * m ) )

	S <- O <- sB.M2^2 *Q

	par.M2 <- c(Eta.M2, sB.M2)

	assign('par.M2', par.M2, env.in)

	loglik.M2 <- - n * m * log( 2 * pi ) / 2 - n * log( det( sB.M2^2 * Q ) ) / 2 - 

            sum( diag( as.matrix( t( Data[,-1] - Eta.M2 * time ) ) %*% ( Q.inv / sB.M2^2 ) %*% as.matrix( Data[,-1] - Eta.M2 * time ) ) ) / 2

	if(FIM == 1){

		Cov.Matrix.FI.M2 <- solve( FI.M2( sB.M2, Q.inv, time, m ) )

		par.FI.M2 <- par.FICI.M2(alpha, par.M2, Cov.Matrix.FI.M2, n )

		if(MTTFval==1 && Eta.M2>0) MTTF.FI.M2 <- MTTF.FICI.M2( alpha, W, Eta.M2, Cov.Matrix.FI.M2, n )

		if(quanval==1) tq.FI.M2 <- try(tq.FICI.M2.D( q, alpha, W, Eta.M2, sB.M2, Cov.Matrix.FI.M2, n, range(t2) ),T)

		if(cdfplotval==1) CDF.FI.M2 <- cdf.FICI.logit.M2.D( Par.M2 = par.M2, t=t2, Cov.Matrix.FI.M2, alpha=alpha, W=W, n )

	}

	if(OIMA == 1){

		A <- OI.A.M2( Data, Eta.M2, sB.M2, S, Q.inv, time, m )

		Cov.Matrix.OI.M2.A <- solve(A)

		par.OI.M2.A <- par.OICI.M2( alpha, par.M2, Cov.Matrix.OI.M2.A )

		if(MTTFval==1 && Eta.M2>0) MTTF.OIA.M2 <- MTTF.OICI.M2( alpha, W, Eta.M2, Cov.Matrix.OI.M2.A )

		if(quanval==1) tq.OIA.M2 <- try(tq.OICI.M2.D( q, alpha, W, Eta.M2, sB.M2, Cov.Matrix.OI.M2.A, range(t2) ),T)

		if(cdfplotval==1) CDF.OIA.M2 <- cdf.OICI.logit.M2.D( Par.M2 = par.M2, t=t2, Cov.Matrix.OI.M2.A, alpha=alpha, W=W )

	}

	if(OIMB == 1){

		B <- OI.B.M2( Data, Eta.M2, sB.M2, O, S, Q.inv, time, m )

		Cov.Matrix.OI.M2.B <- solve(B)

		par.OI.M2.B <- par.OICI.M2( alpha, par.M2, Cov.Matrix.OI.M2.B )

		if(MTTFval==1 && Eta.M2>0) MTTF.OIB.M2 <- MTTF.OICI.M2( alpha, W, Eta.M2, Cov.Matrix.OI.M2.B )

		if(quanval==1) tq.OIB.M2 <- try(tq.OICI.M2.D( q, alpha, W, Eta.M2, sB.M2, Cov.Matrix.OI.M2.B, range(t2) ),T)

		if(cdfplotval==1) CDF.OIB.M2 <- cdf.OICI.logit.M2.D( Par.M2 = par.M2, t=t2, Cov.Matrix.OI.M2.B, alpha=alpha, W=W )

	}

	if(OIMC == 1){

		A <- OI.A.M2( Data, Eta.M2, sB.M2, S, Q.inv, time, m )

		B <- OI.B.M2( Data, Eta.M2, sB.M2, O, S, Q.inv, time, m )

		Cov.Matrix.OI.M2.C <- solve(A) %*% B %*% solve(A)

		par.OI.M2.C <- par.OICI.M2( alpha, par.M2, Cov.Matrix.OI.M2.C )

		if(MTTFval==1 && Eta.M2>0) MTTF.OIC.M2 <- MTTF.OICI.M2( alpha, W, Eta.M2, Cov.Matrix.OI.M2.C )

		if(quanval==1) tq.OIC.M2 <- try(tq.OICI.M2.D( q, alpha, W, Eta.M2, sB.M2, Cov.Matrix.OI.M2.C, range(t2) ),T)

		if(cdfplotval==1) CDF.OIC.M2 <- cdf.OICI.logit.M2.D( Par.M2 = par.M2, t=t2, Cov.Matrix.OI.M2.C, alpha=alpha, W=W )

	}

	if(MTTFval==1 && Eta.M2>0){

		MTTF.approx.M2 <- MTTF.M2( par.M2[1], par.M2[2], W )[1]

		MTTF.exact.M2 <- MTTF.M2( par.M2[1], par.M2[2], W )[2]

	}

	if(quanval==1){

		fun_tp <- function(p2){

			t.percM2 <- function(t){

				x0 <- 2 * Eta.M2 * W / sB.M2^2

				x2 <- - ( Eta.M2 * t + W ) / sqrt( sB.M2^2 * t )

				if( x2< -38){
      			      pnorm( ( Eta.M2 * t - W ) / sqrt( sB.M2^2 * t ) ) - p2
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
      		      	pnorm( ( Eta.M2 * t - W ) / sqrt( sB.M2^2 * t ) ) + V - p2
				}
			}
  
			tp <- uniroot( t.percM2, range(t2) )$root

		}

		tt <- try(apply( t(q), 2, fun_tp),T)

	}

	pdf_M2 <- function(t){

			sqrt( W^2 / ( 2 * pi * t^3 * sB.M2^2 ) ) *

			exp( - ( W - Eta.M2 * t )^2 / ( 2 * t * sB.M2^2 ) )

			}

	cdf_M2 <- function(t){

		x0 <- 2 * Eta.M2 * W / sB.M2^2

		x2 <- - ( Eta.M2 * t + W ) / sqrt( sB.M2^2 * t )

		if( x2< -38){
      	      pnorm( ( Eta.M2 * t - W ) / sqrt( sB.M2^2 * t ) )
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
      	      pnorm( ( Eta.M2 * t - W ) / sqrt( sB.M2^2 * t ) ) + V
		}

	}

	tq_M2 <- function(p2){

		t.percM2 <- function(t){

			x0 <- 2 * Eta.M2 * W / sB.M2^2

			x2 <- - ( Eta.M2 * t + W ) / sqrt( sB.M2^2 * t )

			if( x2< -38){
	      	      pnorm( ( Eta.M2 * t - W ) / sqrt( sB.M2^2 * t ) ) - p2
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
	      	      pnorm( ( Eta.M2 * t - W ) / sqrt( sB.M2^2 * t ) ) + V - p2
			}
		}
  
		uniroot( t.percM2, range(t2) )$root

	}

	eta.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	sb.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	MTTF.approx.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	MTTF.exact.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	tq.CI.M <- array(0,c(sum(FIM,OIMA,OIMB,OIMC),4,length(q)))

	ind <- 1

	if(FIM==1){
		eta.CI.M[ind,] <- t(c(par.FI.M2$CI.Eta.M2,par.FI.M2$CI.ln.Eta.M2))
		sb.CI.M[ind,] <- t(c(par.FI.M2$CI.sB.M2,par.FI.M2$CI.ln.sB.M2))
		if(MTTFval==1 && Eta.M2>0) MTTF.approx.CI.M[ind,] <- t(c(MTTF.FI.M2$CI.MTTF.approx.M2, MTTF.FI.M2$CI.ln.MTTF.approx.M2))
		if(MTTFval==1 && Eta.M2>0) MTTF.exact.CI.M[ind,] <- t(c(MTTF.FI.M2$CI.MTTF.exact.M2, MTTF.FI.M2$CI.ln.MTTF.exact.M2))
		if(quanval==1 && length(tq.FI.M2)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.FI.M2[i,2:5])
			}
		}
		ind <- ind + 1
	}

	if(OIMA==1){
		eta.CI.M[ind,] <- t(c(par.OI.M2.A$CI.Eta.M2,par.OI.M2.A$CI.ln.Eta.M2))
		sb.CI.M[ind,] <- t(c(par.OI.M2.A$CI.sB.M2,par.OI.M2.A$CI.ln.sB.M2))
		if(MTTFval==1 && Eta.M2>0) MTTF.approx.CI.M[ind,] <- t(c(MTTF.OIA.M2$CI.MTTF.approx.M2, MTTF.OIA.M2$CI.ln.MTTF.approx.M2))
		if(MTTFval==1 && Eta.M2>0) MTTF.exact.CI.M[ind,] <- t(c(MTTF.OIA.M2$CI.MTTF.exact.M2, MTTF.OIA.M2$CI.ln.MTTF.exact.M2))
		if(quanval==1 && length(tq.OIA.M2)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.OIA.M2[i,2:5])
			}
		}
		ind <- ind + 1
	}

	if(OIMB==1){
		eta.CI.M[ind,] <- t(c(par.OI.M2.B$CI.Eta.M2,par.OI.M2.B$CI.ln.Eta.M2))
		sb.CI.M[ind,] <- t(c(par.OI.M2.B$CI.sB.M2,par.OI.M2.B$CI.ln.sB.M2))
		if(MTTFval==1 && Eta.M2>0) MTTF.approx.CI.M[ind,] <- t(c(MTTF.OIB.M2$CI.MTTF.approx.M2, MTTF.OIB.M2$CI.ln.MTTF.approx.M2))
		if(MTTFval==1 && Eta.M2>0) MTTF.exact.CI.M[ind,] <- t(c(MTTF.OIB.M2$CI.MTTF.exact.M2, MTTF.OIB.M2$CI.ln.MTTF.exact.M2))
		if(quanval==1 && length(tq.OIB.M2)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.OIB.M2[i,2:5])
			}
		}
		ind <- ind + 1
	}

	if(OIMC==1){
		eta.CI.M[ind,] <- t(c(par.OI.M2.C$CI.Eta.M2,par.OI.M2.C$CI.ln.Eta.M2))
		sb.CI.M[ind,] <- t(c(par.OI.M2.C$CI.sB.M2,par.OI.M2.C$CI.ln.sB.M2))
		if(MTTFval==1 && Eta.M2>0) MTTF.approx.CI.M[ind,] <- t(c(MTTF.OIC.M2$CI.MTTF.approx.M2, MTTF.OIC.M2$CI.ln.MTTF.approx.M2))
		if(MTTFval==1 && Eta.M2>0) MTTF.exact.CI.M[ind,] <- t(c(MTTF.OIC.M2$CI.MTTF.exact.M2, MTTF.OIC.M2$CI.ln.MTTF.exact.M2))
		if(quanval==1 && length(tq.OIC.M2)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.OIC.M2[i,2:5])
			}
		}
		ind <- ind + 1
	}

	eta.CI.M <- data.frame(eta.CI.M)
	sb.CI.M <- data.frame(sb.CI.M)
	if(MTTFval==1 && Eta.M2>0) MTTF.approx.CI.M <- data.frame(MTTF.approx.CI.M)
	if(MTTFval==1 && Eta.M2>0) MTTF.exact.CI.M <- data.frame(MTTF.exact.CI.M)
	if(quanval==1) tq.CI.M <- data.frame(tq.CI.M)

	row.name <- c("FIM", "OIM(Hessian matrix)", "OIM(score vector)",
			"OIM(robust matrix)")

	name.ind <- c(FIM,OIMA,OIMB,OIMC)
	row.names(eta.CI.M) <- row.name[which(name.ind==1)]
	row.names(sb.CI.M) <- row.name[which(name.ind==1)]
	if(MTTFval==1 && Eta.M2>0) row.names(MTTF.approx.CI.M) <- row.name[which(name.ind==1)]
	if(MTTFval==1 && Eta.M2>0) row.names(MTTF.exact.CI.M) <- row.name[which(name.ind==1)]
	if(quanval==1) row.names(tq.CI.M) <- row.name[which(name.ind==1)]

	names(eta.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	names(sb.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	if(MTTFval==1 && Eta.M2>0) names(MTTF.approx.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	if(MTTFval==1 && Eta.M2>0) names(MTTF.exact.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	if(quanval==1) names(tq.CI.M) <- rep(c("LCI","UCI","LCI.ln","UCI.ln"),length(q))


	cat(paste(rep("#",70),sep="",collapse=""),"\n")
	cat("      Model: Brownian motion    ","\n")
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

	cat("eta =",format(Eta.M2,digits=22),"\n")
	if(ind>1) print(eta.CI.M)
	cat("\n")

	cat("sigma_B =",format(sB.M2,digits=22),"\n")
	if(ind>1) print(sb.CI.M)
	cat("\n")

	cat("log-likelihood =",format(loglik.M2,digits=22),"\n")
	cat("\n")

	if(MTTFval==1){
		if(Eta.M2>0){
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				cat(paste(rep("-",60),sep="",collapse=""),"\n")
				cat(paste("        MTTF and ", 100*(1-alpha),"% confidence interval",sep=""),"\n")
				cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
			}else{
				cat(paste(rep("-",60),sep="",collapse=""),"\n")
				cat("        MTTF","\n")
				cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
			}
			cat("MTTF(exact)=",MTTF.exact.M2,"\n")
			if(ind>1) print(MTTF.exact.CI.M)
			cat("\n")

			cat("MTTF(approx.)=",MTTF.approx.M2,"\n")
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
			plot( t2, pdf_M2(t2), type="l", main="PDF estimation",
				 xlab="t", ylab=expression(f[T](t)), lwd=2, col=1 )
		}
		if(cdfplotval==1){
			x11()
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				plot( t2, apply(as.matrix(t2),1,cdf_M2), type="l", ylim=c(0,1), main='CDF and confidence interval with logit transformation',
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}else{
				plot( t2, apply(as.matrix(t2),1,cdf_M2), type="l", ylim=c(0,1), main="CDF estimation",
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}
			ind2 <- 2
			if(FIM==1){
				lines( t2, CDF.FI.M2$lcl.cdf.M2, lty = ind2, col = ind2 )
				lines( t2, CDF.FI.M2$ucl.cdf.M2, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMA==1){
				lines( t2, CDF.OIA.M2$lcl.cdf.M2, lty = ind2, col = ind2 )
				lines( t2, CDF.OIA.M2$ucl.cdf.M2, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMB==1){
				lines( t2, CDF.OIB.M2$lcl.cdf.M2, lty = ind2, col = ind2 )
				lines( t2, CDF.OIB.M2$ucl.cdf.M2, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMC==1){
				lines( t2, CDF.OIC.M2$lcl.cdf.M2, lty = ind2, col = ind2 )
				lines( t2, CDF.OIC.M2$ucl.cdf.M2, lty = ind2, col = ind2 )
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
			plot( t2, pdf_M2(t2), type="l", main="PDF estimation",
				 xlab="t", ylab=expression(f[T](t)), lwd=2, col=1 )
			if( all(PFT>0) && Eta.M2>0 ){
				points(sort(PFT), rep(0,length(PFT)), pch=1 )
				legend('topright', c("PDF estimation","Pseudo failure time"),
				col=c(1,1),lty=c(1,-1),lwd=c(2,-1),pch=c(-1,1),merge=TRUE,box.col="white",bg="white",inset = .01)
			}
		}
		if(cdfplotval==1){
			x11()
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				plot( t2, apply(as.matrix(t2),1,cdf_M2), type="l", ylim=c(0,1), main='CDF and confidence interval with logit transformation',
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}else{
				plot( t2, apply(as.matrix(t2),1,cdf_M2), type="l", ylim=c(0,1), main="CDF estimation",
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}
			ind2 <- 2
			if(FIM==1){
				lines( t2, CDF.FI.M2$lcl.cdf.M2, lty = ind2, col = ind2 )
				lines( t2, CDF.FI.M2$ucl.cdf.M2, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMA==1){
				lines( t2, CDF.OIA.M2$lcl.cdf.M2, lty = ind2, col = ind2 )
				lines( t2, CDF.OIA.M2$ucl.cdf.M2, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMB==1){
				lines( t2, CDF.OIB.M2$lcl.cdf.M2, lty = ind2, col = ind2 )
				lines( t2, CDF.OIB.M2$ucl.cdf.M2, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMC==1){
				lines( t2, CDF.OIC.M2$lcl.cdf.M2, lty = ind2, col = ind2 )
				lines( t2, CDF.OIC.M2$ucl.cdf.M2, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			text1 <- c(
			paste("FIM ", 100*(1-alpha),"% CI",sep=""),paste("OIM (Hessian matrix) ", 100*(1-alpha),"% CI",sep=""),
			paste("OIM (score vector) ", 100*(1-alpha),"% CI",sep=""),
			paste("OIM (robust matrix) ", 100*(1-alpha),"% CI",sep=""))

			ind3 <- which(c(FIM,OIMA,OIMB,OIMC)==1)

			if( all(PFT>0) && Eta.M2>0 ){
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
		if(all(PFT>0) && Eta.M2>0){
			if(ppplotval==1){
				x11()
				plot((1:length(PFT)-0.5)/length(PFT),apply(as.matrix(sort(PFT)),1,cdf_M2),main="Probability-probability plot",xlab="Theoretical cumulative probabilities",ylab="Sample cumulative probabilities",xlim=c(0,1),ylim=c(0,1))
				abline(0,1)
			}
			if(qqplotval==1){
				qx <- try(apply(as.matrix((1:length(PFT)-0.5)/length(PFT)),1,tq_M2),T)
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
				value.M2 <- NULL
				for(i in 1:n) value.M2 <- c(value.M2, cdf_M2(sort(PFT)[i]))
				cat("Kolmogorov-Smirnov test :","\n")
				cat("statistic =",(STATISTIC <- max(abs(value.M2 - (1:n - 1)/(n)))),"\n")
				cat("p-value =",(1 - .C("pkolmogorov2x", p = as.double(STATISTIC), as.integer(n), PACKAGE = "stats")$p),"\n\n")
			}
			if(adtestval==1){
				value.M2 <- NULL
				for(i in 1:n) value.M2 <- c(value.M2, cdf_M2(sort(PFT)[i]))
				cat("Anderson-Darling test","\n")
				cat("statistic =",(STATISTIC.M2 <- ad.test(value.M2)$statistic),"\n")
				cat("p-value =",(ad.test(value.M2)$p.value),"\n\n")
			}
		}else{
			cat("Warning message: because of the negative PFT or the negative drift rate estimation,","\n",
			"the functions in the Goodness-of-Fit will not work because it is meaningless to do so.","\n")
		}
	}
}


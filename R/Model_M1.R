Model_M1 <-
function(Data, int.eta, int.seta, int.se, W, q, alpha, FIM, OIMA, OIMB, OIMC, t2, h.val, pseudoVal, MTTFval, quanval, pdfplotval, cdfplotval, ppplotval, qqplotval, kstestval, adtestval, dataname, optim.alg, env.in){

	Data <- as.matrix(Data)

	time <- Data[,1]

	m <- length(time); n <- length(Data[1,-1])

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

			stop('Please try different initial values (and/or the optimization algorithm)\n  to obtain more stable and more accurate estimates of the unknown parameters.')

		}
	}

	if(iter.M1 $ convergence == 0){

		Eta.M1 <- iter.M1 $ par[1]

		sEta.M1 <- abs( iter.M1 $ par[2] )

		sE.M1 <- abs( iter.M1 $ par[3] )

	}

	O <- sE.M1^2 * diag(m)

	S <- sEta.M1^2 * time %*% t(time) + O

	par.M1 <- c(Eta.M1, sEta.M1, sE.M1)

	assign('par.M1', par.M1, env.in)

	loglik.M1 <- - n * m * log( 2 * pi ) / 2 - n * log( det(S) ) / 2 - 

	            sum( diag( as.matrix( t( Data[,-1] - Eta.M1 * time ) ) %*% solve(S) %*% as.matrix( Data[,-1] - Eta.M1 * time ) ) ) / 2

	pdf_M1 <- function(t){

			sqrt( W^2 / ( 2 * pi * t^4 * sEta.M1^2 ) ) *

			exp( - ( W - Eta.M1 * t )^2 / ( 2 * t^2 * sEta.M1^2 ) )

			}

	cdf_M1 <- function(t){

		pnorm( ( Eta.M1 * t - W ) / ( sEta.M1 * t ) )

	}

	cdf_M1_INF <- pnorm( Eta.M1 / sEta.M1 )

	tq_M1 <- function(p1){

		t.percM1 <- function(t){

			pnorm( ( Eta.M1 * t - W ) / ( sEta.M1 * t ) ) - p1

		}
  
		uniroot( t.percM1, range(t2) )$root
	}

	if(FIM == 1){

		Cov.Matrix.FI.M1 <- solve( FI.M1( sEta.M1, O, S, time ) )

		par.FI.M1 <- par.FICI.M1(alpha, par.M1, Cov.Matrix.FI.M1, n )

		if(MTTFval==1 && Eta.M1>0){
			MTTF.FI.M1 <- MTTF.FICI.M1.D2( alpha, W, Eta.M1, sEta.M1, Cov.Matrix.FI.M1, n, h.val )
		}
		if(quanval==1 && cdf_M1_INF>max(q)) tq.FI.M1 <- try(tq.FICI.M1.D( q, alpha, W, Eta.M1, sEta.M1, Cov.Matrix.FI.M1, n , range(t2) ),T)

		if(cdfplotval==1) CDF.FI.M1 <- cdf.FICI.logit.M1.D( Par.M1 = par.M1, t=t2, Cov.Matrix.FI.M1, alpha=alpha, W=W, n )

	}

	if(OIMA == 1){

		A <- OI.A.M1( Data, Eta.M1, sEta.M1, O, S, time )

		Cov.Matrix.OI.M1.A <- solve(A)

		par.OI.M1.A <- par.OICI.M1( alpha, par.M1, Cov.Matrix.OI.M1.A )

		if(MTTFval==1 && Eta.M1>0){
			MTTF.OIA.M1 <- MTTF.OICI.M1.D2( alpha, W, Eta.M1, sEta.M1, Cov.Matrix.OI.M1.A, h.val )
		}

		if(quanval==1 && cdf_M1_INF>max(q)) tq.OIA.M1 <- try(tq.OICI.M1.D( q, alpha, W, Eta.M1, sEta.M1, Cov.Matrix.OI.M1.A, range(t2) ),T)

		if(cdfplotval==1) CDF.OIA.M1 <- cdf.OICI.logit.M1.D( Par.M1 = par.M1, t=t2, Cov.Matrix.OI.M1.A, alpha=alpha, W=W )

	}

	if(OIMB == 1){

		B <- OI.B.M1( Data, Eta.M1, sEta.M1, O, S, time )

		Cov.Matrix.OI.M1.B <- solve(B)

		par.OI.M1.B <- par.OICI.M1( alpha, par.M1, Cov.Matrix.OI.M1.B )

		if(MTTFval==1 && Eta.M1>0){
			MTTF.OIB.M1 <- MTTF.OICI.M1.D2( alpha, W, Eta.M1, sEta.M1, Cov.Matrix.OI.M1.B, h.val )
		}

		if(quanval==1 && cdf_M1_INF>max(q)) tq.OIB.M1 <- try(tq.OICI.M1.D( q, alpha, W, Eta.M1, sEta.M1, Cov.Matrix.OI.M1.B, range(t2) ),T)

		if(cdfplotval==1) CDF.OIB.M1 <- cdf.OICI.logit.M1.D( Par.M1 = par.M1, t=t2, Cov.Matrix.OI.M1.B, alpha=alpha, W=W )

	}

	if(OIMC == 1){

		A <- OI.A.M1( Data, Eta.M1, sEta.M1, O, S, time )

		B <- OI.B.M1( Data, Eta.M1, sEta.M1, O, S, time )

		Cov.Matrix.OI.M1.C <- solve(A) %*% B %*% solve(A)

		par.OI.M1.C <- par.OICI.M1( alpha, par.M1, Cov.Matrix.OI.M1.C )

		if(MTTFval==1 && Eta.M1>0){
			MTTF.OIC.M1 <- MTTF.OICI.M1.D2( alpha, W, Eta.M1, sEta.M1, Cov.Matrix.OI.M1.C, h.val )
		}

		if(quanval==1 && cdf_M1_INF>max(q)) tq.OIC.M1 <- try(tq.OICI.M1.D( q, alpha, W, Eta.M1, sEta.M1, Cov.Matrix.OI.M1.C, range(t2) ),T)

		if(cdfplotval==1) CDF.OIC.M1 <- cdf.OICI.logit.M1.D( Par.M1 = par.M1, t=t2, Cov.Matrix.OI.M1.C, alpha=alpha, W=W )

	}

	if(MTTFval==1 && Eta.M1>0){

		MTTF.approx.M1 <- MTTF.M1( par.M1[1], par.M1[2], W )[1]

		MTTF.exact.M1 <- MTTF.M1( par.M1[1], par.M1[2], W )[2]

	}

	if(quanval==1 && cdf_M1_INF>max(q)){

		fun_tp <- function(p1){

			t.percM1 <- function(t){

				pnorm( ( Eta.M1 * t - W ) / ( sEta.M1 * t ) ) - p1

			}
  
			tp <- uniroot( t.percM1, range(t2) )$root

		}

		tt <- try(apply( t(q), 2, fun_tp),T)

	}


	eta.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	seta.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	se.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	MTTF.approx.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	MTTF.exact.CI.M <- matrix(0,sum(FIM,OIMA,OIMB,OIMC),4)
	tq.CI.M <- array(0,c(sum(FIM,OIMA,OIMB,OIMC),4,length(q)))

	ind <- 1

	if(FIM==1){
		eta.CI.M[ind,] <- t(c(par.FI.M1$CI.Eta.M1,par.FI.M1$CI.ln.Eta.M1))
		seta.CI.M[ind,] <- t(c(par.FI.M1$CI.sEta.M1,par.FI.M1$CI.ln.sEta.M1))
		se.CI.M[ind,] <- t(c(par.FI.M1$CI.sE.M1,par.FI.M1$CI.ln.sE.M1))
		if(MTTFval==1 && Eta.M1>0) MTTF.approx.CI.M[ind,] <- t(c(MTTF.FI.M1$CI.MTTF.approx.M1, MTTF.FI.M1$CI.ln.MTTF.approx.M1))
		if(MTTFval==1 && Eta.M1>0) MTTF.exact.CI.M[ind,] <- t(c(MTTF.FI.M1$CI.MTTF.exact.M1, MTTF.FI.M1$CI.ln.MTTF.exact.M1))
		if(quanval==1 && cdf_M1_INF>max(q) && length(tq.FI.M1)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.FI.M1[i,2:5])
			}
		}
		ind <- ind + 1
	}

	if(OIMA==1){
		eta.CI.M[ind,] <- t(c(par.OI.M1.A$CI.Eta.M1,par.OI.M1.A$CI.ln.Eta.M1))
		seta.CI.M[ind,] <- t(c(par.OI.M1.A$CI.sEta.M1,par.OI.M1.A$CI.ln.sEta.M1))
		se.CI.M[ind,] <- t(c(par.OI.M1.A$CI.sE.M1,par.OI.M1.A$CI.ln.sE.M1))
		if(MTTFval==1 && Eta.M1>0) MTTF.approx.CI.M[ind,] <- t(c(MTTF.OIA.M1$CI.MTTF.approx.M1, MTTF.OIA.M1$CI.ln.MTTF.approx.M1))
		if(MTTFval==1 && Eta.M1>0) MTTF.exact.CI.M[ind,] <- t(c(MTTF.OIA.M1$CI.MTTF.exact.M1, MTTF.OIA.M1$CI.ln.MTTF.exact.M1))
		if(quanval==1 && cdf_M1_INF>max(q) && length(tq.OIA.M1)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.OIA.M1[i,2:5])
			}
		}
		ind <- ind + 1
	}

	if(OIMB==1){
		eta.CI.M[ind,] <- t(c(par.OI.M1.B$CI.Eta.M1,par.OI.M1.B$CI.ln.Eta.M1))
		seta.CI.M[ind,] <- t(c(par.OI.M1.B$CI.sEta.M1,par.OI.M1.B$CI.ln.sEta.M1))
		se.CI.M[ind,] <- t(c(par.OI.M1.B$CI.sE.M1,par.OI.M1.B$CI.ln.sE.M1))
		if(MTTFval==1 && Eta.M1>0) MTTF.approx.CI.M[ind,] <- t(c(MTTF.OIB.M1$CI.MTTF.approx.M1, MTTF.OIB.M1$CI.ln.MTTF.approx.M1))
		if(MTTFval==1 && Eta.M1>0) MTTF.exact.CI.M[ind,] <- t(c(MTTF.OIB.M1$CI.MTTF.exact.M1, MTTF.OIB.M1$CI.ln.MTTF.exact.M1))
		if(quanval==1 && cdf_M1_INF>max(q) && length(tq.OIB.M1)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.OIB.M1[i,2:5])
			}
		}
		ind <- ind + 1
	}

	if(OIMC==1){
		eta.CI.M[ind,] <- t(c(par.OI.M1.C$CI.Eta.M1,par.OI.M1.C$CI.ln.Eta.M1))
		seta.CI.M[ind,] <- t(c(par.OI.M1.C$CI.sEta.M1,par.OI.M1.C$CI.ln.sEta.M1))
		se.CI.M[ind,] <- t(c(par.OI.M1.C$CI.sE.M1,par.OI.M1.C$CI.ln.sE.M1))
		if(MTTFval==1 && Eta.M1>0) MTTF.approx.CI.M[ind,] <- t(c(MTTF.OIC.M1$CI.MTTF.approx.M1, MTTF.OIC.M1$CI.ln.MTTF.approx.M1))
		if(MTTFval==1 && Eta.M1>0) MTTF.exact.CI.M[ind,] <- t(c(MTTF.OIC.M1$CI.MTTF.exact.M1, MTTF.OIC.M1$CI.ln.MTTF.exact.M1))
		if(quanval==1 && cdf_M1_INF>max(q) && length(tq.OIC.M1)>1){
			for(i in 1:length(q)){
				tq.CI.M[ind,,i] <- t(tq.OIC.M1[i,2:5])
			}
		}
		ind <- ind + 1
	}


	eta.CI.M <- data.frame(eta.CI.M)
	seta.CI.M <- data.frame(seta.CI.M)
	se.CI.M <- data.frame(se.CI.M)
	if(MTTFval==1 && Eta.M1>0) MTTF.approx.CI.M <- data.frame(MTTF.approx.CI.M)
	if(MTTFval==1 && Eta.M1>0) MTTF.exact.CI.M <- data.frame(MTTF.exact.CI.M)
	if(quanval==1 && cdf_M1_INF>max(q)) tq.CI.M <- data.frame(tq.CI.M)

	row.name <- c("FIM", "OIM(Hessian matrix)", "OIM(score vector)",
			"OIM(robust matrix)")

	name.ind <- c(FIM,OIMA,OIMB,OIMC)
	row.names(eta.CI.M) <- row.name[which(name.ind==1)]
	row.names(seta.CI.M) <- row.name[which(name.ind==1)]
	row.names(se.CI.M) <- row.name[which(name.ind==1)]
	if(MTTFval==1 && Eta.M1>0) row.names(MTTF.approx.CI.M) <- row.name[which(name.ind==1)]
	if(MTTFval==1 && Eta.M1>0) row.names(MTTF.exact.CI.M) <- row.name[which(name.ind==1)]
	if(quanval==1 && cdf_M1_INF>max(q)) row.names(tq.CI.M) <- row.name[which(name.ind==1)]

	names(eta.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	names(seta.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	names(se.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	if(MTTFval==1 && Eta.M1>0) names(MTTF.approx.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	if(MTTFval==1 && Eta.M1>0) names(MTTF.exact.CI.M) <- c("LCI","UCI","LCI.ln","UCI.ln")
	if(quanval==1 && cdf_M1_INF>max(q)) names(tq.CI.M) <- rep(c("LCI","UCI","LCI.ln","UCI.ln"),length(q))

	cat(paste(rep("#",70),sep="",collapse=""),"\n")
	cat("     Model: Random effect + Measurement error    ","\n")
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

	cat("eta =",format(Eta.M1,digits=22),"\n")
	if(ind>1) print(eta.CI.M)
	cat("\n")

	cat("sigma_eta =",format(sEta.M1,digits=22),"\n")
	if(ind>1) print(seta.CI.M)
	cat("\n")

	cat("sigma_epsilon =",format(sE.M1,digits=22),"\n")
	if(ind>1) print(se.CI.M)
	cat("\n")

	cat("log-likelihood =",format(loglik.M1,digits=22),"\n")
      cat("Optimization Algorithm:",optim.alg,"\n")
	cat("\n")

	if(MTTFval==1){
		if(Eta.M1>0){
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				cat(paste(rep("-",60),sep="",collapse=""),"\n")
				cat(paste("        MTTF and ", 100*(1-alpha),"% confidence interval",sep=""),"\n")
				cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
			}else{
				cat(paste(rep("-",60),sep="",collapse=""),"\n")
				cat("        MTTF","\n")
				cat(paste(rep("-",60),sep="",collapse=""),"\n\n")
			}
			cat("MTTF(exact)=",MTTF.exact.M1,"\n")
			if(ind>1) print(MTTF.exact.CI.M)
			cat("\n")

			cat("MTTF(approx.)=",MTTF.approx.M1,"\n")
			if(ind>1) print(MTTF.approx.CI.M)
			cat("\n")
		}else{
			cat("\n\n\n")
			cat('Warning message: the estimated drift rate is negative, so it is meaningless to calculate the MTTF.','\n\n\n\n')
		}
	}

	if(quanval==1 && cdf_M1_INF>max(q)){
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
		cat("limit_{ t \ to infty } F_T(t) = ",cdf_M1_INF,"\n\n")
	}
	if(quanval==1 && cdf_M1_INF<max(q)){
		cat("limit_{ t \ to infty } F_T(t) = ",cdf_M1_INF,"\n\n")
	}


	if( (pdfplotval==1 || cdfplotval==1) && pseudoVal=="0" ){
		if(pdfplotval==1){
			x11()
			plot( t2, pdf_M1(t2), type="l", main="PDF estimation",
				 xlab="t", ylab=expression(f[T](t)), lwd=2, col=1 )
		}
		if(cdfplotval==1){
			x11()
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				plot( t2, apply(as.matrix(t2),1,cdf_M1), type="l", ylim=c(0,1), main='CDF and confidence interval with logit transformation',
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}else{
				plot( t2, apply(as.matrix(t2),1,cdf_M1), type="l", ylim=c(0,1), main="CDF estimation",
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}
			ind2 <- 2
			if(FIM==1){
				lines( t2, CDF.FI.M1$lcl.cdf.M1, lty = ind2, col = ind2 )
				lines( t2, CDF.FI.M1$ucl.cdf.M1, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMA==1){
				lines( t2, CDF.OIA.M1$lcl.cdf.M1, lty = ind2, col = ind2 )
				lines( t2, CDF.OIA.M1$ucl.cdf.M1, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMB==1){
				lines( t2, CDF.OIB.M1$lcl.cdf.M1, lty = ind2, col = ind2 )
				lines( t2, CDF.OIB.M1$ucl.cdf.M1, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMC==1){
				lines( t2, CDF.OIC.M1$lcl.cdf.M1, lty = ind2, col = ind2 )
				lines( t2, CDF.OIC.M1$ucl.cdf.M1, lty = ind2, col = ind2 )
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
			plot( t2, pdf_M1(t2), type="l", main="PDF estimation",
				 xlab="t", ylab=expression(f[T](t)), lwd=2, col=1 )
			if( all(PFT>0) && Eta.M1>0 ){
				points(sort(PFT), rep(0,length(PFT)), pch=1 )
				legend('topright', c("PDF estimation","Pseudo failure time"),
				col=c(1,1),lty=c(1,-1),lwd=c(2,-1),pch=c(-1,1),merge=TRUE,box.col="white",bg="white",inset = .01)
			}
		}
		if(cdfplotval==1){
			x11()
			if(FIM == 1 || OIMA == 1 || OIMB == 1 || OIMC == 1){
				plot( t2, apply(as.matrix(t2),1,cdf_M1), type="l", ylim=c(0,1), main='CDF and confidence interval with logit transformation',
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}else{
				plot( t2, apply(as.matrix(t2),1,cdf_M1), type="l", ylim=c(0,1), main="CDF estimation",
					 xlab="t", ylab=expression(F[T](t)), lwd=2, col=1 )
			}
			ind2 <- 2
			if(FIM==1){
				lines( t2, CDF.FI.M1$lcl.cdf.M1, lty = ind2, col = ind2 )
				lines( t2, CDF.FI.M1$ucl.cdf.M1, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMA==1){
				lines( t2, CDF.OIA.M1$lcl.cdf.M1, lty = ind2, col = ind2 )
				lines( t2, CDF.OIA.M1$ucl.cdf.M1, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMB==1){
				lines( t2, CDF.OIB.M1$lcl.cdf.M1, lty = ind2, col = ind2 )
				lines( t2, CDF.OIB.M1$ucl.cdf.M1, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			if(OIMC==1){
				lines( t2, CDF.OIC.M1$lcl.cdf.M1, lty = ind2, col = ind2 )
				lines( t2, CDF.OIC.M1$ucl.cdf.M1, lty = ind2, col = ind2 )
				ind2 <- ind2 + 1
			}
			text1 <- c(
			paste("FIM ", 100*(1-alpha),"% CI",sep=""),paste("OIM (Hessian matrix) ", 100*(1-alpha),"% CI",sep=""),
			paste("OIM (score vector) ", 100*(1-alpha),"% CI",sep=""),
			paste("OIM (robust matrix) ", 100*(1-alpha),"% CI",sep=""))

			ind3 <- which(c(FIM,OIMA,OIMB,OIMC)==1)

			if( all(PFT>0) && Eta.M1>0 ){
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
		if(all(PFT>0) && Eta.M1>0){
			if(ppplotval==1){
				x11()
				plot((1:length(PFT)-0.5)/length(PFT),apply(as.matrix(sort(PFT)),1,cdf_M1),main="Probability-probability plot",xlab="Theoretical cumulative probabilities",ylab="Sample cumulative probabilities",xlim=c(0,1),ylim=c(0,1))
				abline(0,1)
			}
			if(qqplotval==1 && cdf_M1_INF>max(q)){
				qx <- try(apply(as.matrix((1:length(PFT)-0.5)/length(PFT)),1,tq_M1),T)
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
				value.M1 <- NULL
				for(i in 1:n) value.M1 <- c(value.M1, cdf_M1(sort(PFT)[i]))
				cat("Kolmogorov-Smirnov test :","\n")
				cat("statistic =",(STATISTIC <- max(abs(value.M1 - (1:n - 1)/(n)))),"\n")
				cat("p-value =",(1 - .C("pkolmogorov2x", p = as.double(STATISTIC), as.integer(n), PACKAGE = "stats")$p),"\n\n")
			}
			if(adtestval==1){
				value.M1 <- NULL
				for(i in 1:n) value.M1 <- c(value.M1, cdf_M1(sort(PFT)[i]))
				cat("Anderson-Darling test","\n")
				cat("statistic =",(STATISTIC.M1 <- ad.test(value.M1)$statistic),"\n")
				cat("p-value =",(ad.test(value.M1)$p.value),"\n\n")
			}
		}else{
			cat("Warning message: because of the negative PFT or the negative drift rate estimation,","\n",
			"the functions in the Goodness-of-Fit will not work because it is meaningless to do so.","\n")
		}
	}
}


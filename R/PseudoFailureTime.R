PseudoFailureTime <-
function(Data,W,pftname){
	time <- Data[,1]
	fit <- function(y,...)	lm(y~0+time,data=Data)$coef
	til.eta <- apply( Data[,-1], 2, fit )
	x <- W / til.eta
	cat(paste(rep("-",75),sep="",collapse=""),"\n")
	cat("                Pseudo Failure Time Estimation","\n")
	cat(paste(rep("-",75),sep="",collapse=""),"\n\n")
	cat('Threshold =',W,'\n\n')
	print(x)
	cat('\n')
	if(any(x<0)) cat('Warning message: the negative PFT estimation does not make sense in practice.','\n')
	assign(paste('PFT.',pftname,sep=''), x, envir=.GlobalEnv)
	cat("\n")
}


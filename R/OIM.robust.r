OIM.robust <- function(Data, Eta, sEta, sB, sE){
	Time <- Data[,1]
	if( sEta>=0 & sB>=0 & sE>=0 && all(Time>0) ){
		A <- OIM.hessian(Data, Eta, sEta, sB, sE)
		B <- OIM.score(Data, Eta, sEta, sB, sE)
		C <- solve(A) %*% B %*% solve(A)
		return(C)
	}else{
		cat("The parameter and measuring time should not be negative and should be a positive vector, respectively!","\n")
	}
}
par.FICI.M3 <-
function(alpha, Par.M3, Cov.Mat.M3, n, ...){

    Eta.M3 <- Par.M3[1]

    sEta.M3 <- Par.M3[2]

    sB.M3 <- Par.M3[3]

    CI.lower.eta <- Eta.M3 - qnorm( 1 - alpha / 2 ) * sqrt( Cov.Mat.M3[1,1] / n )

    CI.upper.eta <- Eta.M3 + qnorm( 1 - alpha / 2 ) * sqrt( Cov.Mat.M3[1,1] / n )

    if(Eta.M3>0){
	    CI.ln.lower.eta <- exp( log(Eta.M3) - qnorm( 1 - alpha / 2 ) * sqrt( Eta.M3^(-2) * Cov.Mat.M3[1,1] / n ) )
	    CI.ln.upper.eta <- exp( log(Eta.M3) + qnorm( 1 - alpha / 2 ) * sqrt( Eta.M3^(-2) * Cov.Mat.M3[1,1] / n ) )
    }else{
	    CI.ln.lower.eta <- NA
	    CI.ln.upper.eta <- NA
    }

    var.s.eta <- ( 2 * sEta.M3 )^(-2) * Cov.Mat.M3[2,2]

    CI.lower.s.eta <- sEta.M3 - qnorm( 1 - alpha / 2 ) * sqrt( var.s.eta / n )

    CI.upper.s.eta <- sEta.M3 + qnorm( 1 - alpha / 2 ) * sqrt( var.s.eta / n )

    CI.ln.lower.s.eta <- exp( log(sEta.M3) - qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sEta.M3^2 )^(-2) * Cov.Mat.M3[2,2] / n ) )

    CI.ln.upper.s.eta <- exp( log(sEta.M3) + qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sEta.M3^2 )^(-2) * Cov.Mat.M3[2,2] / n ) )

    var.sB <- (2*sB.M3)^(-2) * Cov.Mat.M3[3,3]

    CI.lower.sB <- sB.M3 - qnorm( 1 - alpha / 2 ) * sqrt( var.sB / n )

    CI.upper.sB <- sB.M3 + qnorm( 1 - alpha / 2 ) * sqrt( var.sB / n )

    CI.ln.lower.sB <- exp( log(sB.M3) - qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sB.M3^2 )^(-2) * Cov.Mat.M3[3,3] / n ) )

    CI.ln.upper.sB <- exp( log(sB.M3) + qnorm( 1 - alpha / 2 ) * sqrt( ( 2 * sB.M3^2 )^(-2) * Cov.Mat.M3[3,3] / n ) )

    list( CI.Eta.M3 = c(CI.lower.eta, CI.upper.eta), 

         CI.ln.Eta.M3 = c(CI.ln.lower.eta, CI.ln.upper.eta), 

         CI.sEta.M3 = c(CI.lower.s.eta, CI.upper.s.eta), 

         CI.ln.sEta.M3 = c(CI.ln.lower.s.eta, CI.ln.upper.s.eta), 

         CI.sB.M3 = c(CI.lower.sB, CI.upper.sB), 

         CI.ln.sB.M3 = c(CI.ln.lower.sB, CI.ln.upper.sB) )

}


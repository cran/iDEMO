FI.M4 <-
function( sE.M4, Time, m , ...){

    FI.M4 <- matrix( 0, 2, 2 )

    tOt <- c( t(Time) %*% diag(m) %*% Time / sE.M4^2 )

    FI.M4[1,1] <- tOt

    FI.M4[2,2] <- m / ( 2 * sE.M4^4 )

    return(FI.M4)

}


FI.M2 <-
function( sB.M2, Qinv, Time, m, ...){

    tQt <- c( t(Time) %*% Qinv %*% Time )

    FI.M2 <- matrix( 0, 2, 2 )

    FI.M2[1,1] <- tQt / sB.M2^2

    FI.M2[2,2] <- m / ( 2 * sB.M2^4 )

    return(FI.M2)

}


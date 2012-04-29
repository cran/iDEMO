FI.M1 <-
function( sEta.M1, O, S, Time,...){

    FI.M1 <- matrix( 0, 3, 3 )

    tSt <- c( t(Time) %*% solve(S) %*% Time )

    tSSt <- c( t(Time) %*% solve(S) %*% solve(S) %*% Time )

    SS <- solve(S) %*% solve(S)

    FI.M1[1,1] <- tSt

    FI.M1[2,2] <- tSt^2 / 2

    FI.M1[2,3] <- FI.M1[3,2] <- tSSt / 2

    FI.M1[3,3] <- sum( diag(SS) ) / 2

    return(FI.M1)

}


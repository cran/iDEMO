FI.M3 <-
function( Q, O, S, Time,...){

    FI.M3 <- matrix( 0, 3, 3 )

    tSt <- c( t(Time) %*% solve(S) %*% Time )

    tSQSt <- c( t(Time) %*% solve(S) %*% Q %*% solve(S) %*% Time )

    QSQS <- Q %*% solve(S) %*% Q %*% solve(S)


    FI.M3[1,1] <- tSt

    FI.M3[2,2] <- (tSt)^2 / 2

    FI.M3[2,3] <- FI.M3[3,2] <- tSQSt / 2

    FI.M3[3,3] <- sum( diag(QSQS) ) / 2

    return(FI.M3)

}


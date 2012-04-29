FI.M5 <-
function( Q, O, Time,...){

    tOt <- c( t(Time) %*% solve(O) %*% Time )

    QOQO <- Q %*% solve(O) %*% Q %*% solve(O)

    QOO <- Q %*% solve(O) %*% solve(O)

    OO <- solve(O) %*% solve(O)

    FI.M5 <- matrix( 0, 3, 3 )

    FI.M5[1,1] <- tOt

    FI.M5[2,2] <- sum( diag(QOQO) ) / 2

    FI.M5[2,3] <- FI.M5[3,2] <- sum( diag(QOO) ) / 2

    FI.M5[3,3] <- sum( diag(OO) ) / 2

    return(FI.M5)

}


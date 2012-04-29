FI.M0 <-
function( sEta.M0, Q, O, S, Time,...){

    tSt <- c( t(Time) %*% solve(S) %*% Time )

    tOt <- c( t(Time) %*% solve(O) %*% Time )

    one_sEta2tOt <- c( 1 + sEta.M0^2 * tOt )

    tOQOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Time )

    tOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Time )

    QOQO <- Q %*% solve(O) %*% Q %*% solve(O)

    tOQOQOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Q %*% solve(O) %*% Time )

    QOO <- Q %*% solve(O) %*% solve(O)

    tOOQOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Q %*% solve(O) %*% Time )

    OO <- solve(O) %*% solve(O)

    tOOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% solve(O) %*% Time )

    FI.M0 <- matrix( 0, 4, 4 )

    FI.M0[1,1] <- tSt

    FI.M0[2,2] <- (tOt)^2 / ( 2 * one_sEta2tOt^2 )

    FI.M0[2,3] <- FI.M0[3,2] <- tOQOt / ( 2 * one_sEta2tOt^2 )

    FI.M0[2,4] <- FI.M0[4,2] <- tOOt / ( 2 * one_sEta2tOt^2 )

    FI.M0[3,3] <- sum( diag(QOQO) ) / 2 + sEta.M0^4 * tOQOt^2 / ( 2 * one_sEta2tOt^2 ) - 

               sEta.M0^2 * tOQOQOt / one_sEta2tOt

    FI.M0[3,4] <- FI.M0[4,3] <- sum( diag(QOO) ) / 2 + sEta.M0^4 * tOQOt * tOOt / ( 2 * one_sEta2tOt^2 ) -

               sEta.M0^2 * tOOQOt / one_sEta2tOt

    FI.M0[4,4] <- sum( diag(OO) ) / 2 + sEta.M0^4 * tOOt^2 / ( 2 * one_sEta2tOt^2 ) - 

               sEta.M0^2 * tOOOt / one_sEta2tOt


    return(FI.M0)

}


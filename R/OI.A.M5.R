OI.A.M5 <-
function( Data, Eta.M5, Q, O, Time,...){

    tOt <- c( t(Time) %*% solve(O) %*% Time )

    tOQOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Time )

    tOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Time )

    OQOQ <- solve(O) %*% Q %*% solve(O) %*% Q

    OQOQO <- solve(O) %*% Q %*% solve(O) %*% Q %*% solve(O)

    QOO <- Q %*% solve(O) %*% solve(O)

    OOQO <- solve(O) %*% solve(O) %*% Q %*% solve(O)

    OQOO <- solve(O) %*% Q %*% solve(O) %*% solve(O)

    OO <- solve(O) %*% solve(O)

    OOO <- solve(O) %*% solve(O) %*% solve(O)


    A.cal <- function( vec.y, ... ){

        tOQOy <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% vec.y )

        tOOy <- c( t(Time) %*% solve(O) %*% solve(O) %*% vec.y )

        yEtat <- c( vec.y - Eta.M5 * Time )

        M.A.temp <- matrix( 0, 3, 3 )

        M.A.temp[1,1] <- - tOt

        M.A.temp[1,2] <- M.A.temp[2,1] <- - tOQOy + Eta.M5 * tOQOt

        M.A.temp[1,3] <- M.A.temp[3,1] <- - tOOy + Eta.M5 * tOOt

        M.A.temp[2,2] <- sum( diag( OQOQ ) ) / 2 - t(yEtat) %*% OQOQO %*% yEtat

        M.A.temp[2,3] <- M.A.temp[3,2] <- sum( diag( QOO ) ) / 2 - t(yEtat) %*% ( OOQO + OQOO ) %*% yEtat / 2

        M.A.temp[3,3] <- sum( diag( OO ) ) / 2 - t(yEtat) %*% OOO %*% yEtat
 
        M.A.temp

    }


    A.cal.2 <- apply( Data[,-1], 2, A.cal )

    OI.A.M5 <- (-1) * matrix( apply( A.cal.2, 1, sum ), 3, 3 )

    return(OI.A.M5)

}


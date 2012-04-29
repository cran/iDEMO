OI.A.M0 <-
function( Data, Eta.M0, sEta.M0, Q, O, S, Time,...){

    tSt <- c( t(Time) %*% solve(S) %*% Time )

    tOt <- c( t(Time) %*% solve(O) %*% Time )

    one_sEta2tOt <- c( 1 + sEta.M0^2 * tOt )

    tOQOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Time )

    tOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Time )

    OttO <- solve(O) %*% Time %*% t(Time) %*% solve(O)

    OQOttO <- solve(O) %*% Q %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

    OttOQO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% Q %*% solve(O)

    OOttO <- solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

    OttOO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O)

    OQOQ <- solve(O) %*% Q %*% solve(O) %*% Q

    tOQOQOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Q %*% solve(O) %*% Time )

    OQOQO <- solve(O) %*% Q %*% solve(O) %*% Q %*% solve(O)

    OQOQOttO <- solve(O) %*% Q %*% solve(O) %*% Q %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

    OQOttOQO <- solve(O) %*% Q %*% solve(O) %*% Time %*% t(Time) %*% solve(O) %*% Q %*% solve(O)

    OttOQOQO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Q %*% solve(O)

    QOO <- Q %*% solve(O) %*% solve(O)

    tOQOOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% solve(O) %*% Time )

    tOOQOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Q %*% solve(O) %*% Time )

    OOQO <- solve(O) %*% solve(O) %*% Q %*% solve(O)

    OQOO <- solve(O) %*% Q %*% solve(O) %*% solve(O)

    OOQOttO <- solve(O) %*% solve(O) %*% Q %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

    OQOOttO <- solve(O) %*% Q %*% solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

    OQOttOO <- solve(O) %*% Q %*% solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O)

    OOttOQO <- solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O) %*% Q %*% solve(O)

    OttOOQO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O) %*% Q %*% solve(O)

    OttOQOO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% Q %*% solve(O) %*% solve(O)

    OO <- solve(O) %*% solve(O)

    tOOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% solve(O) %*% Time )

    OOO <- solve(O) %*% solve(O) %*% solve(O)

    OOOttO <- solve(O) %*% solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

    OOttOO <- solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O)

    OttOOO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O) %*% solve(O)
    

    A.cal <- function( vec.y, ... ){

        tOy <- c( t(Time) %*% solve(O) %*% vec.y )

        tOQOy <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% vec.y )

        tOOy <- c( t(Time) %*% solve(O) %*% solve(O) %*% vec.y )

        yEtat <- c( vec.y - Eta.M0 * Time )

        M.A.temp <- matrix( 0, 4, 4 )

        M.A.temp[1,1] <- - tSt

        M.A.temp[1,2] <- M.A.temp[2,1] <- (- tOt * tOy + Eta.M0 * tOt^2 ) / one_sEta2tOt^2

        M.A.temp[1,3] <- M.A.temp[3,1] <- - tOQOy + Eta.M0 * tOQOt - sEta.M0^4 * tOQOt / one_sEta2tOt^2 * ( tOt * tOy - Eta.M0 * tOt^2 ) +

                                            sEta.M0^2 / one_sEta2tOt * ( tOQOt * tOy + tOt * tOQOy - 2 * Eta.M0 * tOt * tOQOt )

        M.A.temp[1,4] <- M.A.temp[4,1] <- - tOOy + Eta.M0 * tOOt - sEta.M0^4 * tOOt / one_sEta2tOt^2 * ( tOt * tOy - Eta.M0 * tOt^2 ) +

                                            sEta.M0^2 / one_sEta2tOt * ( tOOt * tOy + tOt * tOOy - 2 * Eta.M0 * tOt * tOOt )

        M.A.temp[2,2] <- tOt^2 / ( 2 * one_sEta2tOt^2 ) - tOt / one_sEta2tOt^3 * t(yEtat) %*% OttO %*% yEtat

        M.A.temp[2,3] <- M.A.temp[3,2] <- tOQOt / ( 2 * one_sEta2tOt^2 ) - t(yEtat) %*% ( OQOttO + OttOQO ) %*% yEtat / ( 2 * one_sEta2tOt^2 ) +

                                            sEta.M0^2 * tOQOt / one_sEta2tOt^3 * t(yEtat) %*% OttO %*% yEtat

        M.A.temp[2,4] <- M.A.temp[4,2] <- tOOt / ( 2 * one_sEta2tOt^2 ) - t(yEtat) %*% ( OOttO + OttOO ) %*% yEtat / ( 2 * one_sEta2tOt^2 ) +

                                            sEta.M0^2 * tOOt / one_sEta2tOt^3 * t(yEtat) %*% OttO %*% yEtat

        M.A.temp[3,3] <- sum( diag( OQOQ ) ) / 2 - sEta.M0^2 * tOQOQOt / one_sEta2tOt + sEta.M0^4 * tOQOt^2 / ( 2 * one_sEta2tOt^2 ) -

                          t(yEtat) %*% OQOQO %*% yEtat + sEta.M0^2 / one_sEta2tOt * t(yEtat) %*% ( OQOQOttO + OQOttOQO + OttOQOQO ) %*% yEtat -

                          sEta.M0^4 * tOQOt / one_sEta2tOt^2 * t(yEtat) %*% ( OQOttO + OttOQO ) %*% yEtat -
    
                          sEta.M0^4 * tOQOQOt / one_sEta2tOt^2 * t(yEtat) %*% OttO %*% yEtat + sEta.M0^6 * tOQOt^2 / one_sEta2tOt^3 *

                          t(yEtat) %*% OttO %*% yEtat

        M.A.temp[3,4] <- M.A.temp[4,3] <- sum( diag( QOO ) ) / 2 + sEta.M0^4 * tOOt * tOQOt / ( 2 * one_sEta2tOt^2 ) - sEta.M0^2 / ( 2 * one_sEta2tOt ) *

                                            ( tOQOOt + tOOQOt ) - 0.5 * t(yEtat) %*% ( OOQO + OQOO ) %*% yEtat -

                                            sEta.M0^4 * tOOt / ( 2 * one_sEta2tOt^2 ) * t(yEtat) %*% ( OQOttO + OttOQO ) %*% yEtat +

                                            sEta.M0^2 / ( 2 *one_sEta2tOt ) * t(yEtat) %*% ( OOQOttO + OQOOttO + OQOttOO + OOttOQO + OttOOQO + OttOQOO ) %*% yEtat +

                                            ( sEta.M0^6 * tOOt * tOQOt / one_sEta2tOt^3 - sEta.M0^4 * ( tOOQOt + tOQOOt ) / ( 2 *one_sEta2tOt^2 ) ) *

                                            t(yEtat) %*% OttO %*% yEtat -

                                            sEta.M0^4 * tOQOt / ( 2 * one_sEta2tOt^2 ) * t(yEtat) %*% ( OOttO + OttOO ) %*% yEtat

        M.A.temp[4,4] <- sum( diag( OO ) ) / 2  + sEta.M0^4 * tOOt^2 / ( 2 * one_sEta2tOt^2 ) - sEta.M0^2 * tOOOt / one_sEta2tOt -

                    t(yEtat) %*% OOO %*% yEtat - sEta.M0^4 * tOOt / one_sEta2tOt^2 * t(yEtat) %*% ( OOttO + OttOO ) %*% yEtat +   

                    sEta.M0^2 / one_sEta2tOt * t(yEtat) %*% ( OOOttO + OOttOO + OttOOO ) %*% yEtat +

                    sEta.M0^6 * tOOt^2 / one_sEta2tOt^3 * t(yEtat) %*% OttO %*% yEtat -

                    sEta.M0^4 * tOOOt / one_sEta2tOt^2 * t(yEtat) %*% OttO %*% yEtat
 
        M.A.temp

    }


    A.cal.2 <- apply( Data[,-1], 2, A.cal )

    OI.A.M0 <- (-1) * matrix( apply( A.cal.2, 1, sum ), 4, 4 )

    return(OI.A.M0)

}


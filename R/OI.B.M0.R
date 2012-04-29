OI.B.M0 <-
function( Data, Eta.M0, sEta.M0, Q, O, S, Time,...){

    tSt <- c( t(Time) %*% solve(S) %*% Time )

    tOt <- c( t(Time) %*% solve(O) %*% Time )

    one_sEta2tOt <- c( 1 + sEta.M0^2 * tOt )

    OttO <- solve(O) %*% Time %*% t(Time) %*% solve(O)

    OQ <- solve(O) %*% Q

    tOQOt <- c( t(Time) %*% solve(O) %*% Q %*% solve(O) %*% Time )

    OQO <- solve(O) %*% Q %*% solve(O)

    OQOttO <- solve(O) %*% Q %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

    OttOQO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% Q %*% solve(O)

    tOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Time )

    OO <- solve(O) %*% solve(O)

    OOttO <- solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

    OttOO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O)


    uu <- function( vec.y, ... ){

        u <- matrix(0,4,1)

        tSy <- t(Time) %*% solve(S) %*% vec.y

        yEtat <- ( vec.y - Eta.M0 * Time )

        u[1,1] <- tSy - Eta.M0 * tSt

        u[2,1] <- - tOt / ( 2 * one_sEta2tOt ) +  t(yEtat) %*% OttO %*% yEtat / 

                ( 2 * one_sEta2tOt^2 )

        u[3,1] <- - sum( diag(OQ) ) / 2  +  sEta.M0^2 * tOQOt / ( 2 * one_sEta2tOt )  +

                t(yEtat) %*% OQO %*% yEtat / 2  -  sEta.M0^2 / ( 2 * one_sEta2tOt ) * 

                t(yEtat) %*% ( OQOttO + OttOQO ) %*% yEtat  +

                sEta.M0^4 * tOQOt / ( 2 * one_sEta2tOt^2 ) * 

                t(yEtat) %*% OttO %*% yEtat

        u[4,1] <- - sum( diag( solve(O) ) ) / 2  +  sEta.M0^2 * tOOt / ( 2 * one_sEta2tOt )  +

                t(yEtat) %*% OO %*% yEtat / 2  -  sEta.M0^2 / ( 2 * one_sEta2tOt ) * 

                t(yEtat) %*% ( OOttO + OttOO ) %*% yEtat  +

                sEta.M0^4 * tOOt / ( 2 * one_sEta2tOt^2 ) * 

                t(yEtat) %*% OttO %*% yEtat

        u %*% t(u)
        
    }

    uu.cal <- apply( Data[,-1], 2, uu )

    OI.B.M0 <- matrix( apply( uu.cal, 1, sum ), 4, 4, byrow=TRUE )

    return(OI.B.M0)


}


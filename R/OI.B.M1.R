OI.B.M1 <-
function( Data, Eta.M1, sEta.M1, O, S, Time,...){

    tOt <- t(Time) %*% solve(O) %*% Time 

    one_sEta2tOt <- c( 1 + sEta.M1^2 * tOt )

    tSt <- c( t(Time) %*% solve(S) %*% Time )

    OttO <- solve(O) %*% Time %*% t(Time) %*% solve(O)    

    tOOt <- c( t(Time) %*% solve(O) %*% solve(O) %*% Time )

    OO <- solve(O) %*% solve(O)

    trace.O.inv <- sum( diag( solve(O) ) )

    OOttO <- solve(O) %*% solve(O) %*% Time %*% t(Time) %*% solve(O)

    OttOO <- solve(O) %*% Time %*% t(Time) %*% solve(O) %*% solve(O)

    uu <- function( vec.y, ... ){

        u <- matrix(0,3,1)

        tSy <- t(Time) %*% solve(S) %*% vec.y

        yEtat <- ( vec.y - Eta.M1 * Time )

        u[1,1] <- tSy - Eta.M1 * tSt

        u[2,1] <- - tOt / ( 2 * one_sEta2tOt ) +  t(yEtat) %*% OttO %*% yEtat / 

                ( 2 * one_sEta2tOt^2 )

        u[3,1] <- - trace.O.inv / 2  +  sEta.M1^2 * tOOt / ( 2 * one_sEta2tOt )  +

                t(yEtat) %*% OO %*% yEtat / 2  -  sEta.M1^2 / ( 2 * one_sEta2tOt ) * 

                t(yEtat) %*% ( OOttO + OttOO ) %*% yEtat  +

                sEta.M1^4 * tOOt / ( 2 * one_sEta2tOt^2 ) * 

                t(yEtat) %*% OttO %*% yEtat

        u %*% t(u)
        
    }

    uu.cal <- apply( Data[,-1], 2, uu )

    OI.B.M1 <- matrix( apply( uu.cal, 1, sum ), 3, 3, byrow=TRUE )

    return(OI.B.M1)

}


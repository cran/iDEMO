OI.B.M3 <-
function( Data, Eta.M3, sEta.M3, Q, O.inv, S, Time,...){

    tSt <- c( t(Time) %*% solve(S) %*% Time )
    
    tOt <- t(Time) %*% O.inv %*% Time

    one_sEta2tOt <- c( 1 + sEta.M3^2 * tOt )

    OttO <- O.inv %*% Time %*% t(Time) %*% O.inv
    
    OQ <- O.inv %*% Q

    tOQOt <- c( t(Time) %*% O.inv %*% Q %*% O.inv %*% Time )

    OQO <- O.inv %*% Q %*% O.inv

    OQOttO <- O.inv %*% Q %*% O.inv %*% Time %*% t(Time) %*% O.inv

    OttOQO <- O.inv %*% Time %*% t(Time) %*% O.inv %*% Q %*% O.inv




    uu <- function( vec.y, ... ){

        u <- matrix(0,3,1)

        tSy <- t(Time) %*% solve(S) %*% vec.y

        yEtat <- ( vec.y - Eta.M3 * Time )

        u[1,1] <- tSy - Eta.M3 * tSt

        u[2,1] <- - tOt / ( 2 * one_sEta2tOt ) +  t(yEtat) %*% OttO %*% yEtat / 

                ( 2 * one_sEta2tOt^2 )

        u[3,1] <- - sum( diag(OQ) ) / 2  +  sEta.M3^2 * tOQOt / ( 2 * one_sEta2tOt )  +

                t(yEtat) %*% OQO %*% yEtat / 2  -  sEta.M3^2 / ( 2 * one_sEta2tOt ) * 

                t(yEtat) %*% ( OQOttO + OttOQO ) %*% yEtat  +

                sEta.M3^4 * tOQOt / ( 2 * one_sEta2tOt^2 ) * 

                t(yEtat) %*% OttO %*% yEtat


        u %*% t(u)
        
    }

    uu.cal <- apply( Data[,-1], 2, uu )

    OI.B.M3 <- matrix( apply( uu.cal, 1, sum ), 3, 3, byrow=TRUE )

    return(OI.B.M3)

}


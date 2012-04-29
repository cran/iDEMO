OI.B.M2 <-
function( Data, Eta.M2, sB.M2, O, S, Qinv, Time, m, ...){

    tSt <- t(Time) %*% ( Qinv / sB.M2^2 ) %*% Time

    OQ <- diag(m) / sB.M2^2

    OQO <- Qinv / sB.M2^4

    uu <- function( vec.y, ... ){

        u <- matrix(0,2,1)

        tSy <- t(Time) %*% ( Qinv / sB.M2^2 ) %*% vec.y

        yEtat <- ( vec.y - Eta.M2 * Time )

        u[1,1] <- tSy - Eta.M2 * tSt

        u[2,1] <- - sum( diag(OQ) ) / 2  +  t(yEtat) %*% OQO %*% yEtat / 2

        u %*% t(u)
        
    }

    uu.cal <- apply( Data[,-1], 2, uu )

    OI.B.M2 <- matrix( apply( uu.cal, 1, sum ), 2, 2, byrow=TRUE )

    return(OI.B.M2)

}


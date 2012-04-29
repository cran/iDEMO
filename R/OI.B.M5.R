OI.B.M5 <-
function( Data, Eta.M5, Q, O, Time,...){

    tOt <- c( t(Time) %*% solve(O) %*% Time )

    OQ <- solve(O) %*% Q

    OQO <- solve(O) %*% Q %*% solve(O)

    OO <- solve(O) %*% solve(O)


    uu <- function( vec.y, ... ){

        u <- matrix(0,3,1)

        tOy <- t(Time) %*% solve(O) %*% vec.y

        yEtat <- ( vec.y - Eta.M5 * Time )

        u[1,1] <- tOy - Eta.M5 * tOt

        u[2,1] <- - sum( diag(OQ) ) / 2 + t(yEtat) %*% OQO %*% yEtat / 2

        u[3,1] <- - sum( diag( solve(O) ) ) / 2 + t(yEtat) %*% OO %*% yEtat / 2

        u %*% t(u)
        
    }

    uu.cal <- apply( Data[,-1], 2, uu )

    OI.B.M5 <- matrix( apply( uu.cal, 1, sum ), 3, 3, byrow=TRUE )

    return(OI.B.M5)

}


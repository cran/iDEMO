OI.A.M2 <-
function( Data, Eta.M2, sB.M2, S, Qinv, Time, m, ...){

    tQt <- c( t(Time) %*% Qinv %*% Time )

    SttS <- solve(S) %*% Time %*% t(Time) %*% solve(S)

    SSttS <- solve(S) %*% solve(S) %*% Time %*% t(Time) %*% solve(S)

    SttSS <- solve(S) %*% Time %*% t(Time) %*% solve(S) %*% solve(S)

    SSS <- solve(S) %*% solve(S) %*% solve(S)
    
    A.cal <- function( vec.y, ... ){

        tQy <- c( t(Time) %*% Qinv %*% vec.y )

        yEtat <- c( vec.y - Eta.M2 * Time )

        M.A.temp <- matrix( 0, 2, 2 )

        M.A.temp[1,1] <- - tQt / sB.M2^2

        M.A.temp[1,2] <- M.A.temp[2,1] <- ( - tQy + Eta.M2 * tQt ) / sB.M2^4

        M.A.temp[2,2] <- m / ( 2 * sB.M2^4 ) - t(yEtat) %*% Qinv %*% yEtat / sB.M2^6
 
        M.A.temp

    }


    A.cal.2 <- apply( Data[,-1], 2, A.cal )

    OI.A.M2 <- (-1) * matrix( apply( A.cal.2, 1, sum ), 2, 2 )

    return(OI.A.M2)

}


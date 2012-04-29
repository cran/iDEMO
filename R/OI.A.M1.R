OI.A.M1 <-
function( Data, Eta.M1, sEta.M1, O, S, Time,...){

    tSt <- c( t(Time) %*% solve(S) %*% Time )

    tSSt <- c( t(Time) %*% solve(S) %*% solve(S) %*% Time )

    SS <- solve(S) %*% solve(S)

    SttS <- solve(S) %*% Time %*% t(Time) %*% solve(S)

    SSttS <- solve(S) %*% solve(S) %*% Time %*% t(Time) %*% solve(S)

    SttSS <- solve(S) %*% Time %*% t(Time) %*% solve(S) %*% solve(S)

    SSS <- solve(S) %*% solve(S) %*% solve(S)
    
    A.cal <- function( vec.y, ... ){

        tSy <- c( t(Time) %*% solve(S) %*% vec.y )

        tSSy <- c( t(Time) %*% solve(S) %*% solve(S) %*% vec.y )

        yEtat <- c( vec.y - Eta.M1 * Time )

        M.A.temp <- matrix( 0, 3, 3 )

        M.A.temp[1,1] <- - tSt

        M.A.temp[1,2] <- M.A.temp[2,1] <- - tSt * tSy + Eta.M1 * tSt^2

        M.A.temp[1,3] <- M.A.temp[3,1] <- - tSSy + Eta.M1 * tSSt

        M.A.temp[2,2] <- tSt^2 / 2 - tSt * t(yEtat) %*% SttS %*% yEtat

        M.A.temp[2,3] <- M.A.temp[3,2] <- tSSt / 2 - t(yEtat) %*% ( SSttS + SttSS ) %*% yEtat / 2

        M.A.temp[3,3] <- sum( diag( SS ) ) / 2 - t(yEtat) %*% SSS %*% yEtat
 
        M.A.temp

    }


    A.cal.2 <- apply( Data[,-1], 2, A.cal )

    OI.A.M1 <- (-1) * matrix( apply( A.cal.2, 1, sum ), 3, 3 )

    return(OI.A.M1)

}


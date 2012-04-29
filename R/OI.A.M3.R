OI.A.M3 <-
function( Data, Eta.M3, Q, O, S, Time,...){

    tSt <- c( t(Time) %*% solve(S) %*% Time )

    tSQSt <- c( t(Time) %*% solve(S) %*% Q %*% solve(S) %*% Time )

    SttS <- solve(S) %*% Time %*% t(Time) %*% solve(S)

    SQSttS <- solve(S) %*% Q %*% solve(S) %*% Time %*% t(Time) %*% solve(S)

    SttSQS <- solve(S) %*% Time %*% t(Time) %*% solve(S) %*% Q %*% solve(S)

    SQSQ <- solve(S) %*% Q %*% solve(S) %*% Q

    SQSQS <- solve(S) %*% Q %*% solve(S) %*% Q %*% solve(S)
    
    A.cal <- function( vec.y, ... ){

        tSy <- c( t(Time) %*% solve(S) %*% vec.y )

        tSQSy <- c( t(Time) %*% solve(S) %*% Q %*% solve(S) %*% vec.y )

        yEtat <- c( vec.y - Eta.M3 * Time )

        M.A.temp <- matrix( 0, 3, 3 )

        M.A.temp[1,1] <- - tSt

        M.A.temp[1,2] <- M.A.temp[2,1] <- - tSt * tSy + Eta.M3 * tSt^2

        M.A.temp[1,3] <- M.A.temp[3,1] <- - tSQSy + Eta.M3 * tSQSt

        M.A.temp[2,2] <- tSt^2 / 2 - tSt * t(yEtat) %*% SttS %*% yEtat

        M.A.temp[2,3] <- M.A.temp[3,2] <- tSQSt / 2 - t(yEtat) %*% ( SQSttS + SttSQS ) %*% yEtat / 2

        M.A.temp[3,3] <- sum( diag( SQSQ ) ) / 2 - t(yEtat) %*% SQSQS %*% yEtat
 
        M.A.temp

    }


    A.cal.2 <- apply( Data[,-1], 2, A.cal )

    OI.A.M3 <- (-1) * matrix( apply( A.cal.2, 1, sum ), 3, 3 )

    return(OI.A.M3)

}


OI.A.M4 <-
function( Data, Eta.M4, sE.M4, Time, m, ...){

    tOt <- c( t(Time) %*% diag(m) %*% Time / sE.M4^2 )

    tOOt <- c( t(Time) %*% diag(m) %*% Time / sE.M4^4 )

    OOO <- diag(m) / sE.M4^6
    
    A.cal <- function( vec.y, ... ){

        tOOy <- c( t(Time) %*% diag(m) %*% vec.y / sE.M4^4 )

        yEtat <- c( vec.y - Eta.M4 * Time )

        M.A.temp <- matrix( 0, 2, 2 )

        M.A.temp[1,1] <- - tOt

        M.A.temp[1,2] <- M.A.temp[2,1] <- - tOOy + Eta.M4 * tOOt

        M.A.temp[2,2] <- m / ( 2 * sE.M4^4 ) - t(yEtat) %*% OOO %*% yEtat
 
        M.A.temp

    }


    A.cal.2 <- apply( Data[,-1], 2, A.cal )

    OI.A.M4 <- (-1) * matrix( apply( A.cal.2, 1, sum ), 2, 2 )

    return(OI.A.M4)

}


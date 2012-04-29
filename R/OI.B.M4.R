OI.B.M4 <-
function( Data, Eta.M4, sE.M4, Time, m, ...){

    tOt <- c( t(Time) %*% diag(m) %*% Time / sE.M4^2 )

    uu <- function( vec.y, ... ){

        u <- matrix(0,2,1)

        tOy <- c( t(Time) %*% diag(m) %*% vec.y / sE.M4^2 )

        yEtat <- ( vec.y - Eta.M4 * Time )

        u[1,1] <- tOy - Eta.M4 * tOt

        u[2,1] <- - m / ( 2 * sE.M4^2 ) + c( t(yEtat) %*% diag(m) %*% yEtat /

			( 2 * sE.M4^4 ) )

        u %*% t(u)
        
    }

    uu.cal <- apply( Data[,-1], 2, uu )

    OI.B.M4 <- matrix( apply( uu.cal, 1, sum ), 2, 2, byrow=TRUE )

    return(OI.B.M4)

}


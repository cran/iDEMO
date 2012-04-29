idemo.box <- function(data.name,data.in){

        data.in <- rbind(0,data.in)

        col_1 <- as.numeric(data.matrix(t(data.in))[-1,])

        col_2 <- rep( data.in[,1], each = (length(data.in[1,])-1) )

        newdata <- data.frame( col_1, col_2 )

        x11()

        boxplot(col_1 ~ col_2, data = newdata, xlab = 't', ylab = 'Y(t)', at = data.in[,1], boxwex = max(data.in[,1])/(2.5*nrow(data.in)), xlim = c(0, max(data.in[,1])),axes=FALSE)

        box()

        axis(1, data.in[,1], round(data.in[,1], digits=2))

        axis(2)

        points(data.in[,1], apply(data.in[,-1], 1, mean, na.rm = T), pch=2, cex=0.8, col=2 )

        lines(data.in[,1], apply(data.in[,-1], 1, mean, na.rm = T), lty=2,col=2 )

        title(paste('Box plot for ',data.name,' data',sep=""))

        legend("topleft","mean", pch=2, lty=2, merge=TRUE, col=2, cex=0.8, inset = .02, bg='white')

}
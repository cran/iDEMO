plot.degrad.path <- function(dataname, Data, ...){

        Data.new <- rbind(0, Data)

        x11()
        plot(Data.new[, 1],Data.new[, 2], xlab = "t", ylab = "Y(t)", type = "l", lty = 1, col = 1,
                ylim = c(min(Data.new[,-1], na.rm = T), max(Data.new[,-1], na.rm = T)),
                main = paste("Degradation path for ", dataname, " data", sep=""))

        addline.fun <- function(obs){
                ind <- which(!is.na(obs))
                lines(Data.new[ind,1],obs[ind])
        }
        apply(Data.new[,-1],2,addline.fun)

}


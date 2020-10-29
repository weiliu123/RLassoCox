rw <-
function(W,p0,gamma){
    p0 <- t(p0)
    p0 <- p0/sum(p0)
    PT <- p0 
    k <- 0
    delta <- 1  
    Ng <- dim(W)[2]
    for (i in seq(from = 1, to = Ng)){
        sumr<-sum(W[i,])
        if(sumr==0){
            W[i,] <-numeric(length=length(W[i,]));
        }
        if(sumr>0){
            W[i,] <- W[i,]/sumr
        }
    }
    W <- t(W)
    while(delta > 1e-10){
        PT1 <- (1-gamma)*W
        PT2 <- PT1 %*% t(PT)
        PT3 <- (gamma*p0)
        PT4 <- t(PT2) + PT3
        delta <- sum(abs(PT4 - PT))
        PT <- PT4
        k <- k + 1
    }  
    PT<-t(PT)
    rownames(PT)<-NULL
    return(PT)
}

library('smoothmest')
simu <- function(n,a, b,  sdis){
    x <- rbeta(n, a, b) 
    if(sdis == 0){
        eps <- rnorm(n, 0, sig)
    }else if(sdis == 1){
        eps <- rdoublex(n, 0, sig)
    }else{
        eps <- runif(n, -upper, upper)
    }
    y <- x + eps
    return(y)
}
rdoublex <- function (n, mu = 0, lambda = 1) 
{
    x <- rexp(n)
    s <- sample(c(-1, 1), size = n, replace = TRUE)
    y <- s * x
    lambda * y + mu
}
ddoublex <- function (x, mu = 0, lambda = 1){
	 exp(-abs(x - mu)/lambda)/(2 * lambda)
}
mden <- vector('list')
nsim <- 200
for(n in seq(500, 2000, 100)){
for(sdis in 1){
    mden[[sdis + 1]] <- matrix(NA, nsim, 100)
if(sdis == 0){
sig <- 0.5#sqrt(1/3)#norm 1
}else if(sdis == 1){
sig <- 0.5
}else{
    upper <- 0.25
    sig <-  sqrt((upper *2)^2/12)
}

    a <- 4
    b <- 4
    x <- seq(0, 1, length.out = 100)
    
for(itr in 1 : nsim){
    y <- simu(n, a, b, sdis)
   
    
    bdw <-  PI_deconvUknownth4(y,'norm',sig^2,sig)#CVdeconv(y, 'norm', sig)
    if(sdis == 1){
         bdw <-  PI_deconvUknownth4(y,'Lap',2 * sig^2,sqrt(2 * sig^2))
    }
    if(sdis == 0){
        density <- fdecUknown(x, y, bdw, 'norm', sig)
    }else if(sdis == 1){
         density <- fdecUknown(x, y, bdw, 'Lap', sqrt(2 * sig^2))
     }else{
         density <- fdecUknown(x, y, bdw, 'unif', sig)
     }
    mden[[sdis + 1]][itr, ] <- density
  #  mtesty[itr, ] <- testxy$y
   # mtestx[itr, ] <- testxy$x
    
   # covden[itr] <- max(abs(denw - testxy$y))
   # myden[itr] <- max(abs(estden - testxy$y))
    print(itr)
}
}
save(mden, file = paste('mden', n, sep= '_'))
}
f <-   dbeta(x, a, b)
n <- 2000
load(paste('mden', n, sep= '_'))
for(sdis in 1){
    pdf(paste('den', n, sdis, '.pdf', sep = ''))
    fest <- apply(mden[[sdis + 1]][, ],2, median, na.rm = T)
    qfest <- apply(mden[[sdis + 1]][, ],2,quantile, c(0.05, 0.95),  na.rm = T)
    plot(f~x, type = 'l', ylim = c(0, 10))
    lines((fest ~x), lty = 2)

    lines(x, qfest[1, ], lty = 2)
    lines(x, qfest[2, ], lty = 2)
    dev.off()
}
mf <- matrix(f, nsim, length(x), byrow = T)
ln <- seq(500, 2000, 100)
l <- length(ln)
nh <- sdiff1 <- matrix(NA, 3, l)
for(i in 1:l){
    n <- ln[i]
    load(paste('mden', n, sep = '_'))
for(sdis in 1){
    
  
    sdiff1[sdis+1, i] <- mean(apply(abs(mden[[sdis+ 1]] - mf), 1,max, na.rm = T), na.rm = T)
    nh[sdis, i] <- bdw
    }
}
for(sdis in 1){
    if(sdis ==1){
        plot(sdiff[sdis+ 1, -c(1:3)]~ln, type= 'l', ylab = 'maximum absolute error', ylim = range(c(sdiff[sdis + 1, ], sdiff1[sdis + 1, ])), xlab = 'sample size')
        lines(sdiff1[sdis + 1, ] ~ln, lty = 2)
        
    }
}


simu <- function(n,a, b,  sdis){
    x <-  rbeta(n, a, b) 
    if(sdis == 0){
        eps <- rnorm(n, 0, sigw)
    }else if(sdis == 1){
        eps <- rdoublex(n, 0, sigw)
    }else{
        eps <- runif(n, ul, up)
    }
    w <- x + eps
    mx <- sin(x * 2 * pi)#1 - 28.76 * x + 251.97 * x^2 - 1028.50 * x^3 + 2069.72 * x^4 - 1975.95 * x^5 + 712.78 * x^6#  sin(x * 2 * pi)# x  * (x - 0.5) * (x - 1) #sin(x * 2  * pi) #dbeta(x, 2, 2)#2* sin((x - 0.5)  * 12) * exp(-((x - 0.5) * 12)^2/10)##
    y <- mx  + rnorm(n, 0, sig) 
    return(cbind(y, w))
}


a <- 4
b <- 4
x <- seq(0, 1, length.out = 100)
mmean <- vector('list')
for(n in seq(500, 2000, 100)){
for(sdis in 0:2){
    mmean[[sdis + 1]] <- matrix(NA, nsim, 100)
    upper = 1
    up <- 1/4
    ul <- -1/4
    sig = 0.5
    sigw <- 0.5
for(itr in 1:nsim){
    ywdata <- simu(n, a, b, sdis)
     bdw <- 0.005#CVdeconv(y, 'norm', sig)sd(ywdata[, 2]) * n^(-1/5) * 1.06
    if(sdis == 0){
        bdw <-  PI_deconvUknownth4(ywdata[, 2],'norm',sigw^2,sigw)
         m <- NWDec(x, ywdata[, 2], ywdata[, 1], 'norm',  sigw, bdw, 1e-6 )
    }else if(sdis == 1){
        bdw <-  PI_deconvUknownth4(ywdata[, 2],'Lap',2 * sigw^2,sqrt(2 * sigw^2))
          m <- NWDec(x, ywdata[, 2], ywdata[, 1], 'Lap',  sqrt(2 * sigw^2), bdw, 1e-6 )
    }else{
        bdw <-  PI_deconvUknownth4(ywdata[, 2],'unif', (up - ul)^2/12, sqrt((up - ul)^2/12))
          m <- NWDec(x, ywdata[, 2], ywdata[, 1], 'unif',  sqrt((up - ul)^2/12), bdw, 1e-6 )
     }
    mmean[[sdis + 1]][itr, ] <- m
     print(itr)
}
   
}
save(mmean, file = paste('mmean', n, sep= '_'))
}




f <- sin(2 * pi * x)#1 - 28.76 * x + 251.97 * x^2 - 1028.50 * x^3 + 2069.72 * x^4 - 1975.95 * x^5 + 712.78 * x^6
for(sdis in 0:2){
    pdf(paste('deconvmean', n, sdis, '.pdf', sep = ''))
    fest <- apply(mmean[[sdis + 1]],2,mean, na.rm = T)
    qfest <- apply(mmean[[sdis + 1]],2,quantile, c(0.025, 0.975),  na.rm = T)
    plot(f~x, type = 'l', ylim = c(-2, 3), ylab= 'm')
    lines((fest ~x), lty = 2)

    lines(x, qfest[1, ], lty = 2)
    lines(x, qfest[2, ], lty = 2)
    dev.off()
}

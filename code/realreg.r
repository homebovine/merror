library(splines)
library(smoothmest)
library(optimx)
library(statmod)
library(numDeriv)
library(fda)
library(minpack.lm)
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
    mx <- 1 - 28.76 * x + 251.97 * x^2 - 1028.50 * x^3 + 2069.72 * x^4 - 1975.95 * x^5 + 712.78 * x^6#  sin(x * 2 * pi)# x  * (x - 0.5) * (x - 1) #sin(x * 2  * pi) #dbeta(x, 2, 2)#2* sin((x - 0.5)  * 12) * exp(-((x - 0.5) * 12)^2/10)##
    y <- mx  + rnorm(n, 0, sig) 
    return(cbind(y, w))
}
ddnorm <- function(x){
    1/sqrt(2 * pi * sig^2) * exp(-x^2/(2 *sig^2)) * (-x/sig^2)
}
sdlap <- function(x, m, s){
      1/(2 * s) * exp( -(x- m)/s) * pnorm((x - m)/h, 0, 1)  + 1/(2 * s) * exp( -(m- x)/s) * (1 - pnorm((x - m)/h, 0, 1) )

}
ahmatrix <- function(theta, sdis, bs, ywxc, ywxxc, tol){
    bst1<- (basis1 %*% theta)
    bst2 <- (basis2 %*% theta)
    bst <- (basis %*% theta)
    bstah <- basisah%*% theta
    suml <- function(c, x, y, w, x2,  sdis){
#    	 oy <- log((1 + y)/(1 - y))
#    	 ow <- log((1 + w)/(1 - w))	
        if(sdis ==0 ){
            deps <- dnorm(w * sigw + x2 - x, 0, sigw) * dnorm(y * sig * sqrt(2) + bstah - bst, 0, sig) * c #* -2/(w ^2 - 1) * -2/(y^2 - 1)
        }else if (sdis == 1 ){
            deps <- ddoublex(w * sigw + x2 -x, 0, sigw) * dnorm(y* sig * sqrt(2) + bstah - bst, 0, sig) * c 
        }else {
             deps <- dunif(w - x, ul, up) * dnorm(y * sig * sqrt(2) + bstah - bst, 0, sig) * c        }
    }
    sumL <- suml(ywxc[, 1], ywxc[, 2], ywxc[, 3], ywxc[, 4], ywxc[, 5],  sdis)
     sumL <- apply(matrix(sumL, ncol = llx, byrow = T), 1, sum)
    sumL <- rep(sumL,  llx)
    
    sumL[sumL == 0] = sumL[sumL == 0] + 1e-16
    integrand1 <- function(x1, x2,c,  y, w, sdis){
    	       #oy <- log((1 + y)/(1 - y))
    	       #ow <- log((1 + w)/(1 - w))	
    if(sdis == 0 ){
        deps <-  dnorm(w* sigw + x1 - x2, 0, sigw) *c *  dnorm(y * sig * sqrt(2)+ bst1 - bst2,  0, sig)/sumL# * dnorm(w, 0, 1) *   dnorm(y,  0, 1) 
        ddeps <-  dnorm(w* sigw + x1 - x2, 0, sigw) *c *  ddnorm(y * sig * sqrt(2)+ bst1 - bst2)/sumL# * dnorm(w, 0, 1) *   dnorm(y,  0, 1)
    }
    else if(sdis == 1){
         deps <-  ddoublex(w* sigw + x1 - x2, 0, sigw) *c *  dnorm(y * sig * sqrt(2)+ bst1 - bst2,  0, sig)/sumL# * ddoublex(w, 0, 1) *   dnorm(y,  0, 1)
        ddeps <- ddoublex(w* sigw + x1 - x2, 0, sigw) *c *  ddnorm(y * sig * sqrt(2)+ bst1 - bst2)/sumL #* ddoublex(w - x1, 0, sigw) *   dnorm(y - bst1,  0, sig)
        
    }else{
        deps <-  dunif(w - x2, ul, up) *c *  dnorm(y * sig * sqrt(2) + bst1 -  bst2,  0, sig)/sumL * dunif(w - x1, ul, up) #*   dnorm(y - bst1,  0, sig)
        ddeps <- dunif(w - x2, ul, up) *c *  ddnorm(y * sig * sqrt(2) + bst1 - bst2)/sumL * dunif(w - x1, ul, up) #*   dnorm(y - bst1,  0, sig) 
       
    }
    return(cbind(deps, ddeps))
}
    temp <- integrand1(ywxxc[, 3], ywxxc[, 4], ywxxc[, 5], ywxxc[, 1], ywxxc[, 2], sdis)
    mA <- matrix(temp[, 1], ncol = lly * llw, byrow = T) %*% Delta
    mH <- matrix(temp[, 2], ncol = lly * llw, byrow = T) %*% Delta
    mA <- matrix(mA,  llx, llx)
    mH <- -matrix(mH,  llx, llx)
    mb <- matrix(NA, lb, llx)
    
        mb <- t(mH %*% basisx)
    
    ma = mb %*% mA %*% ginv(t(mA) %*% mA, tol)#### ##mb %*% ginv(t(mA), tol)#
   return(ma)
    
}

ahmatrix1 <- function(theta, sdis, bs, ywxc, ywxxc, tol){
    bst1<- (basis1 %*% theta)
    bst2 <- (basis2 %*% theta)
    bst <- (basis %*% theta)
    bstah <- basisah%*% theta
    suml <- function(c, x, y, w, x2,  sdis){
#    	 oy <- log((1 + y)/(1 - y))
#    	 ow <- log((1 + w)/(1 - w))	
        if(sdis ==0 ){
            deps <- dnorm(w * sigw + x2 - x, 0, sigw) * dnorm(y * sig * sqrt(2) + bstah - bst, 0, sig) * c #* -2/(w ^2 - 1) * -2/(y^2 - 1)
        }else if (sdis == 1 ){
            deps <- ddoublex(w * sigw + x2 -x, 0, sigw) * dnorm(y* sig * sqrt(2)+ bstah - bst, 0, sig) * c 
        }else {
             deps <- dunif(w - x, ul, up) * dnorm(y * sig * sqrt(2) + bstah  - bst, 0, sig) * c        }
    }
    sumL <- suml(ywxc[, 1], ywxc[, 2], ywxc[, 3], ywxc[, 4], ywxc[, 5],  sdis)
     sumL <- apply(matrix(sumL, ncol = llx, byrow = T), 1, sum)
    sumL <- rep(sumL,  llx)
    
    sumL[sumL == 0] = sumL[sumL == 0] + 1e-16
    integrand1 <- function(x1, x2,c,  y, w, sdis){
    	       #oy <- log((1 + y)/(1 - y))
    	       #ow <- log((1 + w)/(1 - w))	
    if(sdis == 0 ){
        deps <-  dnorm(w* sigw + x1 - x2, 0, sigw) *c *  dnorm(y * sig* sqrt(2) + bst1 - bst2,  0, sig)/sumL #* dnorm(w, 0, 1) *   dnorm(y,  0, 1) 
        ddeps <-  dnorm(w* sigw + x1 - x2, 0, sigw) *c *  ddnorm(y * sig * sqrt(2)+ bst1 - bst2)/sumL #* dnorm(w, 0, 1) *   dnorm(y,  0, 1)
    }
    else if(sdis == 1){
         deps <-  ddoublex(w* sigw + x1 - x2, 0, sigw) *c *  dnorm(y * sig * sqrt(2)+ bst1 - bst2,  0, sig)/sumL# * ddoublex(w, 0, 1) *   dnorm(y,  0, 1)
        ddeps <- ddoublex(w* sigw + x1 - x2, 0, sigw) *c *  ddnorm(y * sig * sqrt(2)+ bst1 - bst2)/sumL #* ddoublex(w - x1, 0, sigw) *   dnorm(y - bst1,  0, sig)
        
    }else{
        deps <-  dunif(w - x2, ul, up) *c *  dnorm(y * sig * sqrt(2) + bst1 - bst2,  0, sig)/sumL * dunif(w - x1, ul, up) #*   dnorm(y - bst1,  0, sig)
        ddeps <- dunif(w - x2, ul, up) *c *  ddnorm(y * sig * sqrt(2) + bst1 - bst2)/sumL * dunif(w - x1, ul, up) #*   dnorm(y - bst1,  0, sig) 
       
    }
    return(cbind(deps, ddeps))
}
    temp <- integrand1(ywxxc[, 3], ywxxc[, 4], ywxxc[, 5], ywxxc[, 1], ywxxc[, 2], sdis)
    mA <- matrix(temp[, 1], ncol = lly * llw, byrow = T) %*% Delta
    mH <- matrix(temp[, 2], ncol = lly * llw, byrow = T) %*% Delta
    mA <- matrix(mA,  llx, llx)
    mH <- -matrix(mH,  llx, llx)
    mb <- matrix(NA, lb, llx)
    
        mb <- t(mH %*% basisx)
    
    ma = mb %*% ginv(t(mA), tol)#
   return(ma)
    
}





meanest <- function(theta,  my, mw, mc, sdis, tol, ahmatrix){
 #   print(theta)
    ma <- ahmatrix(theta, sdis, bs, ywxc, ywxxc, tol)
    bst <- matrix(basisx %*% theta, n, llx, byrow = T)
    if(sdis == 0 ){
        deps <-  dnorm(mw - mlx, 0, sigw) *  dnorm(my - bst,  0, sig) * mc
        ddeps <-  dnorm(mw -mlx, 0, sigw) *  ddnorm(my - bst) * mc
       
    }else if(sdis == 1){
          deps <-  ddoublex(mw - mlx, 0, sigw) *  dnorm(my - bst,  0, sig) * mc
          ddeps <-  ddoublex(mw - mlx, 0, sigw) *  ddnorm(my - bst) * mc
         
        
    }else{
	
        deps <-   dunif(mw - mlx, ul, up) *  dnorm(my - bst,  0, sig) * mc 
        ddeps <-   dunif(mw - mlx, ul, up) *  ddnorm(my - bst) * mc 
       
    }
    numa <- deps %*% t(ma)
    dom <- matrix(apply(deps, 1, sum), ncol = lb, nrow = n)
    dom[dom == 0] = dom[dom == 0] + 1e-16
    nums <- ddeps %*% basisx
    S <- -nums/dom - numa/dom	
    #wt <<-ginv( t(S) %*% S/n, 1e-16)
    S <- apply(S, 2, mean)
    as.numeric(t(S) %*% wt %*% S)
}

smeanest <- function(theta,  my, mw, mc, sdis, tol, ahmatrix){
 #   print(theta)
# theta <- pnorm(theta) * 4 - 2
    ma <- ahmatrix(theta, sdis, bs, ywxc, ywxxc, tol)
    bst <- matrix(basisx %*% theta, n, llx, byrow = T)
    if(sdis == 0 ){
        deps <-  dnorm(mw - mlx, 0, sigw) *  dnorm(my - bst,  0, sig) * mc
        ddeps <-  dnorm(mw -mlx, 0, sigw) *  ddnorm(my - bst) * mc
       
    }else if(sdis == 1){
          deps <-  ddoublex(mw - mlx, 0, sigw) *  dnorm(my - bst,  0, sig) * mc
          ddeps <-  ddoublex(mw - mlx, 0, sigw) *  ddnorm(my - bst) * mc
         
        
    }else{
	
        deps <-   dunif(mw - mlx, ul, up) *  dnorm(my - bst,  0, sig) * mc
        ddeps <-   dunif(mw - mlx, ul, up) *  ddnorm(my - bst) * mc 
       
    }
    numa <- deps %*% t(ma)
    dom <- matrix(apply(deps, 1, sum), ncol = lb, nrow = n)
    dom[dom == 0] = dom[dom == 0] + 1e-16
    nums <- ddeps %*% basisx
    S <- apply(-nums/dom - numa/dom, 2, mean)
}


ssmeanest <- function(theta,  my, mw, mc, sdis, tol, ahmatrix){
 #   print(theta)
# theta <- pnorm(theta) * 4 - 2
    ma <- ahmatrix(theta, sdis, bs, ywxc, ywxxc, tol)
    bst <- matrix(basisx %*% theta, n, llx, byrow = T)
    if(sdis == 0 ){
        deps <-  dnorm(mw - mlx, 0, sigw) *  dnorm(my - bst,  0, sig) * mc
        ddeps <-  dnorm(mw -mlx, 0, sigw) *  ddnorm(my - bst) * mc
       
    }else if(sdis == 1){
          deps <-  ddoublex(mw - mlx, 0, sigw) *  dnorm(my - bst,  0, sig) * mc
          ddeps <-  ddoublex(mw - mlx, 0, sigw) *  ddnorm(my - bst) * mc
         
        
    }else{
	
        deps <-   dunif(mw - mlx, ul, up) *  dnorm(my - bst,  0, sig) * mc
        ddeps <-   dunif(mw - mlx, ul, up) *  ddnorm(my - bst) * mc 
       
    }
    numa <- deps %*% t(ma)
    dom <- matrix(apply(deps, 1, sum), ncol = lb, nrow = n)
    dom[dom == 0] = dom[dom == 0] + 1e-16
    nums <- ddeps %*% basisx
    S <- -nums/dom - numa/dom
    mS <- t(S) %*% S /n
}


myjacobian <- function(func, x, method="Richardson", side=NULL,
      method.args=list(), ...){
  f <- func(x, ...)
  n <- length(x)	 #number of variables.

  if (is.null(side)) side <- rep(NA, n)
  else {
       if(n != length(side)) 
          stop("Non-NULL argument 'side' should have the same length as x")
       if(any(1 != abs(side[!is.na(side)]))) 
          stop("Non-NULL argument 'side' should have values NA, +1, or -1.")
       }

  if(method=="simple"){
    #  very simple numerical approximation
    args <- list(eps=1e-4) # default
    args[names(method.args)] <- method.args

    side[is.na(side)] <- 1
    eps <- (args$eps) * side

    df <-matrix(NA, length(f), n)
    for (i in 1:n) {
      dx <- x
      dx[i] <- dx[i] + eps[i] 
      df[,i] <- (func(dx, ...) - f)/eps[i]
     }
    return(df)
    } 
  else if(method=="complex"){ # Complex step gradient
    if (any(!is.na(side))) stop("method 'complex' does not support non-NULL argument 'side'.")
    # Complex step Jacobian
    eps <- .Machine$double.eps
    h0  <-  rep(0, n)
    h0[1] <- eps * 1i
    v <- try(func(x+h0, ...))
    if(inherits(v, "try-error")) 
      stop("function does not accept complex argument as required by method 'complex'.")
    if(!is.complex(v)) 
      stop("function does not return a complex value as required by method 'complex'.")
  
    h0[1]  <- 0
    jac <- matrix(NA, length(v), n)
    jac[, 1] <- Im(v)/eps
    if (n == 1) return(jac)
    for (i in 2:n) {
      h0[i] <- eps * 1i
      jac[, i] <- Im(func(x+h0, ...))/eps 
      h0[i]  <- 0
      }
    return(jac)
    } 
  else if(method=="Richardson"){
    args <- list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), 
                r=4, v=2, show.details=FALSE) # default
    args[names(method.args)] <- method.args
    d <- args$d
    r <- args$r
    v <- args$v		  
    a <- array(NA, c(length(f),r, n) )
  
    h <- abs(d*x) + args$eps * (abs(x) < args$zero.tol)
    pna <- (side == 1)  & !is.na(side) # double these on plus side
    mna <- (side == -1) & !is.na(side) # double these on minus side

    for(k in 1:r)  { # successively reduce h	       
       ph <- mh <- h
       ph[pna] <- 2 * ph[pna] 
       ph[mna] <- 0           
       mh[mna] <- 2 * mh[mna] 
       mh[pna] <- 0           

       for(i in 1:n)  {
        a[,k,i] <- (func(x + ph*(i==seq(n)), ...) -  
     		     func(x - mh*(i==seq(n)), ...))/(2*h[i])
    		      #if((k != 1)) a[,(abs(a[,(k-1),i]) < 1e-20)] <- 0 #some func are unstable near zero
    		       }
       h <- h/v     # Reduced h by 1/v.
       }     

   for(m in 1:(r - 1)) {	  
       a <- (a[,2:(r+1-m),,drop=FALSE]*(4^m)-a[,1:(r-m),,drop=FALSE])/(4^m-1)
     }
  # drop second dim of a, which is now 1 (but not other dim's even if they are 1
  return(array(a, dim(a)[c(1,3)]))  
  } else stop("indicated method ", method, "not supported.")
}



sdis <- 0

nsim <- 100
lly <- 15
llw <- 15

load('framdata')
n <- ceiling(nrow(framdata))

sig <- sqrt(mean(lm(framdata[, 1] ~ framdata[, 2])$residuals^2))
#sigw <- 0.0646
       
       
gy <- gauss.quad(lly,'hermite', alpha=0,beta=0)
gw <- gauss.quad(llw,'hermite', alpha=0,beta=0)
ly <- gy$nodes
lw <- gw$nodes
Delta <-   as.vector((gy$weights ) %o%  (gw$weights ))
       
x <-  rbeta(n, 1, 1) 
o <- order(x)
x <- x[o]
lux <- pbeta(x, 1, 1)
lx <-  qbeta(seq(0, 1, length.out = 15), 1, 1) ##seq(0.01, 1, length.out = 15)#
cdom <- sum(dbeta(lx, 1, 1))
olc <- dbeta(lx, 1, 1)/cdom
ywxx <- expand.grid(ly, lw, lx, lx)
lc <- dbeta(ywxx[, 4],1, 1)/cdom
ywxxc <- cbind(ywxx, lc)
ywx <- expand.grid(lx, ly, lw, lx)
lc <- dbeta(ywx[, 1] , 1, 1)/cdom
ywxc <- cbind(lc, ywx)
       
lb = 9# 4 * (n <= 500) + 6 * (n>500 & n <= 800) + 7 * (n >800 & n <= 1000)+ 8 * (n >1000 & n <= 1500)  + 9 * (n > 1500)

b2 <- create.bspline.basis(range(x, lx), nbasis = lb, norder = 4)
basis2 <- eval.basis(ywxxc[, 4], b2)
basis1 <-eval.basis(ywxxc[, 3], b2)
basis <- eval.basis(ywxc[, 2], b2) 
basisx <-  eval.basis(lx, b2)
basisah <- eval.basis(ywxc[, ncol(ywxc)], b2)
llx <- length(lx)
lly <- length(ly)
llw <- length(lw)
f <- framdata[framdata[, 2] <=1, 1]
tempba <-eval.basis(framdata[framdata[, 2] <=1, 2], b2) 
tempf <- function(theta){
    sum((tempba %*% theta - f)^2)
}
theta0 <- optim(rep(0.5, lb), tempf, gr = NULL,  method = 'BFGS',  hessian = FALSE) $par
resmean <- matrix(NA, nsim, lb)

for(itr in 1:nsim){
    print(itr)
    set.seed(itr + 2017)
    train<- sample(1:nrow(framdata), n, replace = T)
    
   
    ywdata <- framdata[train, ]
    #testy <- framdata[-train, ]

    my <- matrix(ywdata[, 1], n, llx)
    mw <- matrix(ywdata[, 2], n, llx)
    mlx <- matrix(lx, nrow = n, ncol = llx, byrow = T)
    mc <- matrix(olc, nrow = n, ncol = llx, byrow = T)
    wt <- diag(1, lb, lb)#ginv(ssmeanest(theta0, my, mw, mc, sdis, 1e-8, ahmatrix1), 1e-8)### 

    tryres <- optimx(theta0, meanest, gr = NULL, hess = NULL, lower = -Inf, upper = Inf, method = 'newuoa',  itnmax=500, hessian=FALSE, control=list(), my = my, mw = mw, mc = mc, sdis = sdis, tol = 1e-8, ahmatrix = ahmatrix1)
    tryres$par <- coef(tryres)
    tryres$convergence <- tryres$convcode
    print(tryres$convergence)
            if(class(tryres) != 'try-error'){
                resmean[itr, ] <- tryres$par
                                        
            }
            print(rbind(resmean[itr, ],theta0))
        }
        save(resmean,  tempba, x, lb, train, file = paste('realreg2'))


resmean[itr, ] %*% t(tempba)









nl <- seq(200, 2000, 100)# c(seq(200, 1100, 100), seq(1200, 3000, 200))
l <- length(nl)
nhn <- sdiff <- matrix(NA, 3, l)

for(sdis in 0:2){
for(i in 1:l){
n <- nl[i]
load(paste('resmean', n, sdis, sep = '_'))
 nsim <- 200
x <- seq(0, 1, length.out = 100)
b2 <- create.bspline.basis(c(0, 1), nbasis = lb, norder = 4)
tempba <- eval.basis(x, b2)
#lb = ceiling(1.3 * n^{1/5})#4 * (n <= 500) + 5 * (n>500 & n <= 800) + 6 * (n >800)#ceiling(1.3 * n^(1/5))
nhb = 1/lb

f <- 1 - 28.76 * x + 251.97 * x^2 - 1028.50 * x^3 + 2069.72 * x^4 - 1975.95 * x^5 + 712.78 * x^6#sin(2 * pi * x)

fest <- apply(abs(resmean[, ] %*% t(tempba) - matrix(f, nsim, length(x), byrow = T)), 1, max)# apply(resmean[, ] %*% t(tempba), 2, mean, na.rm = T)
sdiff[sdis+ 1, i]<- median(fest) 
nhn[sdis+ 1, i] <- sqrt(nhb * n)
}
} 
size <- (seq(200, 2000, 100))
sdiff1 <- sdiff * nhn
pdf('temp1.pdf')
plot(sdiff1[1, ]~ size, type = 'l', ylim = c(min(sdiff1),max(sdiff1)))
lines(sdiff1[2, ] ~size, lty = 4)
lines(sdiff1[3, ]~size, lty = 6)
legend('topleft', c('normal error', 'laplace error', 'uniform error'), lty = c(1, 4, 6))
dev.off()


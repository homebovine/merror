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
    mx <-  sin(x * 2 * pi)# x  * (x - 0.5) * (x - 1) #sin(x * 2  * pi) #dbeta(x, 2, 2)#2* sin((x - 0.5)  * 12) * exp(-((x - 0.5) * 12)^2/10)##1 - 28.76 * x + 251.97 * x^2 - 1028.50 * x^3 + 2069.72 * x^4 
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
    suml <- function(c, x, y, w,  sdis){
#    	 oy <- log((1 + y)/(1 - y))
#    	 ow <- log((1 + w)/(1 - w))	
        if(sdis ==0 ){
            deps <- dnorm(w - x, 0, sigw) * dnorm(y - bst, 0, sig) * c #* -2/(w ^2 - 1) * -2/(y^2 - 1)
        }else if (sdis == 1 ){
            deps <- ddoublex(w - x, 0, sigw) * dnorm(y - bst, 0, sig) * c 
        }else {
             deps <- dunif(w - x, ul, up) * dnorm(y - bst, 0, sig) * c        }
    }
    sumL <- suml(ywxc[, 1], ywxc[, 2], ywxc[, 3], ywxc[, 4],  sdis)
     sumL <- apply(matrix(sumL, ncol = llx, byrow = T), 1, sum)
    sumL <- rep(sumL,  llx^2)
    sumL[sumL == 0] = sumL[sumL == 0] + 1e-16
    integrand1 <- function(x1, x2,c,  y, w, sdis){
    	       #oy <- log((1 + y)/(1 - y))
    	       #ow <- log((1 + w)/(1 - w))	
    if(sdis == 0 ){
        deps <-  dnorm(w - x2, 0, sigw) *c *  dnorm(y - bst2,  0, sig)/sumL * dnorm(w - x1, 0, sigw) *   dnorm(y - bst1,  0, sig) 
        ddeps <-  dnorm(w - x2, 0, sigw) *c *  ddnorm(y - bst2)/sumL * dnorm(w - x1, 0, sigw) *   dnorm(y - bst1,  0, sig)
    }
    else if(sdis == 1){
         deps <-  ddoublex(w - x2, 0, sigw) *c *  dnorm(y - bst2,  0, sig)/sumL * ddoublex(w - x1, 0, sigw) *   dnorm(y - bst1,  0, sig)
        ddeps <- ddoublex(w - x2, 0, sigw) *c *  ddnorm(y - bst2)/sumL * ddoublex(w - x1, 0, sigw) *   dnorm(y - bst1,  0, sig)
        
    }else{
        deps <-  dunif(w - x2, ul, up) *c *  dnorm(y - bst2,  0, sig)/sumL * dunif(w - x1, ul, up) *   dnorm(y - bst1,  0, sig)
        ddeps <- dunif(w - x2, ul, up) *c *  ddnorm(y - bst2)/sumL * dunif(w - x1, ul, up) *   dnorm(y - bst1,  0, sig) 
       
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
    suml <- function(c, x, y, w,  sdis){
#    	 oy <- log((1 + y)/(1 - y))
#    	 ow <- log((1 + w)/(1 - w))	
        if(sdis ==0 ){
            deps <- dnorm(w - x, 0, sigw) * dnorm(y - bst, 0, sig) * c #* -2/(w ^2 - 1) * -2/(y^2 - 1)
        }else if (sdis == 1 ){
            deps <- ddoublex(w - x, 0, sigw) * dnorm(y - bst, 0, sig) * c 
        }else {
             deps <- dunif(w - x, ul, up) * dnorm(y - bst, 0, sig) * c        }
    }
    sumL <- suml(ywxc[, 1], ywxc[, 2], ywxc[, 3], ywxc[, 4],  sdis)
     sumL <- apply(matrix(sumL, ncol = llx, byrow = T), 1, sum)
    sumL <- rep(sumL,  llx^2)
    sumL[sumL == 0] = sumL[sumL == 0] + 1e-16
    integrand1 <- function(x1, x2,c,  y, w, sdis){
    	       #oy <- log((1 + y)/(1 - y))
    	       #ow <- log((1 + w)/(1 - w))	
    if(sdis == 0 ){
        deps <-  dnorm(w - x2, 0, sigw) *c *  dnorm(y - bst2,  0, sig)/sumL * dnorm(w - x1, 0, sigw) *   dnorm(y - bst1,  0, sig) 
        ddeps <-  dnorm(w - x2, 0, sigw) *c *  ddnorm(y - bst2)/sumL * dnorm(w - x1, 0, sigw) *   dnorm(y - bst1,  0, sig)
    }
    else if(sdis == 1){
         deps <-  ddoublex(w - x2, 0, sigw) *c *  dnorm(y - bst2,  0, sig)/sumL * ddoublex(w - x1, 0, sigw) *   dnorm(y - bst1,  0, sig)
        ddeps <- ddoublex(w - x2, 0, sigw) *c *  ddnorm(y - bst2)/sumL * ddoublex(w - x1, 0, sigw) *   dnorm(y - bst1,  0, sig)
        
    }else{
        deps <-  dunif(w - x2, ul, up) *c *  dnorm(y - bst2,  0, sig)/sumL * dunif(w - x1, ul, up) *   dnorm(y - bst1,  0, sig)
        ddeps <- dunif(w - x2, ul, up) *c *  ddnorm(y - bst2)/sumL * dunif(w - x1, ul, up) *   dnorm(y - bst1,  0, sig) 
       
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
mn<- seq(200, 1000, 100)
sdis <- 0
a <-  4
b <- 4

nsim <- 500
lly <- 10
llw <- 10
up = 1#sqrt(0.025 * 12)/2
ul = -1#-sqrt(0.025 * 12)/2
u = (up + ul)/2
v = abs(up - ul)/2
for(sdis in c(0, 1, 2)){
for( i in 1:length(mn)){
n <- mn[i]
h <- n^{-1/5}
sig <- 1
sigw <- 1.5
if(sdis == 1){
sig <- 0.5#sqrt((0.5)^2/2) #lap = 1 #norm 1
sigw <- 1#sqrt((0.5)^2/2)#norm sqrt((0.5)^2/2) lap 1.48#norm 1.5
}
if(sdis == 0 |sdis == 1 ){
gy <- gauss.quad(lly,alpha=0,beta=0)
gw <- gauss.quad(llw,alpha=0,beta=0)
ly <- log((1 + gy$nodes)/(1 - gy$nodes)) #sqrt( pi) *pnorm(gy$nodes, 0, 1/sqrt(2))
lw <- log((1 + gw$nodes)/(1 - gw$nodes))#sqrt( pi) *pnorm(gw$nodes, 0, 1/sqrt(2))
Delta <-   as.vector((gy$weights * 2/(1 - gy$nodes^2)) %o%  (gw$weights * 2/(1 - gw$nodes^2)))
}else if (sdis == 1){
myGrid <- createNIGrid(dim=2, type=c("oNC1", "oNC1"), level=c(6, 6))
lyw <- getNodes(myGrid)
ly <- log(lyw[, 1]/(1- lyw[, 1]))
lw <- log(lyw[, 2]/(1- lyw[, 2]))
ly <-  unique(ly)
lw <- unique(lw)
Delta <- getWeights(myGrid) * 1/(lyw[, 1] * (1 - lyw[, 1])) * 1/(lyw[, 2] * (1 - lyw[, 2]))
}else{
gy <- gauss.quad(lly,alpha=0,beta=0)
gw <- gauss.quad(llw)
ly <- log((1 + gy$nodes)/(1 - gy$nodes)) #sqrt( pi) *pnorm(gy$nodes, 0, 1/sqrt(2))
lw <- v * (gw$nodes) + u
Delta <-   as.vector((gy$weights * 2/(1 - gy$nodes^2)) %o%  (gw$weights * v))
}
x <-  rbeta(n, a, b) 
o <- order(x)
x <- x[o]
lux <- pbeta(x, a, b)
lx <-  qbeta(seq(0, 1, length.out = 10), a, b) ##seq(0.01, 1, length.out = 15)#
cdom <- sum(dbeta(lx, a, b))
olc <- dbeta(lx, a, b)/cdom
ywxx <- expand.grid(ly, lw, lx, lx)
lc <- dbeta(ywxx[, 4], a, b)/cdom
ywxxc <- cbind( ywxx, lc)
ywx <- expand.grid(lx, ly, lw)
lc <- dbeta(ywx[, 1] , a, b)/cdom
ywxc <- cbind(lc, ywx)
knots = c(0.25,  0.5,  0.75)#c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
lb = round(1.3 * n^{1/5})#4 * (n <= 500) + 5 * (n>500 & n <= 800) + 6 * (n >800)#ceiling(1.3 * n^(1/5))
b2 <- create.bspline.basis(range(x, lx), nbasis = lb, norder = 4)
basis2 <- eval.basis(ywxxc[, 4], b2)#bs(ywxxc[, 4], 3, knots, Boundary.knots = c(0, 1))
basis1 <-eval.basis(ywxxc[, 3], b2)# bs(ywxxc[, 3], 3, knots, Boundary.knots = c(0, 1))
basis <- eval.basis(ywxc[, 2], b2) # bs(ywxc[, 2], 3, knots, Boundary.knots = c(0, 1))
basisx <-  eval.basis(lx, b2) #bs(lx, 3, knots, Boundary.knots = c(0, 1))
llx <- length(lx)
lly <- length(ly)
llw <- length(lw)
f = sin(x * 2  * pi)
tempba <-eval.basis(x, b2) #bs(x, 3, knots, Boundary.knots = c(0, 1))
tempf <- function(theta){
    sum((tempba %*% theta - f)^2)
}
theta0 <- optim(rep(0.5, lb), tempf, gr = NULL,  method = 'BFGS',  hessian = FALSE) $par
resmean <- matrix(NA, nsim, lb)
resvar <- array(0, c(lb, lb, nsim))
for(itr in 1:nsim){
    print(itr)
set.seed(itr + 2017)
ywdata <- simu(n, a, b, sdis)
my <- matrix(ywdata[, 1], n, llx)
mw <- matrix(ywdata[, 2], n, llx)
mlx <- matrix(lx, nrow = n, ncol = llx, byrow = T)
mc <- matrix(olc, nrow = n, ncol = llx, byrow = T)
wt <- ginv(ssmeanest(theta0, my, mw, mc, sdis, 1e-16, ahmatrix), 1e-16)##diag(1, lb, lb)# 

tryres <- optimx(theta0, meanest, gr = NULL, hess = NULL, lower = -Inf, upper = Inf, method = 'newuoa',  itnmax=500, hessian=FALSE, control=list(), my = my, mw = mw, mc = mc, sdis = sdis, tol = 1e-16, ahmatrix = ahmatrix)
tryres$par <- coef(tryres)
tryres$convergence <- tryres$convcode
print(tryres$convergence)
    if(class(tryres) != 'try-error'){
        resmean[itr, ] <- tryres$par
# jh <- myjacobian(smeanest, resmean[itr, ], method="simple", side=NULL, method.args=list(eps =  c(1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4)), my, mw, mc, sdis, 1e-16, ahmatrix)
	# sh <- ssmeanest(resmean[itr, ], my, mw, mc, sdis, 1e-8, ahmatrix)
	# jwh <- ginv(t(jh) %*% wt %*% jh, 1e-14)
	 resvar[, , itr] <- 0# jwh %*% t(jh) %*% wt %*% sh %*% t(wt) %*% jh %*% t(jwh)
    }
print(rbind(resmean[itr, ],theta0,  sqrt(diag(resvar[,  ,itr]/n))))
}
save(resmean, resvar, tempba, x, lb, file = paste('resmean', n, sdis, sep = '_'))
}}

load(paste('resmean', n, sdis, sep = '_'))
f <- sin(2 * pi * x)
fest <- apply(resmean[, ] %*% t(tempba), 2, median, na.rm = T)

mvf <- fupper <- flower <- fest1 <- matrix(NA, nsim, length(x))
estv <- (cov.rob(resmean)$cov)
for(itr in 1 : nsim){
	print(itr)
    theta <- resmean[itr, ]
fest1[itr, ] <- tempba %*% theta
temp <- tempba
#mvf[itr, ] <- diag(temp %*% resvar[, , itr]%*% t(temp ))/n #
#vf = mvf[itr, ]
#fupper[itr, ] <- fest1[itr, ] + 1.96 * sqrt(vf)
#flower[itr, ] <- fest1[itr, ] - 1.96 * sqrt(vf)

}

# cptheta <- matrix(NA, nsim, lb)
# for(itr in 1 : nsim){
# print(itr)
#     theta <- resmean[itr, ]

# mvf  <- diag(resvar[, , itr]/n)
# #cptheta[itr, ] <- (theta0 <= theta + 1.96 * sqrt(mvf) & theta0 >= theta - 1.96 * sqrt(mvf)  )

# }

# cp <- rep(NA, length(x))
# for(ix in 1:length(x)){

# cp[ix] <- mean(fupper[1:nsim, ix] >= f[ix] & flower[1:nsim, ix] <= f[ix] , na.rm = T)
# }



#fupper <- apply(fupper, 2, median, na.rm = T)
#flower <- apply(flower, 2, median, na.rm = T)
fest <- apply(fest1, 2, median, na.rm = T)
festupl <- apply(fest1, 2, quantile, c(0.025, 0.975), na.rm = T)

pdf('temp.pdf')#pdf(paste(paste('mest', n, sdis, sep = '_'), '.pdf', sep = ''))
plot( f~x, type = 'l', lty = 1, ylim = c(-2.5, 2.5))
lines(fest~x, lty = 2)
#lines(density(x), col = 2)
#lines(x, fupper, lty = 2, col = 2)
#lines(x, flower, lty = 2, col = 2)
lines(festupl[1, ] ~x, lty = 2, col = grey(0.3))
lines(festupl[2, ] ~x, lty = 2, col = grey(0.3))

dev.off()
temp1 <- apply(fest1, 2, sd)
temp2 <- sqrt(apply(((fupper - fest1)/1.96)^2, 2, mean))
pdf('varlap.pdf')
plot(temp2 ~ temp1, ylab= 'estimated variance', xlab = 'emprical variance', main = '')
abline(0, 1)
dev.off()

pdf('mcpunif.pdf')
plot(cp ~x, xlab = 'X', ylab = 'coverage probability', main = '', ylim = c(0.6, 1))
abline(h = 0.95)
dev.off()



# try(spg(theta0, meanest, gr = NULL,  method = 2, lower = theta0 - invAS, upper = theta0 + invAS, project = NULL, projectArgs = NULL, control = list(checkGrad = TRUE, ftol = 1e-7, gtol = 1e-7, M = 20, maxfeval = 1000, eps = 1e-8), quiet = FALSE, alertConvergence=TRUE,  my = my, mw = mw, mc = mc, sdis = sdis))#try(optim(theta0, meanest, gr =  NULL, my, mw, mc, sdis,  method = 'L-BFGS-B', control = list(factr = 1e7, pgtol = 1e-7, parscale = c(0.1, rep(1, 4), 0.1) , ndeps = rep(1e-10, lb)), lower = theta0-  invAS, upper = theta0 +  invAS,  hessian = FALSE))#try(nlm(meanest, theta0, my, mw, mc, sdis))####### try(dfsane(theta0, smeanest, method=3, control=list(), quiet=FALSE, alertConvergence=TRUE, my, mw, mc, sdis))# try(nls.lm(theta0, lower=NULL, upper=NULL, smeanest, jac = NULL,control = nls.lm.control(ftol = 1e-8, ptol = 1e-8, epsfcn = 1e-4),my, mw, mc, sdis))# multiroot(smeanest, theta0, maxiter = 100, rtol = 1e-6, atol = 1e-6, ctol = 1e-6,useFortran = TRUE, positive = FALSE, jacfunc = NULL, jactype = "fullint", verbose = FALSE, bandup = 1, banddown = 1, parms = NULL, my, mw, mc, sdis)



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


sdiff <- matrix(NA, 3, 9)
for(sdis in 0:2){
for(n in seq(200, 1000, 100)){
load(paste('resmean', n, sdis, sep = '_'))
#lb = ceiling(1.3 * n^{1/5})#4 * (n <= 500) + 5 * (n>500 & n <= 800) + 6 * (n >800)#ceiling(1.3 * n^(1/5))
nhb = 1/lb

f <- sin(2 * pi * x)
ix <- 500
fest <- apply(resmean[, ] %*% t(tempba), 2, median, na.rm = T)
f <- sin( 2 * pi * x)
sdiff[sdis+ 1, (n - 200)/100 + 1]<- mean(abs(fest - f)) * sqrt(nhb) 
}
}
size <- sqrt(seq(200, 1000, 100))^{-1}

pdf('temp.pdf')
plot(sdiff[1, ]~ size, type = 'l', ylim = c(min(sdiff), max(sdiff)))
lines(sdiff[2, ] ~size, lty = 4)
lines(sdiff[3, ]~size, lty = 6)
legend('topleft', c('normal error', 'laplace error', 'uniform error'), lty = c(1, 4, 6))
dev.off()


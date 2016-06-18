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
    mx <-   x  * (x - 0.5) * (x - 1) #sin(x * 2  * pi) #dbeta(x, 2, 2)#2* sin((x - 0.5)  * 12) * exp(-((x - 0.5) * 12)^2/10)##1 - 28.76 * x + 251.97 * x^2 - 1028.50 * x^3 + 2069.72 * x^4 
    y <- mx  + rnorm(n, 0, sig) 
    return(cbind(y, w))
}
ddnorm <- function(x){
    1/sqrt(2 * pi * sig^2) * exp(-x^2/(2 *sig^2)) * (-x/sig^2)
}
sdlap <- function(x, m, s){
      1/(2 * s) * exp( -(x- m)/s) * pnorm((x - m)/h, 0, 1)  + 1/(2 * s) * exp( -(m- x)/s) * (1 - pnorm((x - m)/h, 0, 1) )

}
ahmatrix <- function(theta, sdis, bs, ywxc, ywxxc){
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
    
    ma = mb %*% mA %*% ginv(t(mA) %*% mA, 1e-16)#mb %*% ginv(t(mA), tol = 1e-16)## #mb %*% ginv(t(mA), tol = 1e-16)#
   return(ma)
    
}
meanest <- function(theta,  my, mw, mc, sdis){
 #   print(theta)
    ma <- ahmatrix(theta, sdis, bs, ywxc, ywxxc)
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

smeanest <- function(theta,  my, mw, mc, sdis){
 #   print(theta)
# theta <- pnorm(theta) * 4 - 2
    ma <- ahmatrix(theta, sdis, bs, ywxc, ywxxc)
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


ssmeanest <- function(theta,  my, mw, mc, sdis){
 #   print(theta)
# theta <- pnorm(theta) * 4 - 2
    ma <- ahmatrix(theta, sdis, bs, ywxc, ywxxc)
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

a <-  2
b <- 2
#theta0 <- rep(-0.5, lb)
nsim <- 1100
n <- 2000
h <- n^{-1/5}
sdis <- 0

sig <- 0.5 #  sqrt((0.5)^2/2) lap = 1
sigw <- 0.5 #  sqrt((0.5)^2/2)#norm sqrt((0.5)^2/2) lap 1.48
yw<- simu(1000, a, b, sdis)
lly <- 10
llw <- 10
up = 1#sqrt(0.025 * 12)/2
ul = -1#-sqrt(0.025 * 12)/2
u = (up + ul)/2
v = abs(up - ul)/2

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
myGrid <- createNIGrid(dim=2, type=c("GLe", "GLe"), level=c(6, 6))
#lyw <- getNodes(myGrid)
#ly <- log(lyw[, 1]/(1- lyw[, 1]))
#lw <- v * (lyw[, 2] - 1/2) + u
#ly <-  unique(ly)
#lw <- unique(lw)
#Delta <- getWeights(myGrid) * 1/(lyw[, 1] * (1 - lyw[, 1])) * v
gy <- gauss.quad(lly,alpha=0,beta=0)
gw <- gauss.quad(llw)
ly <- log((1 + gy$nodes)/(1 - gy$nodes)) #sqrt( pi) *pnorm(gy$nodes, 0, 1/sqrt(2))
lw <- v * (gw$nodes) + u
Delta <-   as.vector((gy$weights * 2/(1 - gy$nodes^2)) %o%  (gw$weights * v))

}

x <-  rbeta(1000, a, b) 

o <- order(x)
x <- x[o]
lux <- pbeta(x, a, b)
lx <-  qbeta(seq(0.0001, 0.9999, length.out = 15), a, b) ##seq(0.01, 1, length.out = 15)#
cdom <- sum(dbeta(lx, a, b))

olc <- dbeta(lx, a, b)/cdom
ywxx <- expand.grid(ly, lw, lx, lx)
lc <- dbeta(ywxx[, 4], a, b)/cdom
ywxxc <- cbind( ywxx, lc)
ywx <- expand.grid(lx, ly, lw)
lc <- dbeta(ywx[, 1] , a, b)/cdom
ywxc <- cbind(lc, ywx)

knots = c(0.25,  0.5,  0.75)#c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
lb = 6
b2 <- create.bspline.basis(range(x, lx), breaks = c(0,  0.2, 0.8, 1), norder = 4)
basis2 <- eval.basis(ywxxc[, 4], b2)#bs(ywxxc[, 4], 3, knots, Boundary.knots = c(0, 1))
basis1 <-eval.basis(ywxxc[, 3], b2)# bs(ywxxc[, 3], 3, knots, Boundary.knots = c(0, 1))
basis <- eval.basis(ywxc[, 2], b2) # bs(ywxc[, 2], 3, knots, Boundary.knots = c(0, 1))
basisx <-  eval.basis(lx, b2) #bs(lx, 3, knots, Boundary.knots = c(0, 1))
llx <- length(lx)

lly <- length(ly)
llw <- length(lw)




f =  x *  (x - 0.5) * (x - 1)#sin(x * 2  * pi) #dbeta(x, 2, 2)# 2* sin((x - 0.5)  * 12) * exp(-((x - 0.5) * 12)^2/10)#2000 * ((x- 0.5)* 10) - 251.97 * ((x -0.5)*10)^2 - 102 * ((x -0.5)*10)^3 +20* ((x -0.5) * 10)^4 #1 - 28.76 * x + 251.97 * x^2 - 1028.50 * x^3 +2069.72 * x^4 #-1975.95 * x^5 +712.78 * x^6
tempba <-eval.basis(x, b2) #bs(x, 3, knots, Boundary.knots = c(0, 1))
tempf <- function(theta){
    sum((tempba %*% theta - f)^2)
}
theta0 <- optim(rep(0.5, lb), tempf, gr = NULL,  method = 'BFGS',  hessian = FALSE) $par
resmean <- matrix(NA, nsim, lb)
resvar <- array(NA, c(lb, lb, nsim))


for(itr in 1:nsim){
    print(itr)
set.seed(itr + 2017)
ywdata <- simu(n, a, b, sdis)
my <- matrix(ywdata[, 1], n, llx)
mw <- matrix(ywdata[, 2], n, llx)
mlx <- matrix(lx, nrow = n, ncol = llx, byrow = T)
mc <- matrix(olc, nrow = n, ncol = llx, byrow = T)
wt <-  ginv(ssmeanest(theta0, my, mw, mc, sdis), 1e-16)#diag(1, lb, lb)#
#invA <-  ginv(jacobian(smeanest, theta0, method="simple", side=NULL, method.args=list(eps = 1e-4), my, mw, mc, sdis), tol = 1e-16)
#S <- smeanest(theta0, my, mw, mc, sdis)
#invAS <- as.numeric(abs(invA %*% S)) * c(1, rep(1, lb - 1))  * 6
tryres <- optimx(theta0, meanest, gr = NULL, hess = NULL, lower = -Inf, upper = Inf, method = 'BFGS',  itnmax=NULL, hessian=FALSE, control=list(factr = 1e7, pgtol = 1e-7, parscale = c(1e-2, rep(1, 4), 1e-2) , ndeps = rep(1e-7, lb)), my = my, mw = mw, mc = mc, sdis = sdis)
tryres$par <- coef(tryres)
tryres$convergence <- tryres$convcode


print(tryres$convergence)
    if(class(tryres) != 'try-error'){
        resmean[itr, ] <- tryres$par

	jh <- jacobian(smeanest, resmean[itr, ], method="simple", side=NULL, method.args=list(eps =  4e-5), my, mw, mc, sdis)
	sh <- ssmeanest(resmean[itr, ], my, mw, mc, sdis)
	jwh <- ginv(t(jh) %*% wt %*% jh, 1e-16)
	
	resvar[, , itr] <- jwh %*% t(jh) %*% wt %*% sh %*% t(wt) %*% jh %*% t(jwh)
    }




print(rbind(resmean[itr, ],theta0,  sqrt(diag(resvar[,  ,itr]/n))))
}
fest <- apply(resmean[, ] %*% t(tempba), 2, median, na.rm = T)
#fest <- tempba %*% apply(resmean, 2, median, na.rm = T) 



mvf <- fupper <- flower <- fest1 <- matrix(NA, nsim, length(x))
estv <- (cov.rob(resmean)$cov)
for(itr in 1 : nsim){
print(itr)
    theta <- resmean[itr, ]
fest1[itr, ] <- tempba %*% theta
temp <- tempba
mvf[itr, ] <- diag(temp %*% resvar[, , itr]%*% t(temp ))/n #
vf = mvf[itr, ]
fupper[itr, ] <- fest1[itr, ] + 1.96 * sqrt(vf)
flower[itr, ] <- fest1[itr, ] - 1.96 * sqrt(vf)

}

cptheta <- matrix(NA, nsim, lb)
for(itr in 1 : nsim){
print(itr)
    theta <- resmean[itr, ]

mvf  <- diag(resvar[, , itr]/n)
cptheta[itr, ] <- (theta0 <= theta + 1.96 * sqrt(mvf) & theta0 >= theta - 1.96 * sqrt(mvf)  )

}

cp <- rep(NA, length(x))
for(ix in 1:length(x)){

cp[ix] <- mean(fupper[1:nsim, ix] >= f[ix] & flower[1:nsim, ix] <= f[ix] , na.rm = T)
}
fupper <- apply(fupper, 2, median, na.rm = T)
flower <- apply(flower, 2, median, na.rm = T)
fest <- apply(fest1, 2, median, na.rm = T)
festupl <- apply(fest1, 2, quantile, c(0.025, 0.975), na.rm = T)
pdf('est2.pdf')
plot(x, f, type = 'l', lty = 1)
lines(x, fest, lty = 2, col = 2)
#lines(density(x), col = 2)
lines(x, fupper, lty = 2, col = 2)
lines(x, flower, lty = 2, col = 2)
#lines(x, festupl[1, ], lty = 2, col = 3)
#lines(x, festupl[2, ], lty = 2, col = 3)

dev.off()
pdf('cpest.pdf')
plot(cp~x)
dev.off()


pdf('estmean.pdf')
plot(x, f)
points(x, fest, col = 3)
dev.off()


# try(spg(theta0, meanest, gr = NULL,  method = 2, lower = theta0 - invAS, upper = theta0 + invAS, project = NULL, projectArgs = NULL, control = list(checkGrad = TRUE, ftol = 1e-7, gtol = 1e-7, M = 20, maxfeval = 1000, eps = 1e-8), quiet = FALSE, alertConvergence=TRUE,  my = my, mw = mw, mc = mc, sdis = sdis))#try(optim(theta0, meanest, gr =  NULL, my, mw, mc, sdis,  method = 'L-BFGS-B', control = list(factr = 1e7, pgtol = 1e-7, parscale = c(0.1, rep(1, 4), 0.1) , ndeps = rep(1e-10, lb)), lower = theta0-  invAS, upper = theta0 +  invAS,  hessian = FALSE))#try(nlm(meanest, theta0, my, mw, mc, sdis))####### try(dfsane(theta0, smeanest, method=3, control=list(), quiet=FALSE, alertConvergence=TRUE, my, mw, mc, sdis))# try(nls.lm(theta0, lower=NULL, upper=NULL, smeanest, jac = NULL,control = nls.lm.control(ftol = 1e-8, ptol = 1e-8, epsfcn = 1e-4),my, mw, mc, sdis))# multiroot(smeanest, theta0, maxiter = 100, rtol = 1e-6, atol = 1e-6, ctol = 1e-6,useFortran = TRUE, positive = FALSE, jacfunc = NULL, jactype = "fullint", verbose = FALSE, bandup = 1, banddown = 1, parms = NULL, my, mw, mc, sdis)

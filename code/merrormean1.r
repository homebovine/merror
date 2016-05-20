library(splines)
library(smoothmest)
library(optimx)
library(statmod)
library(numDeriv)
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
    mx <- 1 - 28.76 * x + 251.97 * x^2 - 1028.50 * x^3 + 2069.72 * x^4 - 1975.95 * x^5 + 712.78 * x^6
    y <- mx  + rnorm(n, 0, sig) 
    return(cbind(y, w))
}
ddnorm <- function(x){
    1/sqrt(2 * pi * sig^2) * exp(-x^2/(2 *sig^2)) * (-x/sig^2)
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
             deps <- dunif(w - x, ul, up) * dnorm(y - bst, 0, sig) * c * v        }
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
        deps <-  dunif(w - x2, ul, up) *c *  dnorm(y - bst2,  0, sig)/sumL * dunif(w - x1, ul, up) *   dnorm(y - bst1,  0, sig) * v
        ddeps <- dunif(w - x2, ul, up) *c *  ddnorm(y - bst2)/sumL * dunif(w - x1, ul, up) *   dnorm(y - bst1,  0, sig) * v
       
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
    
    ma = mb %*% ginv(t(mA))#mb %*% mA %*% ginv(t(mA) %*% mA)# #
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
    nums <- deps %*% basisx
    S <- sum(apply(-nums/dom - numa/dom, 2, mean)^2)
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
    nums <- deps %*% basisx
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
    nums <- deps %*% basisx
    S <- -nums/dom - numa/dom
    mS <- t(S) %*% S /n
}

a <-  1
b <- 1
#theta0 <- rep(-0.5, lb)
nsim <- 1100
n <- 400
sdis <- 2
sig <- 0.5
sigw <- sqrt((0.5)^2/2)##normal 0.5
yw<- simu(1000, a, b, sdis)
lly <- 10
llw <- 10
up = 1#sqrt(0.025 * 12)/2
ul = -1#-sqrt(0.025 * 12)/2
u = (up + 1 + ul)/2
v = abs(up + 1 - u)

if(sdis == 0 | sdis == 1){
gy <- gauss.quad(lly,kind="hermite",alpha=0,beta=0)
gw <- gauss.quad(llw,kind="hermite",alpha=0,beta=0)
ly <- sqrt( pi) *pnorm(gy$nodes, 0, 1/sqrt(2))
lw <- sqrt( pi) *pnorm(gw$nodes, 0, 1/sqrt(2))
}else{
gy <- gauss.quad(lly,kind="hermite",alpha=0,beta=0)
gw <- gauss.quad(llw)
ly <- sqrt( pi) *pnorm(gy$nodes, 0, 1/sqrt(2))
lw <- v * (gw$nodes) + u
}
Delta <-   as.vector(gy$weights %o%  gw$weights)
lx <-  qbeta(seq(0, 1, length.out = 20), a, b) 
cdom <- sum(dbeta(lx, a, b))
olc <- dbeta(lx, a, b)/cdom
ywxx <- expand.grid(ly, lw, lx, lx)
lc <- dbeta(ywxx[, 4], a, b)/cdom
ywxxc <- cbind( ywxx, lc)
ywx <- expand.grid(lx, ly, lw)
lc <- dbeta(ywx[, 1] , a, b)/cdom
ywxc <- cbind(lc, ywx)

knots = c(0.25,  0.5,  0.75)#c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
lb = length(knots) + 3
basis2 <- bs(ywxxc[, 3], 3, knots, Boundary.knots = c(0, 1))
basis1 <- bs(ywxxc[, 4], 3, knots, Boundary.knots = c(0, 1))
basis <-  bs(ywxc[, 2], 3, knots, Boundary.knots = c(0, 1))
basisx <-  bs(lx, 3, knots, Boundary.knots = c(0, 1))
llx <- length(lx)

lly <- length(ly)
llw <- length(lw)


x <-  rbeta(1000, a, b)  

f = 1 - 28.76 * x + 251.97 * x^2 - 1028.50 * x^3 +2069.72 * x^4 -1975.95 * x^5 +712.78 * x^6
tempba <- bs(x, 3, knots, Boundary.knots = c(0, 1))
tempf <- function(theta){
    sum((tempba %*% theta - f)^2)
}
theta0 <- optim(rep(0.5, lb), tempf, gr = NULL,  method = 'BFGS',  hessian = FALSE) $par
resmean <- matrix(NA, nsim, lb)
resvar <- array(NA, c(lb, lb, nsim))
for(itr in 84:nsim){
    print(itr)
set.seed(itr + 2015)
ywdata <- simu(n, a, b, sdis)
my <- matrix(ywdata[, 1], n, llx)
mw <- matrix(ywdata[, 2], n, llx)
mlx <- matrix(lx, nrow = n, ncol = llx, byrow = T)
mc <- matrix(olc, nrow = n, ncol = llx, byrow = T)
tryres <- try(nls.lm(theta0, lower=NULL, upper=NULL, smeanest, jac = NULL,control = nls.lm.control(ftol = 1e-8, ptol = 1e-8, epsfcn = 1e-4),my, mw, mc, sdis))# try(dfsane(theta0, smeanest, method=2, control=list(), quiet=FALSE, alertConvergence=TRUE, my, mw, mc, sdis))# try(optim(theta0 * 1.01, meanest, gr =  NULL, my, mw, mc, sdis,  method = 'L-BFGS-B', control = list(), lower = -Inf, upper = Inf,  hessian = FALSE))
    if(class(tryres) != 'try-error'){
        resmean[itr, ] <- tryres$par
	invh <- ginv(jacobian(smeanest, resmean[itr, ], method="simple", side=NULL, method.args=list(eps = 1e-4), my, mw, mc, sdis))
	resvar[, , itr] <- invh %*% (ssmeanest(resmean[itr, ], my, mw, mc, sdis)/n) %*% t(invh)
    }




                                        #spg(theta0 * 1.11, meanest, gr = smeanest,  method = 2, lower = -20, upper = 20, project = NULL, projectArgs = NULL, control = list(checkGrad = FALSE), quiet = FALSE, alertConvergence=TRUE,  my, mw, mc, sdis)$par
print(resmean[itr, ])
}

fest <- apply(resmean[, ] %*% t(tempba), 2, median, na.rm = T)
#fest <- tempba %*% apply(resmean, 2, median, na.rm = T) 

pdf('estmean.pdf')
plot(x, f)
points(x, fest, col = 3)
dev.off()



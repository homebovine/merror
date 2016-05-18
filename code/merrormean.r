library(splines)
library(smoothmest)
library(optimx)
library(statmod)
simu <- function(n,a, b,  sdis){
    x <- 2 * rbeta(n, a, b)  - 1
    if(sdis == 0){
        eps <- rnorm(n, 0, 1)
    }else if(sdis == 1){
        eps <- rdoublex(n, 0, 1)
    }else{
        eps <- runif(n, -1, 1)
    }
    w <- x + eps
    y <- x^3 + cos(x) + rnorm(n, 0, 1)
    return(cbind(y, w))
}
ddnorm <- function(x){
    1/sqrt(2 * pi) * exp(-x^2/2) * (-x)
}
ahmatrix <- function(theta, sdis, bs, ywxc, ywxxc){
    bst1<- (basis1 %*% theta)
    bst2 <- (basis2 %*% theta)
    bst <- (basis %*% theta)
    suml <- function(c, x, y, w,  sdis){
        if(sdis ==0 ){
            deps <- dnorm(w - x) * dnorm(y - bst, 0, 1) * c
        }else if (sdis == 1 ){
            deps <- ddoublex(w - x) * dnorm(y - bst, 0, 1) * c
        }else {
             deps <- ddoublex(w - x) * dnorm(y - bst, 0, 1) * c
        }
    }
    sumL <- suml(ywxc[, 1], ywxc[, 2], ywxc[, 3], ywxc[, 4],  sdis)
    sumL <- apply(matrix(sumL, ncol = llx, byrow = T), 1, sum)
    sumL <- rep(sumL,  llx^2)
    integrand1 <- function(x1, x2,c,  y, w, sdis){
    if(sdis == 0 ){
        deps <-  dnorm(w - x2) *c *  dnorm(y - bst2,  0, 1)/sumL * dnorm(w - x2) *   dnorm(y - bst1,  0, 1)
        ddeps <-  dnorm(w - x2) *c *  ddnorm(y - bst2)/sumL * dnorm(w - x2) *   dnorm(y - bst1,  0, 1)
    }
    else if(sdix == 1){
         deps <-  ddoublex(w - x2) *c *  dnorm(y - bst2,  0, 1)/sumL * ddoublex(w - x2) *   dnorm(y - bst1,  0, 1)
        ddeps <- ddoublex(w - x2) *c *  ddnorm(y - bst2)/sumL1 * ddoublex(w - x2) *   dnorm(y - bst1,  0, 1)
        
    }else{
        deps <-  dunif(w - x2, -1, 1) *c *  dnorm(y - bst2,  0, 1)/sumL * dunif(w - x2, -1, 1) *   dnorm(y - bst1,  0, 1)
        ddeps <- dunif(w - x2, -1, 1) *c *  ddnorm(y - bst2)/sumL * dunif(w - x2, -1, 1) *   dnorm(y - bst1,  0, 1)
       
    }
    return(cbind(deps, ddeps))
}
    temp <- integrand1(ywxxc[, 3], ywxxc[, 4], ywxxc[, 5], ywxxc[, 1], ywxxc[, 2], sdis)
    mA <- apply(matrix(temp[, 1], ncol = lly * llw, byrow = T), 1, sum) * Delta
    mH <- apply(matrix(temp[, 2], ncol = lly * llw, byrow = T), 1, sum) * Delta
    mA <- matrix(mA,  llx, llx)
    mH <- -matrix(mH,  llx, llx)
    mb <- matrix(NA, lb, llx)
    for (i in 1: llx){
        mb[, i] <- mH[i, ] %*% basisx
    }
    ma = mb %*% t(ginv(mA) )
   return(ma)
    
}
meanest <- function(theta,  my, mw, mc, sdis){
 #   print(theta)
    ma <- ahmatrix(theta, sdis, bs, ywxc, ywxxc)
    bst <- matrix(basisx %*% theta, n, llx, byrow = T)
    if(sdis == 0 ){
        deps <-  dnorm(mw - mlx) *  dnorm(my - bst,  0, 1) * mc
        ddeps <-  dnorm(mw -mlx) *  ddnorm(my - bst) * mc
       
    }else if(sdix == 1){
          deps <-  ddoublex(mw - mlx) *  dnorm(my - bst,  0, 1) * mc
          ddeps <-  ddoublex(mw - mlx) *  ddnorm(my - bst) * mc
         
        
    }else{
        deps <-   dunif(mw - mlx, -1, 1) *  dnorm(my - bst,  0, 1) * mc
        ddeps <-   dunif(mw - mlx, -1, 1) *  ddnorm(my - bst) * mc
       
    }
    numa <- deps %*% t(ma)
    dom <- matrix(apply(deps, 1, sum), ncol = lb, nrow = n)
    nums <- deps %*% basisx
    S <- sum(apply(-nums/dom - numa/dom, 2, mean)^2)
}

smeanest <- function(theta,  my, mw, mc, sdis){
 #   print(theta)
# theta <- pnorm(theta) * 4 - 2
    ma <- ahmatrix(theta, sdis, bs, ywxc, ywxxc)
    bst <- matrix(basisx %*% theta, n, llx, byrow = T)
    if(sdis == 0 ){
        deps <-  dnorm(mw - mlx) *  dnorm(my - bst,  0, 1) * mc
        ddeps <-  dnorm(mw -mlx) *  ddnorm(my - bst) * mc
       
    }else if(sdix == 1){
          deps <-  ddoublex(mw - mlx) *  dnorm(my - bst,  0, 1) * mc
          ddeps <-  ddoublex(mw - mlx) *  ddnorm(my - bst) * mc
         
        
    }else{
        deps <-   dunif(mw - mlx, -1, 1) *  dnorm(my - bst,  0, 1) * mc
        ddeps <-   dunif(mw - mlx, -1, 1) *  ddnorm(my - bst) * mc
       
    }
    numa <- deps %*% t(ma)
    dom <- matrix(apply(deps, 1, sum), ncol = lb, nrow = n)
    nums <- deps %*% basisx
    S <- apply(-nums/dom - numa/dom, 2, mean)
}

a <-  2
b <- 2
theta0 <- rep(-0.5, lb)
nsim <- 1100
n <- 400
sdis <- 0
yw<- simu(1000, a, b, sdis)
lly <- 15
llw <- 15
ly <-  seq(min(yw[, 1]), max(yw[, 1]), length.out = lly)
lw <- seq(qnorm(0.01), qnorm(0.99), length.out = llw)
Delta <-  diff(ly)[1] * diff(lw)[1]
lx <- 2* qbeta(seq(0.01, 0.99, length.out = 15), a, b) - 1
lc <- dbeta((lx + 1)/2, a, b)/sum(dbeta((lx + 1)/2, a, b))
ywxx <- expand.grid(ly, lw, lx, lx)
ywxxc <- cbind( ywxx, rep(lc, each = nrow(ywxx)/15))
ywx <- expand.grid(lx, ly, lw)
ywxc <- cbind(rep(lc, nrow(ywx)/15), ywx)

knots = c(-0.5,   0,  0.5)#c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
lb = length(knots) + 3
basis2 <- bs(ywxxc[, 3], 3, knots, Boundary.knots = c(-1, 1))
basis1 <- bs(ywxxc[, 4], 3, knots, Boundary.knots = c(-1, 1))
basis <-  bs(ywxc[, 2], 3, knots, Boundary.knots = c(-1, 1))
basisx <-  bs(lx, 3, knots, Boundary.knots = c(-1, 1))
llx <- length(lx)

lly <- length(ly)
llw <- length(lw)



resmean <- matrix(NA, nsim, lb)
for(itr in 1:nsim){
    print(itr)
set.seed(itr)
ywdata <- simu(n, a, b, sdis)
my <- matrix(ywdata[, 1], n, llx)
mw <- matrix(ywdata[, 2], n, llx)
mlx <- matrix(lx, nrow = n, ncol = llx, byrow = T)
mc <- matrix(lc, nrow = n, ncol = llx, byrow = T)
tryres <-  try(optim(theta0 * 1.01, meanest, gr = smeanest, my, mw, mc, sdis,  method = 'BFGS',  hessian = FALSE))
    if(class(tryres) != 'try-error'){
        resmean[itr, ] <- tryres$par
    }
    
                                        #spg(theta0, meanest, gr = NULL,  method = 2, lower = -20, upper = 20, project = NULL, projectArgs = NULL, control = list(), quiet = FALSE, alertConvergence=TRUE,  my, mw, mc, sdis)$par
}


x <- 2 * rbeta(1000, a, b)  - 1
f = x^3 + cos(x)
tempba <- bs(x, 3, knots, Boundary.knots = c(-1, 1))
tempf <- function(theta){
    sum((tempba %*% theta - f)^2)
}
fest <- apply(resmean %*% t(tempba), 2, median, na.rm = T)
fest <- tempba %*% apply(resmean, 2, median, na.rm = T) 
 theta0 <- optim(theta0, tempf, gr = NULL,  method = 'BFGS',  hessian = FALSE) $par
pdf('estmean.pdf')
plot(x, f)
points(x, fest, col = 3)
dev.off()



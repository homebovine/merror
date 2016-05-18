library(splines)
library(smoothmest)
library(optimx)
library(statmod)
library(numDeriv)
library(minpack.lm)
library(fda)
library(msm)
simu <- function(n,a, b,  sdis){
    x <- rbeta(n, a, b) 
    if(sdis == 0){
        eps <- rnorm(n, 0, sig)
    }else if(sdis == 1){
        eps <- rdoublex(n, 0, sig)
    }else{
        eps <- runif(n, -sqrt(3), sqrt(3))
    }
    y <- x + eps
    return(y)
}

denest <- function(theta, sdis, bs, y){
    bst<- matrix(exp(basis %*% theta), n, ng, byrow = T)
    integrand1 <- function(x, y, sdis){
    if(sdis == 0 ){
        deps <-  bst * dnorm(y - x,  0, sig) * weights
    }
    else if(sdis == 1){
        deps <- bst * ddoublex(y - x, 0, sig) * weights
    }else{
        deps <- bst * dunif(y - x, -sqrt(3), sqrt(3)) * weights
    }
}
    
    
    intres1 <- apply(integrand1(nodes, y, sdis), 1, sum)# apply(sapply(1 : 10, temp), 1, sum)
    intres2 <- sum(bst[1, ] * weights[1, ])

    
    object <- (-mean(log(intres1) ) +   log(intres2))
    
}

derivdenest <- function(theta, sdis, bs, y){
	       # bst<- matrix(((dbeta(nodes, a, b))), n, ng)#matrix(exp(basis %*% theta), n, ng, byrow = T)
    bst<- matrix(exp(basis %*% theta), n, ng, byrow = T)
    integrand1 <- function(x, y, sdis){
    if(sdis == 0 ){
        deps <-  bst * dnorm(y - x,  0, sig) * weights
    }
    else if(sdis == 1){
        deps <- bst * ddoublex(y - x, 0, sig) * weights
    }else{
        deps <- bst * dunif(y - x, -sqrt(3), sqrt(3)) * weights
    }
    ddeps <- deps %*% basis
    return(list(deps, ddeps))
}
  
    
    temp <- integrand1(nodes, y, sdis)
deps <- temp[[1]]	      
ddeps <- temp[[2]]
    
    intres1 <- matrix(apply(deps, 1, sum), n, lb)# apply(sapply(1 : 10, temp), 1, sum)
    intres1num <- ddeps/intres1
    temp <- bst[1, ] * weights[1, ]
    intres2 <- sum(temp) 
    intres2num <- matrix(t(temp) %*% basis, n, lb, byrow = T)
    temp <- t(intres1num - intres2num/intres2)


    (temp %*% t(temp))/n

    
    
}

score <- function(theta, sdis, bs, y){
    bst<- matrix(exp(basis %*% theta), n, ng, byrow = T)
    integrand1 <- function(x, y, sdis){
    if(sdis == 0 ){
        deps <-  bst * dnorm(y - x,  0, sig) * weights
    }
    else if(sdis == 1){
        deps <- bst * ddoublex(y - x, 0, sig) * weights
    }else{
        deps <- bst * dunif(y - x, -sqrt(3), sqrt(3)) * weights
    }
    ddeps <- deps %*% basis
    return(list(deps, ddeps))
}
 
    
    temp <- integrand1(nodes, y, sdis)
    deps <- temp[[1]]	      
ddeps <- temp[[2]]
    
    intres1 <- matrix(apply(deps, 1, sum), n, lb)# apply(sapply(1 : 10, temp), 1, sum)
    intres1num <- ddeps/intres1
    temp <- bst[1, ] * weights[1, ]
    intres2 <- sum(temp) 
    intres2num <- matrix(t(temp) %*% basis, n, lb, byrow = T)
    temp <- (intres1num - intres2num/intres2)
    -apply(temp, 2, mean)

    
    
}

hessian <- function(theta, sdis, bs, y){
    bst<- matrix(exp(basis %*% theta), n, ng, byrow = T)
    integrand1 <- function(x, y, sdis){
    if(sdis == 0 ){
        deps <-  bst * dnorm(y - x,  0, sig) * weights
    }	
    else if(sdis == 1){
        deps <- bst * ddoublex(y - x, 0, sig) * weights
    }else{
        deps <- bst * dunif(y - x, -sqrt(3), sqrt(3)) * weights
    }
    BBnum <- BBnum1 <- array(NA, c(lb, lb,  n))
    for(i in 1:n){
    	   dnom <- sum(deps[i, ])
    	   BBnum[, , i] <- t(diag(deps[i, ]) %*% basis) %*% (basis)/dnom
	   temp <- (t(deps[i, ]) %*% basis)
	   BBnum1[, , i] <-  t(temp) %*% (temp)/dnom^2
    }   
    BBnum <- apply(BBnum, c(1, 2), mean)
    BBnum1 <- apply(BBnum1, c(1, 2), mean)
    return(list(BBnum, BBnum1))
}
 
    
    temp <- integrand1(nodes, y, sdis)
    R1 <- temp[[1]]	      
    R2<- temp[[2]]
    
    
    temp <- bst[1, ] * weights[1, ]
    intres2 <- sum(temp) 
    intres2num <- t(diag(temp) %*% basis) %*% basis
    R3 <- intres2num/intres2
    R4 <- t(temp) %*% basis
    R4 <- t(R4) %*% (R4)/intres2^2
    
   (R1 -R2 -  R3  +  R4)

    
    
}



a <- 4
b <- 4
sig <- sqrt(1/2)
nsim <- 1100
n <- 400
sdis <- 2
ng <-  30
t <- rbeta(10000, a, b)
lb <- 6

knots = c(0.3, 0.5, 0.7)#c( 0.2, 0.4,  0.6, 0.8)#c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
nh <- 0.5



#theta0 <- rep(-0.5, lb)
rtheta <- matrix(NA, nsim, lb)
rvar <- matrix(NA, nsim, lb)
wn <- gauss.quad(ng,kind="legendre",alpha=0,beta=0)
#wn$nodes <- seq(-1, 1, length.out = ng)
#wn$weights <- diff(c(-1, wn$nodes))#wn$weights * 1 /sqrt(1 - wn$nodes^2)#
wn$nodes <- 1/2 * wn$nodes + 1/2
objbasis <- create.bspline.basis(c(0, 1),  nbasis = lb,  norder = 4)
basis <- eval.basis(wn$nodes, objbasis)#ns(wn$nodes, df = lb,     Boundary.knots= range(0, 1))
wn$weights <- 1/2 * wn$weights
#basis <- predict(basis, wn$nodes)
nodes <-  matrix(wn$nodes, n, ng, byrow = T) 
weights <-  matrix(wn$weights, n, ng, byrow = T)
covflg <- rep(NA, nsim)

o <- order(t)
t <- t[o]
x <- t 
f <-   dbeta(t, a, b)
getini <- function(theta){
sum((exp(eval.basis(x, objbasis) %*% theta)/(sum(exp(basis %*% theta) * weights[1, ]))  - f)^2)
}
    #theta0 <- optim(rep(1, lb), getini, gr = NULL,  method = 'Nelder-Mead', control = list(maxit = 500),   hessian = TRUE)$par * 1.01
   theta0 <- optim(rep(1, lb), getini, gr = NULL,  method = 'BFGS', control = list(maxit = 500),   hessian = TRUE)$par 
for(itr in 1 : nsim){
    print(itr)
    y <- simu(n, a, b, sdis)
    my <- matrix(y, n, ng)
    
#    res <- spg(theta0, denest, gr=score, method=3, lower=rep(-Inf, lb), upper=rep(Inf, lb), project=NULL, projectArgs=NULL, control=list(maxit = 500), quiet=FALSE, alertConvergence=TRUE, sdis = sdis, bs = basis, y = my)


    res <- try(optim(theta0, denest, gr = score, sdis, basis, my, method = 'BFGS',  control = list(ndeps = rep(1e-6,lb),  reltol = 1e-8, abstol = 1e-8,    maxit = 500), hessian = TRUE))
    temp <- try(ginv(hessian(res$par, sdis, basis, my), tol = 1e-4))#ginv(res$hessian)#
	
print(res$par)
	if(class(res) != 'try-error' && class(temp) != 'try-error'){
sand <- derivdenest(res$par, sdis, basis, my)

    rtheta[itr, ]<- res$par
    S<- eigen(sand)
    eig <- eigen(sand)
    eig$value[eig$value<1e-5] <- 0
    eig$value[eig$value>0] <- eig$value[eig$value>0]^{-1}
    
    rvar[ itr, ] <- diag(eig$vector %*% diag(eig$value) %*% t(eig$vector)/n)#diag(ginv(sand, tol= 1e-2)/n)# diag( temp/n)#diag((temp %*% sand %*% t(temp))/n)
    covflg[itr] <- res$convergence
}
    
}




mtheta <- apply(rtheta, 2, median)
fest1 <- matrix(NA, nsim, length(x))
for(itr in 1 : nsim){
    theta <- rtheta[itr, ]
fest1[itr, ] <- exp(eval.basis( x, objbasis) %*% theta) /sum(exp(basis %*% theta)  * weights[1, ])
}
fest <- apply(fest1, 2, median)
pdf('est1.pdf')
plot(x, f)
lines(x, fest)
lines(density(x), col = 2)
dev.off()

fest  <- exp(predict(basis, x) %*% mtheta) /sum(exp(basis %*% mtheta)  * weights[1, ])


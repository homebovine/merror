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
        eps <- runif(n, -upper, upper)
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
        deps <- bst * dunif(y - x, -upper, upper) * weights
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
        deps <- bst * dunif(y - x, -upper, upper) * weights
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
        deps <- bst * dunif(y - x, -upper, upper) * weights
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
        deps <- bst * dunif(y - x, -upper, upper) * weights
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
sdis <- 0
nsim <- 100
ng <-  30
ntest <- 10

load('framdata')
sig <- sigw#0.0646
n <- ceiling(nrow(framdata))
lb <- 7#round(1.3* n^(1/5))
testint <- gauss.quad(ntest,'hermite', alpha=0,beta=0)
testnode <- sqrt(2) * sig * testint$nodes
testweight <- testint$weights /sqrt(pi)
testnode <- matrix(rep(testnode, 512), ncol = ntest, byrow = T)
testweight <- matrix(rep(testweight, 512), ncol = ntest, byrow = T)
wn <- gauss.quad(ng,kind="legendre",alpha=0,beta=0)
#wn$nodes <- seq(-1, 1, length.out = ng)
#wn$weights <- diff(c(-1, wn$nodes))#wn$weights * 1 /sqrt(1 - wn$nodes^2)#
wn$nodes <- 1/2 * wn$nodes + 1/2
rnodes <- range(wn$nodes)
objbasis <- create.bspline.basis(c(0, 1), nbasis = lb,  norder = 4)
basis <- eval.basis(wn$nodes, objbasis)#ns(wn$nodes, df = lb,     Boundary.knots= range(0, 1))
wn$weights <- 1/2 * wn$weights
#basis <- predict(basis, wn$nodes)
nodes <-  matrix(wn$nodes, n, ng, byrow = T) 
weights <-  matrix(wn$weights, n, ng, byrow = T)
dfram <- density(framdata[, 2], from = 0, to = 1)
t <- dfram$x


x <- t
f <-  dfram$y
rtheta <- matrix(NA, nsim, lb)
getini <- function(theta){
mean((exp(eval.basis(x, objbasis) %*% theta)/(sum(exp(basis %*% theta) * weights[1, ]))  - f)^2)
}

#theta0 <- optim(rep(1, lb), getini, gr = NULL,  method = 'BFGS')$par 
#theta0 <- c(-0.3736096, 0.9683432,   2.0193998,   3.5105639,   4.1407997,  -4.2655174)#c(0.7813047, -0.1409199, 2.2438951, 2.4553044, 4.3708748, 2.5468120, -5.2573229)
estden <- matrix(NA, nsim, 512)
for(itr in 1 : nsim){
	set.seed(2017 + itr)
    train<- sample(1:nrow(framdata), n, replace = T)
    
    y <- framdata[train, 2]
    testy <- framdata[-train, 2]
    my <- matrix(y, n, ng)
    res <- try(optimx(theta0, denest, gr = NULL,  hess = NULL, lower = NULL, upper = NULL,  method = 'newuoa',itnmax=1000, hessian=FALSE, control=list(),  sdis = sdis, bs = basis,y = my))

    rtheta[itr, ]<- coef(res)
    # testxy <- density(testy)
    # testy <- testxy$x
    # testden<- testxy$y
    # dnom<- sum(exp(basis %*% rtheta[itr, ])  * weights[1, ])

    
    # testset  <- testy - testnode
    # rindx <- (testset<0 |testset>1)
    # testset[rindx] <- 0
    # estden[itr, ] <- apply(exp(matrix(eval.basis(as.vector(testset) , objbasis) %*% rtheta[itr, ], ncol = ntest)) * (1 - rindx) * testweight, 1, sum)/dnom
    

    
print(rbind(rtheta[itr, ], theta0))    
}

save( rtheta,  file = paste('realfull1'), precheck = 0)


#load(paste('den', n, sdis, sep = '_'))
mtheta <- apply(rtheta, 2, median)
fupper <- flower <- fest1 <- matrix(NA, nsim, length(x))
dfun <- function(theta){
exp(eval.basis( x, objbasis) %*% theta) /sum(exp(basis %*% theta)  * weights[1, ])
}
for(itr in 1 : nsim){
    theta <- rtheta[itr, ]
fest1[itr, ] <- exp(eval.basis( x, objbasis) %*% theta) /sum(exp(basis %*% theta)  * weights[1, ])
temp <- jacobian(dfun, theta, method = 'simple') # (diag(as.vector(exp(eval.basis( x, objbasis) %*% theta))) %*%  eval.basis( x, objbasis)) /(sum(exp(basis %*% theta)  * weights[1, ]))^2 - exp(eval.basis( x, objbasis) %*% theta) %*% t(t(basis) %*% (exp(basis %*% theta)  * weights[1, ]))  /(sum(exp(basis %*% theta)  * weights[1, ]))^2
vf <- diag(temp %*% ((rvar[, , itr])) %*% t(temp )) #
fupper[itr, ] <- fest1[itr, ] + 1.96 * sqrt(vf)
flower[itr, ] <- fest1[itr, ] - 1.96 * sqrt(vf)

}

nl <- seq(200,2000, 100)
l <- length(nl)
ix <- 500
nhn <- sdiff <- matrix(NA, 3, l)

for(sdis in 0:2){
for(i in 1:l){

pdf('den.pdf')#pdf(paste(paste('den_diff', sep = ''), ".pdf", sep = ''))
plot(sdiff1[1, -(1:3) ] ~ nvec, type = 'l', ylim = c(min(sdiff1), max(sdiff1)+ 5),  ylab = 'Root nh maximum absolute error', xlab ='Sample size')
lines(sdiff1[2,-(1:3)] ~ nvec, lty = 3)
lines(sdiff1[3, -(1:3)] ~ nvec, lty = 6)
legend('topleft',c('normal error', 'laplace error', 'uniform error'),  lty = c(1, 3, 6) )

dev.off()
estf <- exp(eval.basis( testy, objbasis) %*% theta) /sum(exp(basis %*% theta)  * weights[1, ])
plot(estf~x, type = 'l', lty = 2, ylab = 'estimated density')
lines(density(y))

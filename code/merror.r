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
nsim <- 200
ng <-  30
upper <- 1/8
sig <- 0.05
for(sdis in 0){
for(n in c(500, 1000, 2000)){
a <- 4
b <- 4
if(sdis == 1){
sig <- 0.05/sqrt(2) #sqrt(1/3)#norm 1
}else{
sig <- 0.05
}

lb <- round(1.3* n^(1/5))
        #lb = 4 * (n <= 500) + 5 * (n>500 & n <= 800) + 6 * (n >800 & n <= 1000)+ 8 * (n >1000 & n <= 1500)  + 9 * (n > 1500)#

knots = c(0.3, 0.5, 0.7)#c( 0.2, 0.4,  0.6, 0.8)#c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
nh <- 0.5



#theta0 <- rep(-0.5, lb)
rtheta <- matrix(NA, nsim, lb)
rvar <- array(NA, c(lb, lb, nsim))#matrix(NA, nsim, lb)
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
t <- rbeta(n, a, b)
covflg <- rep(NA, nsim)

o <- order(t)
t <- t[o]
x <- t 
f <-   dbeta(t, a, b)
getini <- function(theta){
mean((exp(eval.basis(x, objbasis) %*% theta)/(sum(exp(basis %*% theta) * weights[1, ]))  - f)^2)
}

   theta0 <- optim(rep(1, lb), getini, gr = NULL,  method = 'BFGS')$par 
for(itr in 1 : nsim){
    print(itr)
    y <- simu(n, a, b, sdis)
    my <- matrix(y, n, ng)
    


    res <- try(optimx(theta0, denest, gr = NULL,  hess = NULL, lower = -Inf, upper = Inf,  method = 'BFGS',itnmax=500, hessian=FALSE, control=list(),  sdis = sdis, bs = basis,y = my))#,  control = list(ndeps = rep(1e-6,lb),  reltol = 1e-8, abstol = 1e-8,    maxit = 500), hessian = TRUE))
    #temp <- try(ginv(hessian(res$par, sdis, basis, my), tol = 1e-4))#ginv(res$hessian)#
	
#print(res$par)
#	if(class(res) != 'try-error' && class(temp) != 'try-error'){
#sand <- derivdenest(res$par, sdis, basis, my)

    rtheta[itr, ]<- coef(res)#res$par
 #   S<- eigen(sand)
  #  eig <- eigen(sand)
   # eig$value[eig$value<1e-5] <- 0 ##unif 8e-6
    #eig$value[eig$value>0] <- eig$value[eig$value>0]^{-1}
    
    rvar[, , itr ] <- 0#(eig$vector %*% diag(eig$value) %*% t(eig$vector)/n)#diag(ginv(sand, tol= 1e-2)/n)# diag( temp/n)#diag((temp %*% sand %*% t(temp))/n)
    #covflg[itr] <- res$convergence
#}
print(rbind(rtheta[itr, ], theta0))    
}

save(rvar, rtheta, x, objbasis, basis, weights, lb,   file = paste('densmall', n, sdis, sep = '_'), precheck = 0)
}
}

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
cp <- rep(NA, n)
for(ix in 1:n){

cp[ix] <- mean(fupper[, ix] > dbeta(x[ix], a, b) & flower[, ix] <= dbeta(x[ix], a, b) )
}
fupper <- apply(fupper, 2, median, na.rm = T)
flower <- apply(flower, 2, median, na.rm  = T)
fest <- apply(fest1, 2, median, na.rm = T)
festupl <- apply(fest1, 2, quantile, c(0.025, 0.975), na.rm = T)
pdf(paste(paste('den', n, sdis, sep = '_'), '.pdf', sep = ''))
plot(x, f, type = 'l', lty = 1, ylim = c(-0.5, 5))
lines(x, fest, lty = 2)
#lines(density(x), col = 2)
#lines(x, fupper, lty = 2, col = 2)
#lines(x, flower, lty = 2, col = 2)
lines(x, festupl[1, ], lty = 2, col = grey(0.5))
lines(x, festupl[2, ], lty = 2, col = grey(0.5))

dev.off()
temp1 <- apply(fest1, 2, sd)
temp2 <- sqrt(apply(((fupper - fest1)/1.96)^2, 2, mean))
pdf('varlap.pdf')
plot(temp2 ~ temp1, ylab= 'estimated variance', xlab = 'emprical variance', main = '')
abline(0, 1)
dev.off()

pdf('cpunif.pdf')
plot(cp ~x, xlab = 'X', ylab = 'coverage probability', main = '', ylim = c(0.6, 1))
abline(h = 0.95)
dev.off()

fest  <- exp(predict(basis, x) %*% mtheta) /sum(exp(basis %*% mtheta)  * weights[1, ])

emprical uniform  1.7686831 0.8835746 1.0391104 0.9973378 0.8165267 1.7591435
estimated unifrom 1.614624 1.023110 1.189571 1.178298 1.031448 1.621196


emprical norm 1.1208815 0.9302755 1.9416143 1.9962105 0.9787644 1.1324871
estimated norm 1.708348 1.452574 2.825442 2.827463 1.484432 1.660037

estimated lap 1.701376 1.262628 1.965442 1.977015 1.291913 1.700157
emprical lap 1.641232 1.500220 2.092283 1.753365 1.371908 1.800199
nl <- seq(200,2000, 100)
l <- length(nl)
ix <- 500
nhn <- sdiff <- matrix(NA, 3, l)

for(sdis in 0:2){
for(i in 1:l){
n <- nl[i]

fest1 <- matrix(NA, nsim , 100)
load(paste('den', n, sdis, sep = '_'))
#        lb = 4 * (n <= 500) + 6 * (n>500 & n <= 800) + 7 * (n >800 & n <= 1000)+ 8 * (n >1000 & n <= 1500)  + 9 * (n > 1500)#lb <- round(1.3* n^(1/5)) 
x <- seq(0, 1, length.out = 100)
  nh <- 1/(lb)
for(itr in 1 : nsim){
    theta <- rtheta[itr, ]
f <-   dbeta(x, a, b)
fest1[itr, ] <-abs(exp(eval.basis( x, objbasis) %*% theta) /sum(exp(basis %*% theta)  * weights[1, ]) -f)


}
fest <- apply(fest1,1, max, na.rm = T)

sdiff[sdis + 1, i] = median(fest)# * sqrt(n * nh)#mean(abs(fest - f)) * sqrt(nh)
nhn[sdis+ 1, i] <- sqrt(n* nh)
}
}
nvec <- seq(500, 2000, 100)
sdiff1 <- sdiff * nhn
nhn1 <- nvec
sdiff1 <- sdiff1
pdf('den.pdf')#pdf(paste(paste('den_diff', sep = ''), ".pdf", sep = ''))
plot(sdiff1[1, -(1:3) ] ~ nvec, type = 'l', ylim = c(min(sdiff1), max(sdiff1)+ 5),  ylab = 'Root nh maximum absolute error', xlab ='Sample size')
lines(sdiff1[2,-(1:3)] ~ nvec, lty = 3)
lines(sdiff1[3, -(1:3)] ~ nvec, lty = 6)
legend('topleft',c('normal error', 'laplace error', 'uniform error'),  lty = c(1, 3, 6) )

dev.off()
mse <- matrix(NA, 3, 3)
for(sdis in 0:2){
j <- 0
    for(n in c(500, 1000, 2000)){
j <- j + 1
    fest1 <- matrix(NA, nsim , 100)
    load(paste('densmall', n, sdis, sep = '_'))
#        lb = 4 * (n <= 500) + 6 * (n>500 & n <= 800) + 7 * (n >800 & n <= 1000)+ 8 * (n >1000 & n <= 1500)  + 9 * (n > 1500)#lb <- round(1.3* n^(1/5)) 
	 x <- seq(0, 1, length.out = 100)
f <-   dbeta(x, a, b)
mf <- matrix(f, nrow = nsim, ncol = length(x), byrow = T)
x <- seq(0, 1, length.out = 100)
  nh <- 1/(lb)
for(itr in 1 : nsim){
    theta <- rtheta[itr, ]
f <-   dbeta(x, a, b)
fest1[itr, ] <- (exp(eval.basis( x, objbasis) %*% theta) /sum(exp(basis %*% theta)  * weights[1, ]))


}
mse[sdis + 1, j] <- round(mean(apply(abs(fest1 - mf), 1, max)), 3)


	 pdf(paste('den_', n, sdis, '.pdf', sep=''))

	 x <- x
	 fest <- apply(fest1,2, median, na.rm = T)
	 qfest <- apply(fest1,2,quantile, c(0.05, 0.95),  na.rm = T)
	 plot(f~x, type = 'l', ylim = c(0, 4))
	 lines((fest ~x), lty = 2, col= 1)
	 if(sdis == FALSE){
	 	 lines(x, qfest[1, ], lty = 2)
	 	 lines(x, qfest[2, ], lty = 2)
	 }else{
		lines(x, qfest[1, ], lty = 2)
	 	lines(x, qfest[2, ], lty = 2)
		}
#legend('topright', paste('MSE = ', mse, sep = ''))
		dev.off()
		}
}




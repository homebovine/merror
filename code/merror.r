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
    y <- x + eps
    return(y)
}
denest <- function(theta, sdis, bs, y){
    bst<- matrix(exp(basis %*% theta), n, ng, byrow = T)
    integrand1 <- function(x, y, sdis){
    if(sdis == 0 ){
        deps <-  bst * dnorm(y - x,  0, 1) * weights
    }
    else if(sdix == 1){
        deps <- bst * ddoublex(y - x, 0, 1) * weights
    }else{
        deps <- bst * dunif(y - x, -1, 1) * weights
    }
}
    integrand2 <- function(x, sdis){
    if(sdis == 0 ){
        deps <- exp(predict(bs, x) %*% theta) 
    }
    else if(sdix == 1){
        deps <- exp(predict(bs, x) %*% theta) 
    }else{
        deps <- exp(predict(bs, x) %*% theta) 
    }
}
 #   integrand1 <- Vectorize(integrand1)
  #  integrand2 <- Vectorize(integrand2)
    intres1 <- intres2 <- rep(NA, n)
    #temp <- function(i){
    #     integrate(integrand1, 0, 1, y[i], sdis)$value
     #   
    #}
    temp <- function(i){
        integrand1(wn$nodes[i], y, sdis)* wn$weights[i]
    }
    intres1 <- apply(integrand1(nodes, y, sdis), 1, sum)# apply(sapply(1 : 10, temp), 1, sum)
    intres2 <- sum(bst[1, ] * weights[1, ])
    object <- -mean(log(intres1) ) +  log(intres2)
    
}

a <-  2
b <- 5
nsim <- 500
n <- 400
sdis <- 0
ng <-  50
knots = c(-0.5,   0,  0.5)#c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
lb = length(knots) + 3
theta0 <- rep(-0.5, lb)
rtheta <- matrix(NA, nsim, lb)
rvar <- vector('list')
wn <- gauss.quad(ng,kind="legendre",alpha=0,beta=0)
#wn$nodes <- seq(-0.998, 0.998, length.out = ng)
#wn$weights <- wn$weights * 1 /sqrt(1 - wn$nodes^2)#diff(c(-1, wn$nodes))
basis <- bs(wn$nodes, 3, knots, Boundary.knots = c(-1, 1))
nodes <- matrix(wn$nodes, n, ng, byrow = T)
weights <-  matrix(wn$weights, n, ng, byrow = T)
for(itr in 1 : nsim){
    print(itr)
    y <- simu(n, a, b, 0)
    my <- matrix(y, n, ng)
    
    #spg(theta0, denest, gr=NULL, method=3, lower=rep(-10, lb), upper=rep(10, lb), project=NULL, projectArgs=NULL, control=list(), quiet=FALSE, alertConvergence=TRUE, sdis, basis, my)

    res <- optim(theta0, denest, gr = NULL, sdis, basis, my, method = 'BFGS',  hessian = TRUE)
    rtheta[itr, ]<- res$par
    rvar[[itr]] <- -res$hessian/n
    
}


t <- rbeta(400, a, b)
o <- order(t)
t <- t[o]
x <- 2 * t - 1

f <-  1/2 * dbeta(t, a, b)
mtheta <- apply(rtheta, 2, mean)
fest1 <- matrix(NA, nsim, length(x))
for(itr in 1 : nsim){
    theta <- rtheta[itr, ]
fest1[itr, ] <- exp(predict(basis, x) %*% theta) /sum(exp(basis %*% theta)  * weights[1, ])
}
fest <- apply(fest1, 2, median)
fest  <- exp(predict(basis, x) %*% mtheta) /sum(exp(basis %*% mtheta)  * weights[1, ])

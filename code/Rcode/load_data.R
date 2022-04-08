# See doc: https://cran.r-project.org/web/packages/reticulate/index.html
library(reticulate)
np <- import("numpy")
# import data
# mat <- np$load("../matrix.npy") # for matrix in R^{T, S, R}
#setwd("~/Downloads/igor")
Y.kf <- np$load("../data/matrix.npy") # for full matrix in R^{T, S+R, S+R}
 

n <- dim(Y.kf)[1] # timesteps
# s, r hardcoded in python
s <- 18
r <- 31

p <- s+r
p.y = s*r
d <- 2
n.iter <- 10

alpha <- 1

# See in hazard function gam.sum.pred.vec
gam.sum.pred.mat = matrix( exp(alpha), n, p.y ) 
apply.zero.rate = TRUE # enable censoring of repeated events
x_0 = rnorm(p*d, 0, 1)
#then run kalman algorithm

######################################################## simulate data without censoring #########################
apply.zero.rate = FALSE
s <- 18
r <- 31
p <- s+r
p.y = s*r
d <- 2
n.iter <- 10
alpha <- 1
gam.sum.pred.mat = matrix( exp(alpha), n, p.y ) 

n=100
logis=function(time, a, b, c) a*(exp(10*(time/n-b))/(1 + exp(10*(time/n-b)))-.5) + c
a = runif(p*d, -1, 1)*1
b = rep(runif(p, 0.2, 0.8), d)
c = runif(p*d,-2,2) 
X=matrix(NA, n, p*d)
for(i in 1:(p*d))X[, i] = logis(1:n, a[i], b[i], c[i]) 
matplot(X, ty="l")
# get simulated data instead
Y.kf = matrix(NA, n, p.y)
lambda.true = matrix(NA, n, p.y)
for(t in 1:n){
  lambda.true[t, ] = h.function(X[t, ], gam.sum.pred.mat[t, ])
  Y.kf[t, ] = rpois(p.y, lambda = lambda.true[t, ])
}
x_0 = X[1, ]#rnorm(p*d, 0, 1)
#then run kalman algorithm

############## simulate from a censored process
alpha=0.5
gam.sum.pred.mat = matrix( exp(alpha), n, p.y ) 

apply.zero.rate = TRUE
Y.kf = matrix(NA, n, p.y)
lambda.true = matrix(NA, n, p.y)
for(t in 1:n){
  lambda.true[t, ] = h.function(X[t, ], gam.sum.pred.mat[t, ])
  if(t>1 & apply.zero.rate) lambda.true[t, ] = lambda.true[t, ]*zero.rate(matrix(Y.kf[1:(t-1), ], ncol=p.y))
  Y.kf[t, ] = rpois(p.y, lambda = lambda.true[t, ])
}
x_0 = X[1, ]#rnorm(p*d, 0, 1)
#then run kalman algorithm



# See doc: https://cran.r-project.org/web/packages/reticulate/index.html
library(reticulate)
np <- import("numpy")
# import data
# mat <- np$load("../matrix_full.npy") # for matrix in R^{T, S, R}
mat <- np$load("../matrix_full.npy") # for full matrix in R^{T, S+R, S+R}
Y.kf <- mat

n <- dim(Y_kf)[1] # timesteps
# s, r hardcoded in python
s <- 18
r <- 31
# My matrix size is p*p and not p*(p-1)/2 though. What am I doing wrong? 
p <- s+r
p.y = p*(p-1)/2
d <- 2
d.true <- d
n.iter <- 1

sender.true = rep(0, s)
receiver.true = rep(0, r)
alpha <- 3
spline.time.true = rep(alpha, n)

# See in hazard function gam.sum.pred.vec
gam.sum.pred.mat = matrix( 2*exp(alpha), n, p*(p-1)/2)  # its twice the rate!!!!!!!

X <- NaN


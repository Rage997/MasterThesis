
#uncomment following 15 lines
# require(Matrix)
# require(mvtnorm)

negBin=F
alpha = 3
# n.iter = 10
p = 10
p.y = p*(p-1)/2
n = 100
d = 2
d.true = d
dd=1
n.clust = 3



time = 1:n
# sender.true = rnorm(p, 0, 0.5)
# receiver.true = rnorm(p, 0, 0.5)
# spline.time.true = 3 + c(-1, 1)[rbinom(1, 1, 0.5) + 1] * (exp(10*time/n - 5)/(1 + exp(10*time/n - 5))-.5)

sender.true = rep(0, p)
receiver.true = rep(0, p)
spline.time.true = rep(alpha, n)



logis=function(time, a, b, c) a*(exp(10*(time/n-b))/(1 + exp(10*(time/n-b)))-.5) + c

# time=seq(0, 10, length.out = n)
a = runif(p*d.true, -1, 1)*1
b = rep(runif(p, 0.2, 0.8), d.true)
c.clust=0
if(n.clust>1){
  print("clustering implemented for dimension 2 only")
  c.clust = 5*matrix(c(0,0,1,1,-1,1,-1,-1,1,-1), ncol=2, byrow=T)[rep(1:n.clust, each=floor(p/n.clust) , length.out=p), ]
  c.clust = as.vector(t(c.clust))
}
c = runif(p*d.true,-dd,dd) + c.clust
# cc.true = cc = 3
X=matrix(NA, n, p*d.true)

for(i in 1:(p*d.true))X[, i] = logis(time, a[i], b[i], c[i]) 

# for(i in 1:n)plot(matrix(X[i,], ncol=2, byrow=T ), main=i)
# Plot trajectories 
matplot(X, ty="l")


lambda.true = array(NA, dim = c(p, p, n))
Y.array = array(NA, dim = c(p, p, n))
for(t in 1:n)for(i in 1:p)for(j in 1:p){
  if( i!=j ){
    lambda.true[i,j,t] = exp(  spline.time.true[t] + sender.true[i] + receiver.true[j] - as.numeric(dist(matrix(X[ t, c((i-1)*d.true + 1:d.true, (j-1)*d.true + 1:d.true) ], ncol = d.true, byrow=T))^2) ) 
    Y.array[ i, j, t] = if(negBin==T){rnbinom( 1, 1, mu=lambda.true[ i, j, t])}else{rpois( 1, lambda.true[ i, j, t])}
  }
}

plot(time, lambda.true[ 1, 2, ], ylim=c(0, max(Y.array, na.rm = T)), ty="l")
# plot( time, Y.array[ 1, 2, ] )
for(i in 1:p)for(j in 1:p){
  if( i!=j ){
    lines(time, lambda.true[ i, j, ], col=i+j)
    # points( time, Y.array[ i, j, ], ylim=c(0, max(Y, na.rm = T)) )
  }
}


Y.kf = matrix( NA, n, p*(p-1)/2 )
lambda.true.sum = matrix(NA, n, p*(p-1)/2) 
for(t in 1:n){
  ii = 1
  for(i in 1:(p-1))for(j in (i+1):p){
    lambda.true.sum[t, ii] = lambda.true[ i, j, t] + lambda.true[ j, i, t]
    Y.kf[t, ii] = Y.array[ i, j, t] + Y.array[ j, i, t]
    ii = ii + 1
  }
}


x_0  = as.vector(matrix( X[ 1, ], nrow=d.true)[1:d, ])
gam.sum.pred.mat = matrix( 2*exp(alpha), n, p*(p-1)/2)  # its twice the rate!!!!!!!

Y.kf.static = apply(Y.kf, 2, sum)
Y.kf.static = matrix(Y.kf.static, nrow=1)
gam.sum.pred.mat.static = matrix( n*2*exp(alpha), 1, p*(p-1)/2)  # its twice the rate!!!!!!!

##########################################################################################################
# n.time.points = n/4
# D.hat = matrix(0, p, p)
# D.hat[lower.tri(D.hat)] = sqrt( pmax(0, log(2) + 3 - log(apply(Y.kf[1:(n.time.points),] + 1e-30, 2, mean) )) )
# D.hat = D.hat + t(D.hat)
# tmp = cmdscale( D.hat, k=d, add=F )
# # par(mfrow=c(1,2))
# # plot(tmp, col=1:p)
# # plot(matrix( X[ 1, ], ncol=d.true, byrow=T), col=1:p)
# # plot(matrix( X[ 1, ], ncol=d.true, byrow=T)[ , 1],  tmp[, 1])
# # plot(matrix( X[ 1, ], ncol=d.true, byrow=T)[ , 1],  tmp[, 1])
# 
# 
# # plot(dist(matrix( X[ 1, ], ncol=d.true, byrow=T)), dist(tmp) )
# # abline(0,1, col=2)
# 
# x_0  = as.vector(t(tmp))
# # x_0 = runif(p*d, -20, 20)
# 
# 
# llik = function(beta, x){
#   lambda =  exp( beta - as.vector(dist(matrix(x, ncol=d, byrow=T))^2))
#   sum(-lambda * n.time.points) + sum(Y.kf[1:(n.time.points),] %*% log(lambda))
# }
# 
# B=100000
# llik.save = -Inf
# x.mcmc = matrix(NA, B, p*d)
# beta.mcmc = rep(NA, B)
# 
# eps1 = 0.01
# eps2 = 0.01
# x.mcmc[1, ] = x_0
# beta.mcmc[1] = log(2) + 3
# 
# 
# #run simulated annealing mcmc
#   
#   for(b in 2:B){
#     x.mcmc[b, ] = x.mcmc[b-1, ]
#     beta.mcmc[b] = beta.mcmc[b-1]
#     x.mcmc.proposed = x.mcmc[b-1, ] + rnorm( p*d, 0, eps1)
#     beta.mcmc.proposed = beta.mcmc[b-1] + rnorm( 1, 0, eps2)
#     
#     llik.proposed = llik(beta.mcmc.proposed, x.mcmc.proposed)
#     if(runif(1) < 1000/b *exp(llik.proposed - llik.save ) ){
#       llik.save = llik.proposed
#       x.mcmc[b, ] = x.mcmc.proposed
#       beta.mcmc[b] = beta.mcmc.proposed
#     }
#   }
# 
# matplot(x.mcmc, ty="l")
# matplot(beta.mcmc, ty="l")
# 
# x_0 = x.mcmc[B, ]
# 
# 
# 







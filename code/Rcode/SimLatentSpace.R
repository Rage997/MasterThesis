# Params
negBin=F
alpha = 3
# n.iter = 10
p = 10
p.y = p*(p-1)/2
n = 100
# TODO d, and d.true are the same?
d = 2
d.true = d
dd=1
n.clust = 3

time = 1:n

# Init vector with 0
sender.true = rep(0, p)
receiver.true = rep(0, p)
spline.time.true = rep(alpha, n)

# TODO: what is this?
logis=function(time, a, b, c) a*(exp(10*(time/n-b))/(1 + exp(10*(time/n-b)))-.5) + c

# Generate a, b, c from uniform distribution using diffrent parameters
# In particular, c will separate the trajectories into clusters.

# Uniform distribution ruinif(x, min, max)
a = runif(p*d.true, -1, 1)*1
b = rep(runif(p, 0.2, 0.8), d.true)
# Clustering (for what means?)
c.clust=0
if(n.clust>1){
  print("clustering implemented for dimension 2 only")
  c.clust = 5*matrix(c(0,0,1,1,-1,1,-1,-1,1,-1), ncol=2, byrow=T)[rep(1:n.clust, each=floor(p/n.clust) , length.out=p), ]
  c.clust = as.vector(t(c.clust))
}
c = runif(p*d.true,-dd,dd) + c.clust


X=matrix(NA, n, p*d.true)

for(i in 1:(p*d.true))X[, i] = logis(time, a[i], b[i], c[i]) 
# Plot trajectories 
matplot(X, ty="l")


# TODO: no idea? estimating lambda?
lambda.true = array(NA, dim = c(p, p, n))
Y.array = array(NA, dim = c(p, p, n))
for(t in 1:n)for(i in 1:p)for(j in 1:p){
  if( i!=j ){
    lambda.true[i,j,t] = exp(  spline.time.true[t] + sender.true[i] + receiver.true[j] - as.numeric(dist(matrix(X[ t, c((i-1)*d.true + 1:d.true, (j-1)*d.true + 1:d.true) ], ncol = d.true, byrow=T))^2) ) 
    Y.array[ i, j, t] = if(negBin==T){rnbinom( 1, 1, mu=lambda.true[ i, j, t])}else{rpois( 1, lambda.true[ i, j, t])}
  }
}

# plot(time, lambda.true[ 1, 2, ], ylim=c(0, max(Y.array, na.rm = T)), ty="l")

for(i in 1:p)for(j in 1:p){
  if( i!=j ){
    lines(time, lambda.true[ i, j, ], col=i+j)
  }
}

# TODO: go over this again?
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

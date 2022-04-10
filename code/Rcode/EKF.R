# require(Matrix)
# require(mgcv)
require(raster)

# Taken from SimLatentSpace
# n.iter = 500
# d = 3
# p = 5
# p.y = p*(p-1)/2
# n = 100

# transition matrix : KF post ~ T*prior 
# transition noise: Q
# observation matrix: taylor expansion hazard funciton
# observation noise: var(poisson)
# initial state prior: X0

# hazard function / rate -> poisson variable
h.function=function(x.vec, gam.sum.pred.vec){
  gam.sum.pred.vec * exp( -as.vector(t(pointDistance(matrix(x.vec[1:(s*d)], ncol=d, byrow=T), matrix(x.vec[(s*d+1):((s+r)*d)], ncol=d, byrow=T), lonlat=F, allpairs = T )^2)))
}

d.h.function = function(x.i, x.j, gam.sum.pred.ij){
  # First derivative (taylor's expansion) -> latent space
  dist.tmp = dist(rbind(x.i, x.j))^2
  2 * gam.sum.pred.ij * exp(- dist.tmp) * c(x.j - x.i, x.i - x.j)
}

# Inizialiation T=0
X.prior = X.post = matrix(NA, n, p*d) # points in latent space (Expected value of X)
P.prior = P.post = array(NA, dim=c(p*d,p*d,n)) # variance of X

#tuning parameters (transition noise)
Q = diag(0.001, p*d) # X_t+1  = X_T + eps_t : sigma is var(eps)

# Inizialization: starting points    TRY WITH DIFFERENT X STARTING POINTS!
P_0.smooth = P_0 = Q
#x_0.smooth = x_0 = runif(p*d, -1, 1)
# TODO generate it from a gaussian or whatever
x_0.smooth = x_0 #= X[1, ]#rnorm(p*d, 0, 1)
  #as.vector(matrix( X[ 1, ], nrow=d.true)[1:d, ])

# Populate H_ij = [derivative function(i,j)] -> sparse matrix
build.H = function(x.prior, gam.sum.pred.vec){
  H = matrix(0, p.y, p*d)
  ii = 1
  for(i in 1:s)for(j in (s+1):p){
    H[ii, c(d*i - (d-1):0, d*j - (d-1):0)] = d.h.function(x.prior[c(d*i - (d-1):0)], x.prior[c(d*j - (d-1):0)] , gam.sum.pred.vec[ii] )  
    ii = ii + 1 
  }
  return(H)
}

zero.rate = function(y)apply(y, 2, function(x)!any(x>0))
eps = 1e-5

for(iter in 1:n.iter){
  for(t in 1:n){ #Kalmann Filter
    if(t==1){ # init (page 9 igor paper)
      X.prior[ t, ] = x_0
      P.prior[ , , t] = P_0 + Q
    }else{ # a-prior distribution
      X.prior[ t, ] = X.post[ t-1, ]
      P.prior[ , , t] = P.post[ , , t-1] + Q
    }
    # See page 9
    H = build.H(X.prior[ t, ], gam.sum.pred.mat[ t, ]) #first derivative (gam.sum is an offset) on X_prior
    # you can do taylor: mu(x) = mu(x0) + d_h(x - x0)
    mu = h.function(X.prior[ t, ], gam.sum.pred.mat[ t, ]) # poisson parameter -> variance
    if(t>1 & apply.zero.rate){
      # Check if at least one interaction happened in the past [1:t-1] between S and R
      mu = mu*zero.rate(matrix(Y.kf[1:(t-1), ], ncol=p.y)) # -> this makes the linear predictor do not consider re-invasion 
      H = H*zero.rate(matrix(Y.kf[1:(t-1), ], ncol=p.y))# see: we set the gradient to zero if an interaction happened in the past
    }
    vv = mu 
    R = diag(vv)
    R.inv = diag(1/vv)
    R.inv[R.inv==Inf] = 0
    # Sherman-Morrison-Woodbury identity (page 31) -> kinda a lie but okay
    K = tcrossprod(solve(solve(P.prior[ , , t], tol = 1e-70) + crossprod(H, crossprod(R.inv, H)), tol = 1e-70),  t( crossprod(H, R.inv)) )
    
    X.post[ t, ] = X.prior[ t, ] + crossprod(t(K), ( Y.kf[ t, ] - mu )) # see algo 1 pg 10: expected value
  
    P.post[ , , t] = crossprod(t(diag(p*d) - crossprod(t(K), H) ), P.prior[ , , t]) # see algo 1 pg 10: variance
    P.post[ , , t] = (P.post[ , , t] + t(P.post[ , , t]))/2
    if(any( eigen(P.post[ , , t])$values<=0 ) ){
      print(paste("P.post not positive defined at iteration ", iter, ", time, ", t, sep = ""))
      P.post[ , , t] = as.matrix(nearPD(P.post[ , , t])$mat)  
    }  
  }
  
  # Smoother step algo 2 pg 11
  X.smooth = X.post
  P.smooth = P.post
  B.save = array(NA, c(p*d, p*d, n) )
  
  for(t in n:1){
    #B = crossprod(t(P.post[ , , t]), solve(P.prior[ , , t]))
    if(t==1){
      B = crossprod(t(P_0), solve(P.prior[ , , t], tol=1e-30))
      x_0.smooth = x_0 + crossprod(t(B), X.smooth[ t, ] - X.prior[ t, ])
      P_0.smooth = P_0 + tcrossprod(crossprod(t(B), P.smooth[ , , t] - P.prior[ , , t]), B)
      P_0.smooth = (P_0.smooth + t(P_0.smooth))/2
      if(any( eigen(P_0.smooth)$values<=0 ) ){
        print(paste("P.smooth not positive defined at iteration ", iter, ", time, ", t, sep = ""))
        P_0.smooth = as.matrix(nearPD(P_0.smooth)$mat)  
      }  
    }else{
      B = crossprod(t(P.post[ , , t-1]), solve(P.prior[ , , t], tol=1e-30))
      X.smooth[ t-1, ] = X.post[ t-1, ] + crossprod(t(B), X.smooth[ t, ] - X.prior[ t, ])
      P.smooth[ , , t-1] = P.post[ , , t-1] + tcrossprod(crossprod(t(B), P.smooth[ , , t] - P.prior[ , , t]), B)
      P.smooth[ , , t-1] = (P.smooth[ , , t-1] + t(P.smooth[ , , t-1]))/2
      if(any( eigen(P.smooth[ , , t-1])$values<=0 ) ){
        print(paste("P.smooth not positive defined at iteration ", iter, ", time, ", t, sep = ""))
        P.post[ , , t] = as.matrix(nearPD(P.post[ , , t])$mat)  
      }  
    }
    B.save[ , , t] = B 
  }
  
  x_0 = x_0.smooth
  P_0 = P_0.smooth
  
  # maximization step pg 12 see star (Q is sigma)
  Q.hat = NULL
  for(t in 1:n){
    if(t==1){
      Q.hat = tcrossprod(X.smooth[t, ] - x_0.smooth, X.smooth[t, ] - x_0.smooth) + P.smooth[ , , t] + P_0.smooth - crossprod(t(B.save[ , , t]), P.smooth[ , , t]) - t(crossprod(t(B.save[ , , t]), P.smooth[ , , t])) 
    }else{
      Q.hat = Q.hat + tcrossprod(X.smooth[t, ] - X.smooth[t-1, ], X.smooth[t, ] - X.smooth[t-1, ]) + P.smooth[ , , t] + P.smooth[ , , t-1] - crossprod(t(B.save[ , , t]), P.smooth[ , , t]) - t(crossprod(t(B.save[ , , t]), P.smooth[ , , t])) 
    }
  }
  
  Q.hat = Q.hat/n
  Q = diag(rep(apply(matrix(diag(Q.hat), ncol=p), 2, mean), each=d))

  #fit additional parameters in linear predictor
  D=Y.kf
  for(t in 1:n){
    x.vec = X.smooth[t, ]
    D[t, ] = -as.vector(t(pointDistance(matrix(x.vec[1:(s*d)], ncol=d, byrow=T), matrix(x.vec[(s*d+1):((s+r)*d)], ncol=d, byrow=T), lonlat=F, allpairs = T )^2))
    if(t>1 & apply.zero.rate){
      D[t, !zero.rate(matrix(Y.kf[1:(t-1), ], ncol=p.y))] = NA
    }
  }

  alpha.est = glm(as.vector(Y.kf[!is.na(D)]) ~ 1 + offset(as.vector(D[!is.na(D)])), family=poisson(link = "log"))$coef
  gam.sum.pred.mat = matrix( exp(alpha.est), n, p.y )
}


#empirical Kullback Leibler
# KL=NULL
# for(k in 1:100){
#   Y.gen = matrix(rpois( length(Y.kf), lambda.true.sum ), nrow(Y.kf), ncol(Y.kf))
#   KL= c(KL, mean( dpois(Y.gen, lambda.true.sum, log=T) - dpois(Y.gen, E1*gam.sum.pred.mat, log=T)  ) )
# }
# KL=mean(KL)

#re-expand data
# if(scale.data==T){
# gam.sum.pred.mat = gam.sum.pred.mat %*% diag(1/res.fact) #rescale columns
# Y.kf = Y.kf %*% diag(1/res.fact) #rescale columns
# }

# ref = fit.gam$coefficients[(length(fit.gam$coefficients)-2*p+1):length(fit.gam$coefficients)] - c(sender.true - mean(sender.true), receiver.true - mean(receiver.true)) 
# ref = sqrt(sum(rss^2)/(2*p))


#save(KL, compTime, file=paste("/home/artico/simulation1/EKF_simId", bbb,"_n", n,"_p", p,"_d", d,"_dTrue", d.true,"_dd", dd,".RData", sep="") )

#res = c( KL, as.numeric(compTime[3]) )


# PLOTS
lambda.est = matrix(NA, n, p.y)
for(t in 1:n){
  lambda.est[ t, ] = h.function(X.smooth[t, ], gam.sum.pred.mat[t, ])
    if(t>1 & apply.zero.rate)lambda.est[ t, ] = lambda.est[ t, ]*zero.rate(matrix(Y.kf[1:(t-1), ], ncol=p.y))
}

i.samp  = (1:p.y)[apply(Y.kf,2, function(x)any(x>0) )]
i.samp  = sample(p.y, 50)
for(i in i.samp){
  plot(Y.kf[, i], ty="p", main=i)
  #lines(lambda.true[, i])
  lines(lambda.est[, i], col=2)
}

#get estimate for alpha and Q
diag(Q)
alpha.est

matplot(X.smooth, ty="l")e
for(i in 1:n)plot(matrix(X.smooth[i,], ncol=2, byrow=T ), main=i, ylim=c(-2, 4), xlim=c(-2, 3))

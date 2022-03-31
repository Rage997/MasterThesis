res=NULL
require(Matrix)
require(mgcv)

# path = "~/Desktop/state space/Patent Data/results/"

# n.iter = 500
# d = 3
# p = 5
# p.y = p*(p-1)/2
# n = 100
# B.mc = 500
# d.true = d
# dd=1




# network.type = "directed" # "undirected"
# 
# 
# Y.gam = matrix( NA, n, p*(p-1) ) 
# sender.gam = matrix( NA, n, p*(p-1) )
# receiver.gam = matrix( NA, n, p*(p-1) )
# time.gam = matrix( NA, n, p*(p-1) )
# 
# Y.kf = matrix( NA, n, p*(p-1)/2 )
# 
# for(t in 1:n){
#   ii = 1  
#   for(i in 1:(p-1))for(j in (i+1):p){
#     Y.kf[t, ii] = Y.array[ i, j, t] + Y.array[ j, i, t]
#     
#     Y.gam[t, ii] = Y.array[ i, j, t]
#     sender.gam[t, ii] = i
#     receiver.gam[t, ii] = j
#     time.gam[t, ii] = t
#     #here you can add additional covariates
#     ii = ii + 1
#   }
#   for(j in 1:(p-1))for(i in (j+1):p){
#     Y.gam[t, ii] = Y.array[ i, j, t]
#     sender.gam[t, ii] = i
#     receiver.gam[t, ii] = j
#     time.gam[t, ii] = t
#     #here you can add additional covariates
#     ii = ii + 1
#   }
# }
# 
# data.gam = as.data.frame(cbind(as.vector(t(Y.gam)), as.vector(t(sender.gam)), as.vector(t(receiver.gam)), as.vector(t(time.gam)), NA ) )
# names(data.gam) = c("Y", "sender", "receiver", "time", "E")
# data.gam$sender = as.factor(data.gam$sender)
# data.gam$receiver = as.factor(data.gam$receiver)
# 
# 
# 
# 
# 
# 
# formula.lin.pred = if(network.type=="directed"){
#   formula(Y ~ offset(log(E)) )#+ s(time) + s(sender, bs="re") + s(receiver, bs="re")) 
# }else{
#   formula(Y ~ offset(log(2) + log(E)) ) #+ s(time) + s(sender, bs="re") + s(receiver, bs="re")) 
# }


# hazard function / rate -> poisson variable
h.function=function(x.vec, gam.sum.pred.vec) gam.sum.pred.vec * exp( -as.vector(dist(matrix(x.vec, ncol=d, byrow=T))^2))
d.h.function = function(x.i, x.j, gam.sum.pred.ij){
  # First derivative (taylor's expansion) -> latent space
  dist.tmp = dist(rbind(x.i, x.j))^2
  2 * gam.sum.pred.ij * exp(- dist.tmp) * c(x.j - x.i, x.i - x.j)
}

# Inizialiation T=0
# H=diag(p)
X.prior = X.post = matrix(NA, n, p*d) # points in latent space (Expected value of X)
P.prior = P.post = array(NA, dim=c(p*d,p*d,n)) # variance of X

#tuning parameters
# cc = 12
Q = diag(0.001, p*d) # X_t+1  = X_T + eps_t : sigma is var(eps)
# R = diag(0.01^2, p.y)


# Inizialization: starting points    TRY WITH DIFFERENT X STARTING POINTS!
P_0.smooth = P_0 = Q
#x_0.smooth = x_0 = runif(p*d, -1, 1)
x_0.smooth = x_0 = as.vector(matrix( X[ 1, ], nrow=d.true)[1:d, ])

# Populate H_ij = [derivative function(i,j)] -> sparse matrix
build.H = function(x.prior, gam.sum.pred.vec){
  H = matrix(0, p.y, p*d)
  ii = 1
  for(i in 1:(p-1))for(j in (i+1):p){
    H[ii, c(d*i - (d-1):0, d*j - (d-1):0)] = d.h.function(x.prior[c(d*i - (d-1):0)], x.prior[c(d*j - (d-1):0)] , gam.sum.pred.vec[ii] )  
    ii = ii + 1 
  }
  return(H)
}




hess = function(x.i, x.j){
  dist.tmp=dist(rbind(x.i, x.j))^2
  res = matrix(NA, 2*d, 2*d)
  
  diag(res)[1:d] = diag(res)[(d+1):(2*d)] = ( 4*((x.i-x.j)^2) - 2 ) * exp(-dist.tmp)
  if(d==1){
    res[1, 2] = res[2, 1] = - res[1, 1]
  }
  if(d==2){
    res[1, 3] = res[3, 1] = - res[1, 1]
    res[2, 4] = res[4, 2] = - res[2, 2]
    
    res[1, 4] = res[4, 1] = res[2, 3] = res[3, 2] = -4*(x.i-x.j)[1]*(x.i-x.j)[2] * exp(-dist.tmp)
    res[1, 2] = res[2, 1] = res[3, 4] = res[4, 3] = - res[1, 4]
  }
  if(d==3){
    res[1, 4] = res[4, 1] = - res[1, 1]
    res[2, 5] = res[5, 2] = - res[2, 2]
    res[3, 6] = res[6, 3] = - res[3, 3]
    
    res[1, 5] = res[5, 1] = res[2, 4] = res[4, 2] = -4*(x.i-x.j)[1]*(x.i-x.j)[2] * exp(-dist.tmp)
    res[1, 6] = res[6, 1] = res[3, 4] = res[4, 3] = -4*(x.i-x.j)[1]*(x.i-x.j)[3] * exp(-dist.tmp)
    res[2, 6] = res[6, 2] = res[3, 5] = res[5, 3] = -4*(x.i-x.j)[2]*(x.i-x.j)[3] * exp(-dist.tmp)
    res[1, 2] = res[2, 1] = res[4, 5] = res[5, 4] = - res[1, 5]
    res[1, 3] = res[3, 1] = res[4, 6] = res[6, 4] = - res[1, 6]
    res[2, 3] = res[3, 2] = res[5, 6] = res[6, 5] = - res[2, 6]
  }
  res
}


# hess2 = function(x1, x2){
#   if( d>3 | d<1 )stop("the hessian is implemented for dimension 2 and 3 only")
#   dist.tmp=dist(rbind(x1, x2))
#   res = matrix(NA, 2*d, 2*d)
#   
#   diag(res)[1:d] = diag(res)[(d+1):(2*d)] = ( (x1-x2)^2  - dist.tmp^2) /(dist.tmp^3)
#   
#   if(d==2){
#     res[1, 3] = res[3, 1] = - res[1, 1]
#     res[2, 4] = res[4, 2] = - res[2, 2]
#     
#     res[1, 4] = res[4, 1] = res[2, 3] = res[3, 2] = -(x1-x2)[1]*(x1-x2)[2]/(dist.tmp^3)
#     res[1, 2] = res[2, 1] = res[3, 4] = res[4, 3] = - res[1, 4]
#   }else if(d==3){
#     res[1, 4] = res[4, 1] = - res[1, 1]
#     res[2, 5] = res[5, 2] = - res[2, 2]
#     res[3, 6] = res[6, 3] = - res[3, 3]
#     
#     res[1, 5] = res[5, 1] = res[2, 4] = res[4, 2] = -(x1-x2)[1]*(x1-x2)[2]/(dist.tmp^3)
#     res[1, 6] = res[6, 1] = res[3, 4] = res[4, 3] = -(x1-x2)[1]*(x1-x2)[3]/(dist.tmp^3)
#     res[2, 6] = res[6, 2] = res[3, 5] = res[5, 3] = -(x1-x2)[2]*(x1-x2)[3]/(dist.tmp^3)
#     res[1, 2] = res[2, 1] = res[4, 5] = res[5, 4] = - res[1, 5]
#     res[1, 3] = res[3, 1] = res[4, 6] = res[6, 4] = - res[1, 6]
#     res[2, 3] = res[3, 2] = res[5, 6] = res[6, 5] = - res[2, 6]
#   }
#   res
# }

# scale.data = T
# 
# if(scale.data==T){res.fact = 1/pmax(apply(Y.kf, 2, max), 1)}else{res.fact = rep(1, ncol(Y.kf))} #rescale all varaibles in [0,1]
# gam.sum.pred.mat = gam.sum.pred.mat %*% diag(res.fact) #rescale columns
# Y.kf = Y.kf %*% diag(res.fact) #rescale columns

# Y.kf = matrix(pmax( as.vector(Y.kf), 1), nrow=n)

#x_0 = 0.5*x_0

repeat.filter=1

hist.list = NULL
gam.save = NULL

#system.time(

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
    mu = h.function(X.prior[ t, ], gam.sum.pred.mat[ t, ]) # poisson parameter -> variance
    vv = mu 
    for(k in repeat.filter){ # not necessary
      # vv = res.fact * vv
      R = diag(vv)
      R.inv = diag(1/vv)
      # K = P.prior[ , , t] %*% t(H) %*% solve( H %*% P.prior[ , , t] %*% t(H) + R )
      # X.post[ t, ] = X.prior[ t, ] + K %*% ( Y[ t, ] - h.function(X.prior[ t, ]) )
      # P.post[ , , t] = (diag(p*d) - K %*% H) %*% P.prior[ , , t]
      # tryCatch(
      #   K <<- crossprod( t(tcrossprod( P.prior[ , , t], H)), solve( tcrossprod(crossprod(t(H), P.prior[ , , t]), H) + R, tol=1e-70) ),
      #   error = function(e) K <<- crossprod( t(tcrossprod( P.prior[ , , t], H)), solve( off + tcrossprod(crossprod(t(H), P.prior[ , , t]), H) + R, tol=1e-70) )
      # )
      # Sherman-Morrison-Woodbury identity (page 31) -> kinda a lie but okay
      K = tcrossprod(solve(solve(P.prior[ , , t], tol = 1e-70) + crossprod(H, crossprod(R.inv, H)), tol = 1e-70),  t( crossprod(H, R.inv)) )
      
      # K = crossprod( t(tcrossprod( P.prior[ , , t], H)), solve( tcrossprod(crossprod(t(H), P.prior[ , , t]), H) + R, tol=1e-70) )
      X.post[ t, ] = X.prior[ t, ] + crossprod(t(K), ( Y.kf[ t, ] - mu )) # see algo 1 pg 10: expected value
      vv = h.function(X.post[ t, ], gam.sum.pred.mat[ t, ]) # useless
    }
    
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

  
  # E1 = matrix( NA, n, p*(p-1)/2 ) 
  # # E2 = matrix( NA, n, p*(p-1)/2 ) 
  # for(t in 1:n){
  #   ii=1
  #   for(i in 1:(p-1)){
  #     for(j in (i+1):p){
  #       i.ind = (i-1)*d + (1:d)
  #       j.ind = (j-1)*d + (1:d)
  #       x.i = X.smooth[t, i.ind]
  #       x.j = X.smooth[t, j.ind]
  #       E1[t, ii] = exp(-dist(rbind(x.i, x.j))^2) + 0.5 * sum(diag(hess(x.i, x.j) %*% P.smooth[ c(i.ind, j.ind), c(i.ind, j.ind), t]) )
  #       # E2[ii, t] = -dist(rbind(x.i, x.j)) + 0.5 * sum(diag(hess2(x.i, x.j) %*% P.smooth[ c(i.ind, j.ind), c(i.ind, j.ind), t]) )
  #       if(E1[t, ii]<0 | E1[t, ii]>1){
  #         print(paste("Nodes too close at iteration ", iter, ", time, ", t,", at pair (", i,", ", j, ")", sep = ""))
  #         E1[t, ii] = exp(-dist(rbind(x.i, x.j))^2)
  #         print(E1[t, ii])
  #         if(E1[t, ii] == 0) print(paste("Error: infinite distance at iteration ", iter, ", time, ", t,", at pair (", i,", ", j, ")", sep = ""))
  #       }
  #       ii = ii + 1
  #     }
  #   }
  # }
  # # #cc = log(sum(Y)) - log(sum(E1))
  # # # E.llik[iter] = -n*log(det.Q)/2 - sum( E1 * exp(cc) ) + sum( as.vector(t(Y)) * (as.vector(t(E2)) + cc)  )
  # # #CHECK SE HESSIANO E QUESTA FORMULA SONO CORRETTI
  #
  #
  # data.gam$E = as.vector(t(cbind(E1, E1)))
  # SEE: https://m-clark.github.io/
  # fit.gam = gam( formula.lin.pred, family = poisson(link = "log"), data = data.gam )
  # gam.sum.pred.mat = matrix( predict( fit.gam, type="response" ), nrow = n, ncol = p*(p-1), byrow=T )
  # 
  # gam.sum.pred.mat = if(network.type=="directed"){ gam.sum.pred.mat/cbind(E1, E1) }else{ gam.sum.pred.mat/(2*cbind(E1, E1)) }
  # gam.sum.pred.mat = gam.sum.pred.mat[ , 1:(p*(p-1)/2)] + gam.sum.pred.mat[ , (p*(p-1)/2 + 1):(p*(p-1))]
  # 
  # gam.save = c(gam.save, list(fit.gam$coefficients))
  

  # cat(floor(100*iter/n.iter), "%" , " \r")
  # flush.console()
  
}

#) -> compTime #system.time end

#empirical Kullback Leibler

E1 = matrix( NA, n, p*(p-1)/2 )
# E2 = matrix( NA, n, p*(p-1)/2 )
for(t in 1:n){
  ii=1
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      i.ind = (i-1)*d + (1:d)
      j.ind = (j-1)*d + (1:d)
      x.i = X.smooth[t, i.ind]
      x.j = X.smooth[t, j.ind]
      E1[t, ii] = exp(-dist(rbind(x.i, x.j))^2) + 0.5 * sum(diag(hess(x.i, x.j) %*% P.smooth[ c(i.ind, j.ind), c(i.ind, j.ind), t]) )
      # E2[ii, t] = -dist(rbind(x.i, x.j)) + 0.5 * sum(diag(hess2(x.i, x.j) %*% P.smooth[ c(i.ind, j.ind), c(i.ind, j.ind), t]) )
      if(E1[t, ii]<0 | E1[t, ii]>1){
        print(paste("Nodes too close at iteration ", iter, ", time, ", t,", at pair (", i,", ", j, ")", sep = ""))
        E1[t, ii] = exp(-dist(rbind(x.i, x.j))^2)
        print(E1[t, ii])
        if(E1[t, ii] == 0) print(paste("Error: infinite distance at iteration ", iter, ", time, ", t,", at pair (", i,", ", j, ")", sep = ""))
      }
      ii = ii + 1
    }
  }
}

#re-expand data
# if(scale.data==T){
# gam.sum.pred.mat = gam.sum.pred.mat %*% diag(1/res.fact) #rescale columns
# Y.kf = Y.kf %*% diag(1/res.fact) #rescale columns
# }

KL=NULL
for(k in 1:100){
  Y.gen = matrix(rpois( length(Y.kf), lambda.true.sum ), nrow(Y.kf), ncol(Y.kf))
  KL= c(KL, mean( dpois(Y.gen, lambda.true.sum, log=T) - dpois(Y.gen, E1*gam.sum.pred.mat, log=T)  ) )
}
KL=mean(KL)

# ref = fit.gam$coefficients[(length(fit.gam$coefficients)-2*p+1):length(fit.gam$coefficients)] - c(sender.true - mean(sender.true), receiver.true - mean(receiver.true)) 
# ref = sqrt(sum(rss^2)/(2*p))


#save(KL, compTime, file=paste("/home/artico/simulation1/EKF_simId", bbb,"_n", n,"_p", p,"_d", d,"_dTrue", d.true,"_dd", dd,".RData", sep="") )

res = c( KL, as.numeric(compTime[3]) )


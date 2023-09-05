## install the following packages, respectively for fast SVD/EVD and KNN
## gtools is for finding the column permutation
## do not need to library them, already using namespace to call their respective functions
## define packages to install
packages <- c("RSpectra","FNN", "gtools")
## install all packages that are not already installed
install.packages(setdiff(packages, rownames(installed.packages())))
library(Matrix) 
# I suppose this package comes with base R? otherwise, install it as well.

## projection onto the simplex, a is the "radius" of the simplex
.proj <- function(y,a=1){
  v = y
  rho = (sum(y)-a)/length(y)
  vt.len = length(v)
  repeat{
    v = v[which(v>rho)]
    v.len = length(v)
    rho = (sum(v)-a)/v.len
    if (v.len == vt.len){
      break
    } else {
      vt.len = v.len
    }
  }
  tau = rho
  K = length(v)
  y = pmax(y-tau,0)
  return(y)
}

Prune <- function(U,r,q,e){
  l = apply(U,1,norm,"2")
  quan.l = quantile(l,1-q)
  P0 = which(l>=quan.l)
  U0 = U[P0,]
  x = rowMeans(FNN::knnx.dist(U,U0,r))
  quan.x = quantile(x,1-e)
  P = which(x>quan.x)
  P = P0[P]
  return(P)
}

## find out the set S corresponding to unpruned matrix
.find_act <- function(s,P){
  s+sum(P<s)
}

## use direct projection
## A is the response matrix, X is the covariate matrix, alpha is the tuning param
gom.svd <- function(A,K,r=10,q=.25,e=.02,eps=1e-3){
  svd.A = RSpectra::svds(A,K)
  P = Prune(svd.A$u,r,q,e)
  IK = diag(rep(1,K))
  S = c()
  U = svd.A$u
  Y0 = U[-P,]
  Y = Y0
  for (k in 1:K) {
    l = apply(Y,1,norm,"2")
    S = c(S,which.max(l))
    u = Y[S[k],]/norm(Y[S[k],],"2")
    u = data.matrix(u)
    Y = Y%*%(IK-u%*%t(u))
  }
  US = Y0[S,]
  Pi = U%*%solve(US)
  Pi = t(apply(Pi,1,.proj))
  Theta = svd.A$v%*%(svd.A$d*t(svd.A$u))%*%(Pi%*%solve(t(Pi)%*%Pi))
  Theta = pmin(Theta,1-eps)
  Theta = pmax(Theta,eps)
  
  S = lapply(S,.find_act,P)
  return(list(Pi=Pi,Theta=Theta,S=S))
}

## covariate assisted spectral clustering
gom.cov <- function(A,X,K,alpha,r=10,q=.25,e=.02,eps=1e-3){
  L = A%*%t(A)+alpha*X%*%t(X)
  L = L-diag(diag(L))
  evd.A = RSpectra::eigs(L,K)
  P = Prune(evd.A$vectors,r,q,e)
  IK = diag(rep(1,K))
  U = evd.A$vectors
  S = c()
  Y0 = U[-P,]
  Y = Y0
  for (k in 1:K) {
    l = apply(Y,1,norm,"2")
    S = c(S,which.max(l))
    u = Y[S[k],]/norm(Y[S[k],],"2")
    u = data.matrix(u)
    Y = Y%*%(IK-u%*%t(u))
  }
  US = Y0[S,]
  Pi = U%*%solve(US)
  Pi = pmax(Pi,0)
  Pi = t(apply(Pi,1,.proj))
  
  svd.A = RSpectra::svds(A,K)
  Theta = svd.A$v%*%(svd.A$d*t(svd.A$u))%*%(Pi%*%solve(t(Pi)%*%Pi))
  Theta = pmin(Theta,1-eps)
  Theta = pmax(Theta,eps)
  
  S = lapply(S,.find_act,P)
  return(list(Pi=Pi,Theta=Theta,S=S))
}

## To find the best parameter using random missing and recovery performance
find_best_alpha <- function(alpha_seq,A,X,K,r=10,q=.25,e=.02,eps=1e-3,nfold=5){
  err = sapply(alpha_seq,.recovery.cv,A,X,K,r,q,e,eps,nfold)
  model.loess = stats::loess(err~alpha_seq)
  idx.best = which.min(model.loess$fitted)
  return(alpha_seq[idx.best])
}

.recovery.cv <- function(alpha,A,X,K,r=10,q=0.25,e=0.02,eps=1e-3,nfold=5){
  ## create K-fold missing matrix
  n = nrow(A)
  idx.A = rep(1:nfold,length.out=nrow(A)*ncol(A))
  idx.X = rep(1:nfold,length.out=nrow(X)*ncol(X))
  idx.A = sample(idx.A)
  idx.X = sample(idx.X)
  err = rep(0,nfold)
  
  for (kk in 1:nfold) {
    A.tmp = A
    X.tmp = X
    bool_mat.A = matrix(idx.A == kk,nrow=n)
    bool_mat.X = matrix(idx.X == kk,nrow=n)
    A.tmp[as.logical(bool_mat.A)] = 0
    X.tmp[as.logical(bool_mat.X)] = 0
    cov.tmp = gom.cov(A.tmp,X.tmp,K,alpha,r,q,e,eps)
    A.est = cov.tmp$Pi%*%t(cov.tmp$Theta)
    err[kk] = mean(abs(A[bool_mat.A]-A.est[bool_mat.A]))
  }
  return(mean(err))
}

## To align the column permutation in simulation, using Pi as input
## return best column permutation
## It is not efficient, but most accurate for simulation
## Hungarian algorithm (in package "clue") can be of help, but returns wrong result sometimes
find_best_idx <- function(K,Pi,Pi.r){
  idx.all = gtools::permutations(K,K)
  err = apply(idx.all,1,.Pi_diff,Pi,Pi.r)
  idx = which.min(err)
  return(idx.all[idx,])
}

.Pi_diff <- function(idx,Pi,Pi.r){
  mean(abs(Pi.r-Pi[,idx]))
}

## self written gom.jml, can be ignored
.gom.self.proc <- function(A,K){
  A.resp = 1-is.na(A)
  N = nrow(A)
  J = ncol(A)
  score = rowSums(A)/rowSums(A.resp)
  theta0 = seq(-1.5,1.5,len=K)
  Pi = t(exp(-(outer(stats::qlogis((score+.1)/1.2),theta0,"-"))^2)) # may have problems
  Pi = Pi/rowSums(Pi)
  p.item = colSums(A)/colSums(A.resp)
  theta0 = seq(-2,2,len=K)
  Theta = stats::plogis(outer(theta0,-stats::qlogis(p.item),"-")) # may need to transpose
  res = list(Pi=Pi,Theta=Theta)
  return(res)
}

.gom.self.deviance <- function(Pi,Theta,K,A){
  prob0 = prob1 = 0
  for (k in 1:K) {
    prob1 = prob1 + outer(Pi[,k],Theta[,k])
    prob0 = prob0 + outer(Pi[,k],1-Theta[,k])
  }
  dev = -2*sum(log(A*t(prob1)+(1-A)*(prob0)))
  return(dev)
}

gom.jml.self <- function(A,K,max.iter=600,eps=1e-3){
  datproc = .gom.self.proc(A,K)
  Pi = datproc$Pi
  Theta = datproc$Theta
  iter = 0
  conv = 10
  devchange = -100
  dev = .gom.self.deviance(Pi,Theta,K,A)
  while (((eps<devchange)|(eps<conv))&(iter<max.iter)) {
    Pi.t = Pi
    Theta.t = Theta
    dev.t = dev
    
    Pi = Pi.t* ((A/(Pi.t%*%t(Theta)+1e-14))%*%Theta.t 
                + ((1-A)/(Pi.t%*%t(1-Theta)+1e-14))%*%(1-Theta.t))/ncol(A)
    Theta.rev = 1-Theta.t
    Theta = Theta.t*((A/(Pi%*%t(Theta.t)+1e-14))%*%Pi)
    Theta.rev = Theta.rev*(((1-A)/(Pi%*%t(Theta.rev)+1e-14))%*%Pi)
    Theta = Theta/(Theta+Theta.rev)
    Theta = pmin(Theta,1-eps)
    Theta = pmax(Theta,eps)
    
    dev = .gom.self.deviance(Pi,Theta,K,A)
    conv = max(max(abs(Pi-Pi.t)),max(abs(Theta-Theta.t)))
    devchange = abs((dev-dev.t)/dev)
    
    iter = iter + 1
  }
  return(list(Pi=Pi,Theta=Theta,iter=iter))
}

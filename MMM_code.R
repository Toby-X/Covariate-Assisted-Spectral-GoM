library(RSpectra)
library(Rfast)
library(Matrix)
library(gtools)
library(sirt)
library(future.apply)
library(FNN)
plan(multisession,workers=7)

Prune <- function(U,r,q,e){
  l = apply(U,1,norm,"2")
  quan.l = quantile(l,1-q)
  P0 = which(l>=quan.l)
  U0 = U[P0,]
  x = rowmeans(knnx.dist(U,U0,r))
  quan.x = quantile(x,1-e)
  P = which(x>quan.x)
  P = P0[P]
  return(P)
}

find_act <- function(s,P){
  s+sum(P<s)
}

null_cluster <- function(A,K,r=10,q=.4,e=.2,eps=1e-3){
  svd.A = svds(A,K)
  P = Prune(svd.A$u,r,q,e)
  IK = Diag.matrix(K,1)
  S = c()
  U = svd.A$u
  Y0 = U[-P,]
  Y = Y0
  for (k in 1:K) {
    l = apply(Y,1,norm,"2")
    S = c(S,which.max(l))
    u = Y[S[k],]/norm(Y[S[k],],"2")
    u = data.matrix(u)
    Y = mat.mult(Y,IK-mat.mult(u,transpose(u)))
  }
  US = Y0[S,]
  Pi = mat.mult(U,solve(US))
  Pi = pmax(Pi,0)
  Pi = mat.mult(diag((rowsums(Pi)+1e-14)^(-1)),Pi)
  Theta = mat.mult(mat.mult(mat.mult(mat.mult(svd.A$v,diag(svd.A$d)),transpose(svd.A$u)),Pi),spdinv(mat.mult(transpose(Pi),Pi)))
  Theta = pmin(Theta,1-eps)
  Theta = pmax(Theta,eps)
  
  S = lapply(S,find_act,P)
  return(list(Pi=Pi,Theta=Theta,S=S))
}

cov_cluster <- function(A,X,K,alpha,r=10,q=.1,e=.2,eps=1e-3){
  L = mat.mult(A,transpose(A))+alpha*mat.mult(X,transpose(X))
  L = Diag.fill(L,0)
  evd.A = eigs(L,K)
  P = Prune(evd.A$vectors,r,q,e)
  IK = Diag.matrix(K,1)
  U = evd.A$vectors
  S = c()
  Y0 = U[-P,]
  Y = Y0
  for (k in 1:K) {
    l = apply(Y,1,norm,"2")
    S = c(S,which.max(l))
    u = Y[S[k],]/norm(Y[S[k],],"2")
    u = data.matrix(u)
    Y = mat.mult(Y,IK-mat.mult(u,transpose(u)))
  }
  US = Y0[S,]
  Pi = mat.mult(U,solve(US))
  Pi = pmax(Pi,0)
  Pi = mat.mult(diag((rowsums(Pi)+1e-14)^(-1)),Pi)
  
  svd.A = svds(A,K)
  Theta = mat.mult(mat.mult(mat.mult(mat.mult(svd.A$v,diag(svd.A$d)),transpose(svd.A$u)),Pi),spdinv(mat.mult(transpose(Pi),Pi)))
  Theta = pmin(Theta,1-eps)
  Theta = pmax(Theta,eps)
  
  S = lapply(S,find_act,P)
  return(list(Pi=Pi,Theta=Theta,S=S))
}

find_best_idx <- function(idx,Pi,Pi.r,type="1"){
  norm(Pi.r-Pi[,idx],type = type)
}

find_best_alpha <- function(alpha,A,X,K,type="2",r=10,q=0.4,e=0.2,eps=1e-3,nfold=5){
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
    cov.tmp = cov_cluster(A.tmp,X.tmp,K,alpha,r,q,e,eps)
    A.est = mat.mult(cov.tmp$Pi,transpose(cov.tmp$Theta))
    err[kk] = median(abs(A[bool_mat.A]-A.est[bool_mat.A]))
  }
  return(mean(err))
}

.gom.self.proc <- function(A,K){
  A.resp = 1-is.na(A)
  N = nrow(A)
  J = ncol(A)
  score = rowsums(A)/rowsums(A.resp)
  theta0 = seq(-1.5,1.5,len=K)
  Pi = transpose(exp(-(Outer(stats::qlogis((score+.1)/1.2),theta0,"-"))^2))
  Pi = Pi/rowsums(Pi)
  p.item = colsums(A)/colsums(A.resp)
  theta0 = seq(-2,2,len=K)
  Theta = stats::plogis(Outer(theta0,-stats::qlogis(p.item),"-"))
  res = list(Pi=Pi,Theta=Theta)
  return(res)
}

.gom.self.deviance <- function(Pi,Theta,K,A){
  prob0 = prob1 = 0
  for (k in 1:K) {
    prob1 = prob1 + Outer(Pi[,k],Theta[,k])
    prob0 = prob0 + Outer(Pi[,k],1-Theta[,k])
  }
  dev = -2*sum(log(A*transpose(prob1)+(1-A)*transpose(prob0)))
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
    
    Pi = Pi.t*(mat.mult(A/(mat.mult(Pi.t,transpose(Theta.t))+1e-14),Theta.t)+
                 mat.mult((1-A)/(mat.mult(Pi.t,transpose(1-Theta.t))+1e-14),1-Theta.t))/ncol(A)
    Pi = pmax(Pi,eps)
    Pi = pmin(Pi,1-eps)
    Theta.rev = 1-Theta.t
    Theta = Theta.t*mat.mult(transpose(A/(mat.mult(Pi,transpose(Theta.t))+1e-14)),Pi)
    Theta.rev = Theta.rev*mat.mult(transpose((1-A)/(mat.mult(Pi,transpose(Theta.rev))+1e-14)),Pi)
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

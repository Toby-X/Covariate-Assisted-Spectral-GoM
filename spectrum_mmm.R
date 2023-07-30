#-*- coding:utf-8 -*-
library(RSpectra)
library(Rfast)
library(Matrix)
library(gtools)
library(sirt)
library(future.apply)
library(FNN)
plan(multisession)

## data generation
K = 3
N = 500
p = N/5
R = 2
idx.all = permutation(1:K)

pi_gen <- function(pi,K){
  pi = rdirichlet(1,rep(1,K))
  return(pi)
}

set.seed(230627)
Pi = matrix(rep(0,N*K),nrow=N)
Pi = t(apply(Pi,1,pi_gen,K=K))
Theta = matrix(runif(p*K),ncol=K)
M = matrix(rnorm(R*K),nrow=R)
Pi[1:K,] = diag(rep(1,K))

A_t = Pi%*%t(Theta)
X_t = Pi%*%t(M)

A = matrix(rbinom(N*p,1,A_t),nrow = N)
X = X_t + matrix(rnorm(N*R,0,.7),nrow=N)

## tuning parameter (The choice is still undecided)
alpha = .1

## functions
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

null_cluster <- function(A,K,r,q,e,eps=1e-3){
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

cov_cluster <- function(A,X,K,alpha,r,q,e,eps=1e-3){
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
    # err[kk] = mean(abs(A[bool_mat.A]-A.est[bool_mat.A]))
    err[kk] = norm(A-A.est,"2")
  }
  return(mean(err))
}

find_best_alpha2 <- function(alpha,A,X,K,type="2",r=10,q=0.4,e=0.2,eps=1e-3,nfold=5){
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
    U.est = svds(A.est,K)$u
    U = svds(A,K)$u
    # err[kk] = mean(abs(A[bool_mat.A]-A.est[bool_mat.A]))
    # err[kk] = norm(U-U.est,"2")
    err[kk] = mean(abs(U-U.est))
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
  datproc = .gom.proc(A,K)
  Pi = datproc$Pi
  Theta = datproc$Theta
  iter = 0
  conv = 10
  devchange = -100
  dev = .gom.deviance(Pi,Theta,K,A)
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
    
    dev = .gom.deviance(Pi,Theta,K,A)
    conv = max(max(abs(Pi-Pi.t)),max(abs(Theta-Theta.t)))
    devchange = abs((dev-dev.t)/dev)
    
    iter = iter + 1
  }
  return(list(Pi=Pi,Theta=Theta,iter=iter))
}

## diag deletion
r = 10
q = .4
e = .2
tic()
res.null = null_cluster(A,K,r,q,e)
toc()
idx.null = apply(idx.all,1,find_best_idx,Pi=res.null$Pi,Pi.r=Pi)
idx.null = which.min(idx.null)
idx.null = idx.all[idx.null,]
l1.null = norm(Pi-res.null$Pi[,idx.null],"1")
l2.null = norm(Pi-res.null$Pi[,idx.null],"2")
linfty.null = max(abs(Pi-res.null$Pi[,idx.null]))
l1.null
l2.null
linfty.null

res.X = null_cluster(X,K,r,q,e)
norm(res.X$Pi[,c(2,3,1)]-Pi,"2")

alpha = .3
res.cov = cov_cluster(A,X,K,alpha,r,q,e)
idx.cov = apply(idx.all,1,find_best_idx,Pi=res.cov$Pi,Pi.r=Pi)
idx.cov = idx.all[which.min(idx.cov),]
l1.cov = norm(Pi-res.cov$Pi[,idx.cov],"1")
l2.cov = norm(Pi-res.cov$Pi[,idx.cov],"2")
linfty.cov = max(abs(Pi-res.cov$Pi[,idx.cov]))
l1.cov
l2.cov
linfty.cov

norm(res.cov$Theta-res.null$Theta,"1")
norm(res.cov$Theta-res.null$Theta,"2")
norm(res.cov$Theta-res.null$Theta,"I")
norm(res.cov$Theta-res.null$Theta,"F")
norm(res.cov$Theta-res.null$Theta,"M")

max(apply(res.cov$Pi-res.null$Pi,1,norm,"2"))
median(abs(res.cov$Pi-res.null$Pi))

res.cov.tmp = cov_cluster(A,X,K,.55,r,q,e)
median(abs(res.cov.tmp$Pi-res.null$Pi))
norm(res.cov.tmp$Theta-res.null$Theta,"m")-norm(res.cov$Theta-res.null$Theta,"m")
norm(res.cov.tmp$Pi-res.null$Pi,"m")-norm(res.cov$Pi-res.null$Pi,"m")
mean(abs(res.cov.tmp$Theta-res.null$Theta))-median(abs(res.cov$Theta-res.null$Theta))

tic()
res.jml = gom.jml.self(A,K)
toc()
idx.jml = apply(idx.all,1,find_best_idx,Pi=res.jml$Pi,Pi.r=Pi)
idx.jml = which.min(idx.jml)
idx.jml = idx.all[idx.jml,]
l1.jml = norm(Pi-res.jml$Pi[,idx.jml],"1")
l2.jml = norm(Pi-res.jml$Pi[,idx.jml],"2")
linfty.jml = max(abs(Pi-res.jml$Pi[,idx.jml]))
# linfty.jml = norm(Pi-res.jml$Pi[,idx.jml],"I")
l1.jml
l2.jml
linfty.jml

# res.jml1 = sirt::gom.jml(data.frame(A),K)
# idx.jml1 = apply(idx.all,1,find_best_idx,Pi=res.jml1$g,Pi.r=Pi)
# idx.jml1 = which.min(idx.jml1)
# idx.jml1 = idx.all[idx.jml1,]
# l1.jml1 = norm(Pi-res.jml1$g[,idx.jml1],"1")
# l2.jml1 = norm(Pi-res.jml1$g[,idx.jml1],"2")
# linfty.jml1 = max(abs(Pi-res.jml1$g[,idx.jml1]))
# l1.jml1
# l2.jml1
# linfty.jml1

tic()
AA = mat.mult(A,transpose(A))
AA = Diag.fill(AA,0)
XX = mat.mult(X,transpose(X))
XX = Diag.fill(XX,0)
AA.evd = eigs(AA,K+1)
XX.evd = eigs(XX,K+1)
alpha_min = (AA.evd$values[K]-AA.evd$values[K+1])/XX.evd$values[1]
if (R<K){
  # alpha_max = (AA.evd$values[1])/(XX.evd$values[R])
  alpha_max = AA.evd$values[K]/XX.evd$values[R]
} else {
  alpha_max = (AA.evd$values[K])/(XX.evd$values[K]-XX.evd$values[K+1])
}

# alpha_max = (AA.evd$values[K]-AA.evd$values[K+1])/XX.evd$values[R]
alpha_seq = seq(from=0,to=alpha_max,length=100)# fail sometimes
alpha_seq1 = seq(from=0,to=alpha_max,length=100)
# alpha_seq = log(seq(from=exp(.00001),to=exp(0.002),length=100))
# alpha_seq1 = seq(from=.05,to=1,length=100)
# within_var = future_sapply(alpha_seq,find_best_alpha,A,X,3,100,future.seed=T)
# within_var = sapply(alpha_seq,find_best_alpha,A,X,3,100)
toc()
plot(alpha_seq,within_var)
alpha = 0.5# K=3
alpha = 0.75#K=8

alpha = 0.90
tic()
res.cov = cov_cluster(A,X,K,alpha,r,q,e)
toc()
l1.cov = apply(idx.all,1,find_best_idx,Pi=res.cov$Pi,Pi.r=Pi)
idx.cov = idx.all[which.min(l1.cov),]
l1.cov = min(l1.cov)
l2.cov = norm(Pi-res.cov$Pi[,idx.cov],"2")
linfty.cov = max(abs(Pi-res.cov$Pi[,idx.cov]))
l1.cov
l2.cov
linfty.cov

com_accuracy_l1 <- function(alpha){
  res.cov = cov_cluster(A,X,K,alpha,r,q,e,eps=1e-3)
  l1.cov = apply(idx.all,1,find_best_idx,Pi=res.cov$Pi,Pi.r=Pi)
  return(min(l1.cov))
}

com_accuracy_l2 <- function(alpha){
  res.cov = cov_cluster(A,X,K,alpha,r,q,e,eps=1e-3)
  l2.cov1 = apply(idx.all,1,find_best_idx,Pi=res.cov$Pi,Pi.r=Pi,type="2")
  # l2.cov2 = norm(Theta-res.cov$Theta[,idx.all[which.min(l2.cov1)]],"2")
  # l2.cov = min(l2.cov1)+l2.cov2
  l2.cov = min(l2.cov1)
  return(l2.cov)
}

com_accuracy_l2_all <- function(alpha){
  res.cov = cov_cluster(A,X,K,alpha,r,q,e,eps=1e-3)
  l2.cov1 = apply(idx.all,1,find_best_idx,Pi=res.cov$Pi,Pi.r=Pi,type="2")
  l2.cov2 = norm(Theta-res.cov$Theta[,idx.all[which.min(l2.cov1)]],"2")
  l2.cov = min(l2.cov1)+l2.cov2
  # l2.cov = min(l2.cov1)
  return(l2.cov)
}

com_accuracy_l2_A <- function(alpha){
  res.cov = cov_cluster(A,X,K,alpha,r,q,e,eps=1e-3)
  A.est = mat.mult(res.cov$Pi,transpose(res.cov$Theta))
  # l2.cov1 = apply(idx.all,1,find_best_idx,Pi=res.cov$Pi,Pi.r=Pi,type="2")
  # l2.cov2 = norm(Theta-res.cov$Theta[,idx.all[which.min(l2.cov1)]],"2")
  l2.cov = norm(A-A.est,"2")
  # l2.cov = min(l2.cov1)
  return(l2.cov)
}

res.cv = future_sapply(alpha_seq,com_accuracy_l2,future.seed=T)
plot(alpha_seq,res.cv,"l")
res.cv1 = future_sapply(alpha_seq,com_accuracy_l2_all,future.seed=T)
plot(alpha_seq,res.cv1,"l")
res.cv2 = future_sapply(alpha_seq,com_accuracy_l2_A,future.seed=T)
plot(alpha_seq,res.cv2,"l")
res.cv2 = future_sapply(alpha_seq1,com_accuracy_l2_all,future.seed=T)
plot(alpha_seq1,res.cv2,"l")
res.cv1 = future_sapply(alpha_seq1,com_accuracy_l2,future.seed = T)
plot(alpha_seq1,res.cv1,"l")
res.cv2 = future_sapply(alpha_seq,com_accuracy_l2,future.seed = T)
plot(alpha_seq,res.cv2,"l")
plot(alpha_seq,within_var,"l")

# unable to determine alpha from existing method !??

dat = cbind(within_var,res.cv,res.cv2)
std_max = function(l){
  return((l-min(l))/(max(l)-min(l)))
}
dat = apply(dat,2,std_max)
dat = data.frame(cbind(alpha_seq,dat))
colnames(dat) = c("alpha","Var","L1","L2")
dat = dat %>% pivot_longer("Var":"L2",names_to = "type", values_to = "value")
ggplot(data=dat)+
  geom_line(aes(x=alpha,y=value,col=type))

err = future_sapply(alpha_seq,find_best_alpha,A,X,K,future.seed=T)
plot(alpha_seq,err,"l",col="red")
err1 = future_sapply(alpha_seq,find_best_alpha,A,X,K,future.seed=T)
err2 = future_sapply(alpha_seq,find_best_alpha,A,X,K,future.seed=T)
err3 = future_sapply(alpha_seq,find_best_alpha,A,X,K,future.seed=T)
plot(alpha_seq,err2,"l",col="red")
library(ggplot2)
ggplot()+
  geom_smooth(aes(alpha_seq,err3))+
  geom_line(aes(alpha_seq,err3))

#-*- coding:utf-8 -*-
library(RSpectra)
library(Rfast)
library(Matrix)
library(gtools)
library(sirt)
library(tictoc)

## data generation
K = 3
N = 500
p = N/5
R = 2
Pi = matrix(rep(0,N*K),nrow=N)
Theta = matrix(rep(0,p*K),ncol=K)
M = matrix(rnorm(R*K),nrow=R)
idx.all = permutation(1:K)

pi_gen <- function(pi,K){
  pi = rdirichlet(1,rep(1,K))
  return(pi)
}

theta_gen <- function(theta,K){
  theta = runif(K)
  return(theta)
}

set.seed(230625)
Pi = t(apply(Pi,1,pi_gen,K=K))
Theta = t(apply(Theta,1,theta_gen,K=K))
Pi[1:K,] = diag(rep(1,K))

A_t = Pi%*%t(Theta)
X_t = Pi%*%t(M)

A = matrix(rbinom(N*p,1,A_t),nrow = N)
X = X_t + matrix(rnorm(N*R,0,1),nrow=N)

## tuning parameter (The choice is still undecided)
alpha = .1

## functions
find_knn_x <- function(vec,U,r){
  U1 = transpose(transpose(U)-vec)
  d = apply(U1,1,norm,"2")
  d = sort(d)
  x = mean(d[1:(r+1)])
  return(x)
}

Prune <- function(U,r,q,e){
  l = apply(U,1,norm,"2")
  quan.l = quantile(l,1-q)
  P0 = which(l>=quan.l)
  U0 = U[P0,]
  x = apply(U0,1,find_knn_x,U=U,r=r)
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
  Pi = mat.mult(diag(rowsums(Pi)^(-1)),Pi)
  Theta = mat.mult(mat.mult(mat.mult(mat.mult(svd.A$v,diag(svd.A$d)),transpose(svd.A$u)),Pi),spdinv(mat.mult(transpose(Pi),Pi)))
  Theta = pmax(Theta,1-eps)
  Theta = pmin(Theta,eps)
  
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
  Pi = mat.mult(diag(rowsums(Pi)^(-1)),Pi)
  
  svd.A = svds(A,K)
  Theta = mat.mult(mat.mult(mat.mult(mat.mult(svd.A$v,diag(svd.A$d)),transpose(svd.A$u)),Pi),spdinv(mat.mult(transpose(Pi),Pi)))
  Theta = pmax(Theta,1-eps)
  Theta = pmin(Theta,eps)
  
  S = lapply(S,find_act,P)
  return(list(Pi=Pi,Theta=Theta,S=S))
}

## diag deletion
r = 10
q = .05
e = .05
tic()
res.null = null_cluster(A,K,r,q,e)
toc()
idx.null = idx.all[1,]
l1.null = norm(Pi-res.null$Pi[,idx.all[1,]],"1")
for (i in 2:nrow(idx.all)) {
  l1.tmp = norm(Pi-res.null$Pi[,idx.all[i,]],"1")
  if (l1.tmp < l1.null){
    l1.null = l1.tmp
    idx.null = idx.all[i,]
  }
}
l2.null = norm(Pi-res.null$Pi[,idx.null],"2")
linfty.null = max(abs(Pi-res.null$Pi[,idx.null]))
l1.null
l2.null
linfty.null

alpha = .1
res.cov = cov_cluster(A,X,K,alpha,r,q,e)
idx.cov = idx.all[1,]
l1.cov = norm(Pi-res.cov$Pi[,idx.all[1,]],"1")
for (i in 2:nrow(idx.all)) {
  l1.tmp = norm(Pi-res.cov$Pi[,idx.all[i,]],"1")
  if (l1.tmp < l1.cov){
    l1.cov = l1.tmp
    idx.cov = idx.all[i,]
  }
}
l2.cov = norm(Pi-res.cov$Pi[,idx.cov],"2")
linfty.cov = max(abs(Pi-res.cov$Pi[,idx.cov]))
l1.cov
l2.cov
linfty.cov

tic()
res.jml = gom.jml(data.frame(A),K)
toc()
idx.jml = idx.all[1,]
l1.jml = norm(Pi-res.jml$g[,idx.all[1,]],"1")
for (i in 2:nrow(idx.all)) {
  l1.tmp = norm(Pi-res.jml$g[,idx.all[i,]],"1")
  if (l1.tmp < l1.jml){
    l1.jml = l1.tmp
    idx.jml = idx.all[i,]
  }
}
l2.jml = norm(Pi-res.jml$g[,idx.jml],"2")
linfty.jml = max(abs(Pi-res.jml$g[,idx.jml]))
l1.jml
l2.jml
linfty.jml

tic()
AA = mat.mult(A,transpose(A))
AA = Diag.fill(AA,0)
XX = mat.mult(X,transpose(X))
XX = Diag.fill(XX,0)
AA.svd = svds(AA,K+1)
XX.svd = svds(XX,K+1)
alpha_min = (AA.svd$d[K]-AA.svd$d[K+1])/XX.svd$d[1]
if (R<K){
  alpha_max = (AA.svd$d[1])/(XX.svd$d[R])
} else {
  alpha_max = (AA.svd$d[1])/(XX.svd$d[K]-XX.svd$d[K+1])
}
alpha_seq = seq(from=alpha_min,to=alpha_max,length=100)# fail sometimes
alpha_seq = seq(from=.05,to=1,length=100)
within_var = rep(0,length(alpha_seq))
for (i in 1:length(alpha_seq)) {
  L.cov = mat.mult(A,transpose(A))+alpha_seq[i]*mat.mult(X,transpose(X))
  L.all = Diag.fill(L.cov,0)
  evd.all = eigs(L.all,k=K)
  kmeans.all = kmeans(evd.all$vectors,K,algorithm = "Lloyd", nstart = 2,iter.max = 100)
  within_var[i] = kmeans.all$tot.withinss
}
toc()
plot(alpha_seq,within_var)
alpha = 0.15

tic()
res.cov = cov_cluster(A,X,K,alpha,r,q,e)
toc()
idx.cov = idx.all[1,]
l1.cov = norm(Pi-res.cov$Pi[,idx.all[1,]],"1")
for (i in 2:nrow(idx.all)) {
  l1.tmp = norm(Pi-res.cov$Pi[,idx.all[i,]],"1")
  if (l1.tmp < l1.cov){
    l1.cov = l1.tmp
    idx.cov = idx.all[i,]
  }
}
l2.cov = norm(Pi-res.cov$Pi[,idx.cov],"2")
linfty.cov = max(abs(Pi-res.cov$Pi[,idx.cov]))
l1.cov
l2.cov
linfty.cov

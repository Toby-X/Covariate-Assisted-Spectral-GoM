#-*- coding:utf-8 -*-
library(Boom)
library(RSpectra)
library(Rfast)

## data generation
K = 3
N = 500
p = N/50
r = 2
Pi = matrix(rep(0,N*K),nrow=N)
Theta = matrix(rep(0,p*K),ncol=K)
M = matrix(rnorm(r*K),nrow=r)
sp_point = (1:K)/K

pi_gen <- function(pi,K){
  u = runif(1)
  for (i in 1:K) {
    if (u < sp_point[i]){
      pi[i] = 1
      break
    }
  }
  return(pi)
}

theta_gen <- function(theta,K){
  theta = runif(K)
  return(theta)
}

set.seed(230625)
Pi = t(apply(Pi,1,pi_gen,K=K))
Theta = t(apply(Theta,1,theta_gen,K=K))

A_t = Pi%*%t(Theta)
X_t = Pi%*%t(M)

A = matrix(rbinom(N*p,1,A_t),nrow = N)
X = X_t + matrix(rnorm(N*r,0,1),nrow=N)

## tuning parameter
alpha = .1

## Null
L.null = mat.mult(A,transpose(A))
svd.null = svds(L.null,k=K)
kmeans.null = kmeans(svd.null$u,K,algorithm = "Lloyd", nstart=10,iter.max = 100)
result.null = apply(Pi,1,which.max)
for (i in 1:N) {
  if (result.null[i] == 1){
    result.null[i] = 2
  } else if (result.null[i] == 2){
    result.null[i] = 1
  }
}
mean(kmeans.null$cluster==result.null)

## Diagonal Deletion
L.del = Diag.fill(L.null)
svd.del = eigs(L.del,k=K)
kmeans.del = kmeans(svd.del$vectors,K,algorithm = "Lloyd", nstart=10,iter.max = 100)
result.del = apply(Pi,1,which.max)
mean(kmeans.del$cluster==result.del)

## Covariate
L.cov = mat.mult(A,transpose(A))+alpha*mat.mult(X,transpose(X))
evd.cov = eigs(L.cov,k=K)
kmeans.cov = kmeans(evd.cov$vectors,K,algorithm = "Lloyd", nstart = 10,iter.max = 100)
result.cov = apply(Pi,1,which.max)
for (i in 1:N) {
  if (result.cov[i] == 2){
    result.cov[i] = 1
  } else if (result.cov[i] == 1){
    result.cov[i] = 2
  } else{
    result.cov[i] = 3
  }
}
mean(kmeans.cov$cluster==result.cov)

## Concatenate
L.conc = cbind(A,X)
svd.conc = svds(L.conc,K)
kmeans.conc = kmeans(svd.conc$u,K,algorithm = "Lloyd",nstart=10,iter.max=100)
result.conc = apply(Pi,1,which.max)
for (i in 1:N) {
  if (result.conc[i] == 2){
    result.conc[i] = 1
  } else if (result.conc[i] == 1){
    result.conc[i] = 3
  } else {
    result.conc[i] = 2
  }
}
mean(kmeans.conc$cluster==result.conc)

## Diagonal Deletion + Covariate
L.cov = mat.mult(A,transpose(A))+alpha*mat.mult(X,transpose(X))
L.all = Diag.fill(L.cov,0)
evd.all = eigs(L.all,k=K)
kmeans.all = kmeans(evd.all$vectors,K,algorithm = "Lloyd", nstart = 10,iter.max = 100)
result.all = apply(Pi,1,which.max)
for (i in 1:N) {
  if (result.all[i] == 2){
    result.all[i] = 2
  } else if (result.all[i] == 1){
    result.all[i] = 3
  } else{
    result.all[i] = 1
  }
}
mean(kmeans.all$cluster==result.all)

# superiority appears when r is small

## Diagonal Deletion for X
LX = mat.mult(X,transpose(X))
LX = Diag.fill(LX,0)
evd.X = eigs(LX,K)
kmeans.X = kmeans(evd.X$vectors,K,algorithm = "Lloyd", nstart = 10, iter.max = 100)
result.X = apply(Pi,1,which.max)
for (i in 1:N) {
  if (result.X[i] == 2){
    result.X[i] = 1
  } else if (result.X[i] == 1){
    result.X[i] = 3
  } else{
    result.X[i] = 2
  }
}
mean(kmeans.X$cluster==result.X)

## parameter tuning
AA = mat.mult(A,transpose(A))
AA = Diag.fill(AA,0)
XX = mat.mult(X,transpose(X))
XX = Diag.fill(XX,0)
AA.svd = svds(AA,K+1)
XX.svd = svds(XX,K+1)
alpha_min = (AA.svd$d[K]-AA.svd$d[K+1])/XX.svd$d[1]
if (r<K){
  alpha_max = (AA.svd$d[1])/(XX.svd$d[r])
} else {
  alpha_max = (AA.svd$d[1])/(XX.svd$d[K]-XX.svd$d[K+1])
}
alpha_seq = seq(from=alpha_min,to=alpha_max,length=500)
within_var = rep(0,length(alpha_seq))
for (i in 1:length(alpha_seq)) {
  L.cov = mat.mult(A,transpose(A))+alpha_seq[i]*mat.mult(X,transpose(X))
  L.all = Diag.fill(L.cov,0)
  evd.all = eigs(L.all,k=K)
  kmeans.all = kmeans(evd.all$vectors,K,algorithm = "Lloyd", nstart = 10,iter.max = 100)
  within_var[i] = kmeans.all$tot.withinss
}
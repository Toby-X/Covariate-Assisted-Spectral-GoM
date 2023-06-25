#-*- coding:utf-8 -*-
library(Boom)

## data generation
K = 3
N = 500
p = N/5
r = 10
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
X = X_t + matrix(rnorm(N*r,0,.1),nrow=N)

## tuning parameter
alpha = .1

## Null
L.null = A%*%t(A)+

## Diagonal Deletion

## Covariate

## Diagonal Deletion + Covariate
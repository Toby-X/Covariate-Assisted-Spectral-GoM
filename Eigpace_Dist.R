library(gtools)
library(Rfast)
library(pracma)

## data generation
K = 3
N = 2000
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

alpha_seq = seq(from=0,to=1.5,length=300)

AA = mat.mult(A,transpose(A))
AA = Diag.fill(AA,0)
XX = mat.mult(X,transpose(X))
XX = Diag.fill(XX,0)

L_gen <- function(alpha){
  mat.mult(A_t,transpose(A_t))+alpha*mat.mult(X_t,transpose(X_t))
}

U_dist <- function(alpha,K){
  L_t = L_gen(alpha)
  L = AA+alpha*XX
  U_t = eigs(L_t,K)$vectors
  U = eigs(L,K)$vectors
  UTU = mat.mult(transpose(U),U_t)
  sigma1 = svd(UTU)$d[1]
  dist = 2*(K-sigma1)
  return(dist)
}

err = future_sapply(alpha_seq,U_dist,K,future.seed=T)
plot(alpha_seq,err)

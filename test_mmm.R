#-*- coding:utf-8 -*-
source("MMM_code.R")

## data generation
K.all = c(3,8)
N.all = c(500,2000)
R.all = c(2,5)

pi_gen <- function(pi,K){
  pi = rdirichlet(1,rep(1,K))
  return(pi)
}

null.time = jml.time = cov.time = array(rep(0,10*length(K.all)*length(R.all)*length(N.all)),c(length(K.all),length(R.all),length(N.all),10))
metric.null = metric.jml = metric.cov = array(rep(0,10*length(K.all)*length(R.all)*length(N.all)),c(length(K.all),length(R.all),length(N.all),10))

for (iter.K in 1:length(K.all)) {
  K = K.all[iter.K]
  cat("K=",K,"\n")
  for (iter.R in 1:length(R.all)) {
    R = R.all[iter.R]
    cat("R=",R,"\n")
    for (iter.N in 1:length(N.all)) {
      N = N.all[iter.N]
      p = N/5
      cat("N=",N,"\n")
      for (l in 1:10) {
        set.seed(l)
        idx.all = permutation(1:K)
        Pi = matrix(rep(0,N*K),nrow=N)
        Pi = t(apply(Pi,1,pi_gen,K=K))
        Theta = matrix(runif(p*K),ncol=K)
        M = matrix(rnorm(R*K),nrow=R)
        Pi[1:K,] = diag(rep(1,K))
        
        A_t = Pi%*%t(Theta)
        X_t = Pi%*%t(M)
        
        A = matrix(rbinom(N*p,1,A_t),nrow = N)
        X = X_t + matrix(rnorm(N*R,0,.7),nrow=N)
        
        start.time = Sys.time()
        res.null = null_cluster(A,K)
        end.time = Sys.time()
        null.time[iter.K,iter.R,iter.N,l] = end.time-start.time
        idx.null = apply(idx.all,1,find_best_idx,Pi=res.null$Pi,Pi.r=Pi)
        idx.null = which.min(idx.null)
        idx.null = idx.all[idx.null,]
        metric.null[iter.K,iter.R,iter.N,l] = mean(abs(Pi-res.null$Pi[,idx.null]))
        
        start.time = Sys.time()
        res.jml = gom.jml(data.frame(A),K)
        end.time = Sys.time()
        jml.time[iter.K,iter.R,iter.N,l] = end.time-start.time
        idx.jml = apply(idx.all,1,find_best_idx,Pi=data.matrix(res.jml$g),Pi.r=Pi)
        idx.jml = which.min(idx.jml)
        idx.jml = idx.all[idx.jml,]
        metric.jml[iter.K,iter.R,iter.N,l] = mean(abs(Pi-res.jml$g[,idx.jml]))
        
        start.time = Sys.time()
        AA = mat.mult(A,transpose(A))
        AA = Diag.fill(AA,0)
        XX = mat.mult(X,transpose(X))
        XX = Diag.fill(XX,0)
        AA.evd = eigs(AA,K+1)
        XX.evd = eigs(XX,K+1)
        alpha_min = (AA.evd$values[K]-AA.evd$values[K+1])/XX.evd$values[1]
        if (R<K){
          alpha_max = AA.evd$values[K]/XX.evd$values[R]
        } else {
          alpha_max = (AA.evd$values[K])/(XX.evd$values[K]-XX.evd$values[K+1])
        }
        alpha_seq = seq(from=alpha_min,to=alpha_max,length=100)# fail sometimes
        err = future_sapply(alpha_seq,find_best_alpha,A,X,K,future.seed=T)
        alpha = alpha_seq[which.min(err)]
        res.cov = cov_cluster(A,X,K,alpha)
        end.time = Sys.time()
        cov.time[iter.K,iter.R,iter.N,l] = end.time-start.time
        l1.cov = apply(idx.all,1,find_best_idx,Pi=res.cov$Pi,Pi.r=Pi)
        idx.cov = idx.all[which.min(l1.cov),]
        metric.cov[iter.K,iter.R,iter.N,l] = mean(abs(Pi-res.cov$Pi[,idx.cov]))
      }
    }
  }
}

save.image("Simulation1.RData")

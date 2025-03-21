GeneData <- function(n=713, subj, time, group, position, Nij, R=2,L=3, X, lambda, A,U1,U2,U3, phi, tBmatrix_theta, miss_rate=0.1){
  
  N_obs <- length(time)
  ncX <- dim(X)[2]
  xidata <- mvrnorm(n=n,mu=rep(0,L),Sigma=diag(lambda,nrow=L,ncol=L)) #713*3
  #mean
  y <- rep(0,N_obs)
  Mt_tnsr <- ttl(as.tensor(A), list(U1=U1,U2=U2,U3=U3), 1:3)
  Mt <- Mt_tnsr@data
  eta0 <- rep(0,N_obs)
  for(s in 1:4){
    for(p in 1:7){
      index <- (group==s & position==p)
      eta0[index] <- tBmatrix_theta[index,] %*% Mt[s,p,]
    }
  }
  
  #variance
  varest <- rep(0,N_obs)
  for(i in 1:n){
    index <- subj==i
    varest[index] <-  phi[index,] %*% xidata[i,]
  }
  
  eta <- log(Nij)+eta0+varest
  mu <- exp(eta)
  y <- sapply(mu, function(x) {rpois(1,x)})
  
  data0 <- data.frame("subj"=subj, "time"=time, "y"=y, "group"=group, "position"=position, "Nij"=Nij)
  data0$subj <- as.factor(data0$subj)
  
  miss <- rep(0, N_obs)
  for(i in 1:n){
    index <- subj==i
    if(sum(index)==1){
      miss[index] <- 0
    }else{
      miss[index] <- rbinom(n=sum(index),size=1,prob=miss_rate)
    }
  }
  
  data0 <- data0[miss==0,]
  
  return(data0)
}
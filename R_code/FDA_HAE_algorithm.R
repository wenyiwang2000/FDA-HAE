FDA_HAE_algorithm <- function(R=2,L=2,nm=7,p=3,iters=1000,data1,tols=1e-200,
                                   tol0 = 1e-3, tol1 = 5e-3, tol2 = 1e-3, tol3 = 1e-7,Miniter = 20){
  tick  <- proc.time()
  n <- length(unique(data1$subj))
  N_obs <- length(data1$y)
  ns <- length(unique(data1$group))
  np <- length(unique(data1$position))
  K <- nm+3
  
  mt1 <- mgcv::gam(y~s(time,bs='ps', m=c(2,2), k=nm+3), drop.intercept=F,method="REML",data=data1, family=poisson)
  
  #BG <- as.matrix(spline.des(knots=construct.knots(seq(0,1,0.00001),7,knots.option="equally-spaced",p), 
  #                           x=seq(0,1,0.00001), outer.ok=TRUE, sparse=TRUE)$design)
  sm <- smoothCon(s(time, bs = "ps", m=c(2,2), k=nm+3), data = data1)[[1]]
  data1t <- data.frame(subj=rep(1,100001),time=seq(0,1,0.00001),l=c(0.9999999,rep(1,100000)),y=rep(1,100001))
  mt2 <- mgcv::gam(y~s(time,bs='ps',by=l, m=c(2,2), k=nm+3), drop.intercept=T,method="REML",data=data1t, family=poisson)
  BG <- predict.gam(mt2,type="lpmatrix")
  #BG <- as.matrix(spline.des(knots=sm$knots, 
  #                           x=seq(0,1,0.00001), outer.ok=TRUE, sparse=TRUE)$design)
  
  G <- crossprod(BG) / nrow(BG)
  eig_G <- eigen(G, symmetric = T)
  G_half <- eig_G$vectors %*% diag(sqrt(eig_G$values)) %*% t(eig_G$vectors)
  G_invhalf <- solve(G_half)
  
  
  #tB <- as.matrix(spline.des(knots=construct.knots(seq(0,1,0.01),7,knots.option="equally-spaced",p), 
  #                           x=seq(0,1,0.01), outer.ok=TRUE, sparse=TRUE)$design)
  #tB <- as.matrix(spline.des(knots=sm$knots, 
  #                           x=seq(0,1,0.01), outer.ok=TRUE, sparse=TRUE)$design)
  
  data1tt <- data.frame(subj=rep(1,101),time=seq(0,1,0.01),l=c(0.9999999,rep(1,100)),y=rep(1,101))
  mt3 <- mgcv::gam(y~s(time,bs='ps',by=l, m=c(2,2), k=nm+3), drop.intercept=T,method="REML",data=data1tt, family=poisson)
  tB <- predict.gam(mt3,type="lpmatrix")
  tB_theta <- tB %*% G_invhalf
  #Basis matrix
  
  sm <- smoothCon(s(time, bs = "ps", m=c(2,2), k=nm+3), data = data1)[[1]]
  tBmatrix <- sm$X   #for alpha
  tBmatrix_theta <- tBmatrix %*% G_invhalf   #for theta
  
  #indicator
  for(g in 1:ns){
    data1[,paste0("g",g)] <- as.numeric(data1$group==as.character(g))
  }
  for(p in 1:np){
    data1[,paste0("p",p)] <- as.numeric(data1$position==as.character(p))
  }
  for(g in 1:ns){
    for(p in 1:np){
      data1[,paste0("gp",g,p)] <- as.numeric(data1$group==as.character(g) & data1$position==as.character(p))
    }
  }
  
  #matrix indicator
  X <- data1[,6+ns+np+(1:(ns*np))]
  ncX <- dim(X)[2]
  
  
  ###################### initial value ########################
  eta0 <- rep(NA,N_obs)
  A <- array(0,c(R,R,R))
  U1 <- matrix(0,nrow=ns,ncol=R)
  U2 <- matrix(0,nrow=np,ncol=R)
  U3 <- matrix(0,nrow=K,ncol=R)
  fmean <- matrix(0,nrow=N_obs,ncol=R)
  
  M <- matrix(0,nrow=K,ncol=ncX)
  data_i <- data1
  data_i$y <- floor(data_i$y/max(data_i$Nij,1))
  data_i$y <- log(data_i$y+1)
  
  for(i in 1:ncX){
    index <- (X[,i]==1)
    data_initial_temp <- data_i[index,]
    m_initial_temp <- gam(y~s(time,bs='ps', m=c(2,2), k=nm+3), drop.intercept=F,method="REML",data=data_initial_temp)
    Kt <- m_initial_temp$smooth[[1]]$bs.dim 
    m_hat <- c(0, m_initial_temp$coefficients[2:Kt])
    qrc <- attr(m_initial_temp$smooth[[1]],"qr")
    m_hat <- m_initial_temp$coefficients[1] + qr.qty(qrc,m_hat)
    m_hat <- G_half %*% m_hat
    M[,i] <- m_hat
    eta0[index] <-tBmatrix_theta[index,] %*% m_hat
  }
  
  #tucker decomposition: HOOI
  Mt <- array(0,c(ns,np,K))
  for(s in 1:ns){
    for(p in 1:np){
      for(k in 1:K){
        Mt[s,p,k] <- M[k,(s-1)*p+p]
      }
    }
  }
  Mt_tnsr <- as.tensor(Mt)
  tucker_decomp0 <- tucker(Mt_tnsr, ranks = c(R, R, R))
  A <- tucker_decomp0$Z
  A <- A@data
  U1 <- tucker_decomp0$U[[1]]
  U2 <- tucker_decomp0$U[[2]]
  U3 <- tucker_decomp0$U[[3]]
  fmean <- tBmatrix_theta %*% U3
  Mt_tnsrnew <- ttl(tucker_decomp0$Z, tucker_decomp0$U, 1:3)
  Mtnew <- Mt_tnsrnew@data
  #fnorm(Mt_tnsr-Mt_tnsrnew) / fnorm(Mt_tnsr)
  #tucker_decomp0$norm_percent
  
  #eta0
  for(s in 1:ns){
    for(p in 1:np){
      index <- (data1$group==s & data1$position==p)
      eta0[index] <- tBmatrix_theta[index,] %*% Mtnew[s,p,]
    }
  }
  
  #################################################svd_0 <- svd(M)
  #################################################U <- matrix(svd_0$u[,1:L1],ncol=L1)
  #################################################Gamma <- diag(svd_0$d[1:L1],nrow=L1,ncol=L1) %*% t(svd_0$v[,1:L1])
  #################################################psi <- tBmatrix_theta %*% U
  
  
  
  #lambda
  tdata <- as.data.frame(cbind("argvals"=data1$time,"subj"=data1$subj,"y"=data_i$y-eta0))
  fit_face <- face.sparse(tdata,argvals.new=seq(0,1,0.01),calculate.scores=TRUE)
  L_face <- length(fit_face$eigenvalues)
  if(L_face >= L){
    lambda <- fit_face$eigenvalues[1:L]
  }else{
    lambda <- c(fit_face$eigenvalues,rep(0,L-L_face))
  }
  
  #theta
  theta <- c()
  if(L_face >= L){
    theta <- as.matrix(fit_face$U[,1:L],ncol=L)
  }else{
    theta <- cbind(fit_face$U,matrix(0,nrow=nm+3,ncol=L-L_face))
  }
  
  #orthonormal
  temp <- tcrossprod(theta %*% diag(lambda,nrow=L,ncol=L), theta)
  temp <- 0.5*(temp+t(temp))
  eig1 <- eigen(temp)
  lambda <- eig1$values[1:L]
  for(k in 1:L){
    if(crossprod(theta[,k], eig1$vectors[,k]) < 0){
      theta[,k] <- - eig1$vectors[,k]
    }else{
      theta[,k] <- eig1$vectors[,k] 
    }
  } 
  phi <- tBmatrix_theta %*% theta
  Lambda <- diag(lambda,L,L)
  
  #xi (score/zeta)
  xi <- matrix(0,n,L)    #n*L
  for(i in 1:L){
    xi[,i] <- rnorm(n,0,sqrt(lambda[i]))
  }
  
  
  keep_A <- array(0,c(iters,R,R,R))
  keep_U1 <- array(0,c(iters,ns,R))
  keep_U2 <- array(0,c(iters,np,R))
  keep_U3 <- array(0,c(iters,K,R))
  keep_theta <- array(0,c(iters,nm+3,L))
  keep_lambda <- matrix(0,iters,L)
  loss <- matrix(0,iters,3)
  ll <- rep(0,iters)
  keep_xi <- array(0,c(iters,n,L))
  
  A_rel_d <- rep(0,iters)
  U1_rel_d <- rep(0,iters)
  U2_rel_d <- rep(0,iters)
  U3_rel_d <- rep(0,iters)
  theta_rel_d <- rep(0,iters)
  lambda_rel_d <- rep(0,iters)
  
  A_abs_d <- rep(0,iters)
  U1_abs_d <- rep(0,iters)
  U2_abs_d <- rep(0,iters)
  U3_abs_d <- rep(0,iters)
  theta_abs_d <- rep(0,iters)
  lambda_abs_d <- rep(0,iters)
  
  keep_A[1,,,] <- A
  keep_U1[1,,] <- U1
  keep_U2[1,,] <- U2
  keep_U3[1,,] <- U3
  keep_theta[1,,] <- theta
  keep_lambda[1,] <- lambda
  keep_xi[1,,] <- xi
  
  
  for(iter in 2:iters){
    print(iter)
    
    #calculate eta_i(Tij), Psi, Sigma_i, R_i 
    eta_list <- list()
    Sigma_list <- list()
    Sigma_list_inv <- list()
    tPSP_list <- list()
    LtPSP_list <- list()
    tPRP_list <- list()
    for(i in 1:n){
      index <- data1$subj==as.character(i)
      phii <- matrix(phi[index,],ncol=L)
      eta_list[[i]] <-  log(data1$Nij)[index]+eta0[index] + phii %*% matrix(xi[i,],ncol=1)
      eta_list[[i]] <-  sapply(eta_list[[i]],function(x){min(x,4*log(10))})
      eta_list[[i]] <-  sapply(eta_list[[i]],function(x){max(x,-1*log(10))})
      temp <- as.vector(exp(-eta_list[[i]]))
      #temp <- sapply(temp,function(x){max(x,1e-4)})
      #temp <- sapply(temp,function(x){min(x,1e2)})
      Sigma_list[[i]] <- diag(temp,nrow=length(as.vector(eta_list[[i]])),ncol=length(as.vector(eta_list[[i]])))
      Sigma_list_inv[[i]] <- diag(1/temp,nrow=length(as.vector(eta_list[[i]])),ncol=length(as.vector(eta_list[[i]])))
      tPSP_list[[i]] <- t(phii) %*% Sigma_list_inv[[i]] %*% phii
      LtPSP_list[[i]] <- solve(diag(1/lambda,L,L)+tPSP_list[[i]],tol=tols)
      tPRP_list[[i]] <- tPSP_list[[i]]-tPSP_list[[i]] %*% LtPSP_list[[i]] %*% tPSP_list[[i]]
    }
    print("eta")
    print(range(eta_list))
    #print(range(lapply(Sigma_list,function(x){diag(x)})))
    #print(range(lapply(Sigma_list_inv,function(x){diag(x)})))
    
    #calculate Y_star
    Y_star <- c()
    for(i in 1:n){
      index <- data1$subj==as.character(i)
      phii <- matrix(phi[index,],ncol=L)
      term1 <- diag(Sigma_list[[i]])*(data1$y[index]-as.vector(diag(Sigma_list_inv[[i]])))
      term2 <- as.vector(eta_list[[i]])
      varpart <- phii %*% Lambda %*% t(phii) - phii %*% Lambda %*% tPRP_list[[i]] %*% Lambda %*% t(phii)
      #+phii %*% Lambda %*% tPRP_list[[i]] %*% solve(Reduce("+",tPRP_list),tol=tols) %*% tPRP_list[[i]]%*% Lambda %*% t(phii)
      term3 <- -(1/2)*diag(varpart)
      Y_star <- c(Y_star, term1+term2+term3-log(data1$Nij)[index])
    }
    
    #calculate conditional expectation and second moment of xi
    #expectation
    Exi_list <- list()
    Covxi_list <- list()
    Exi2_list <- list()
    for(i in 1:n){
      index <- data1$subj==as.character(i)
      phii <- matrix(phi[index,],ncol=L)
      Exi_list[[i]] <-  LtPSP_list[[i]] %*% t(phii) %*% Sigma_list_inv[[i]] %*% (Y_star[index]-eta0[index])
      Covxi_list[[i]] <- Lambda - Lambda %*% tPRP_list[[i]] %*% Lambda
      Exi2_list[[i]] <- Covxi_list[[i]] + Exi_list[[i]] %*% t(Exi_list[[i]])
    }
    
    
    ###################### update parameter ###################
    
    
    #fix_coef
    z0 <-c()
    Omega1 <- c()
    for(i in 1:n){
      index <- data1$subj==as.character(i) 
      Omega1 <- c(Omega1, as.vector(diag(Sigma_list_inv[[i]])))
      z0 <- c(z0,Y_star[index]-matrix(tBmatrix_theta[index,],ncol=nm+3) %*% theta %*% Exi_list[[i]])
    }
    
    eta0_t <- eta0
    
    # U3(fmean)
    df_U3 <- rep(0,R)
    for(jj in 1:10){
      U3_pre <- U3
      byterm <- matrix(0,nrow=N_obs,ncol=R)
      for(i in 1:n){
        index <- data1$subj==as.character(i) 
        ss <- data1$group[index][1]
        pp <- data1$position[index][1]
        for(r3 in 1:R){
          termt <- 0
          for(r1 in 1:R){
            for(r2 in 1:R){
              termt <- termt+A[r1,r2,r3]*U1[ss,r1]*U2[pp,r2]
            }
          }
          byterm[index,r3] <- rep(termt,sum(index))
        }
      }
      
      for(r3 in 1:R){
        term3 <- rep(0,N_obs)
        for(r_prime in 1:R){
          if(r_prime !=r3){
            term3 <- term3 + byterm[,r_prime]*(tBmatrix_theta %*% U3[,r_prime])
          }
        }
        zn1 <- z0 - term3
        newdata <- data1
        newdata$zn1 <- zn1
        newdata$Omega1 <- Omega1
        newdata$byterm <- byterm[,r3]
        m1 <- mgcv::gam(zn1~s(time, by=byterm ,bs='ps', m=c(2,2), k=nm+3), drop.intercept = T,method="REML",data=newdata, weights=Omega1)
        if(length(m1$coefficients)==nm+2){
          U3[,r3] <- G_half %*% c(m1$coefficients,0)
        }else{
          U3[,r3] <- G_half %*% m1$coefficients
        }
        #U3[,r3] <- G_half %*% m1$coefficients
        df_U3[r3] <- pen.edf(m1)
      }
      if(max(abs(U3_pre - U3)) < 1e-6) break
    }
    fmean <- tBmatrix_theta %*% U3
    
    
    #U1
    for(s in 1:ns){
      newdata1 <- data1[data1$group==s,]
      newfmean1 <- matrix(fmean[data1$group==s,],ncol=R)
      byterm <- matrix(0,nrow=dim(newdata1)[1],ncol=R)
      for(i in unique(newdata1$subj)){
        index <- newdata1$subj==as.character(i) 
        pp <- newdata1$position[index][1]
        for(r1 in 1:R){
          termt <- rep(0,sum(index))
          for(r2 in 1:R){
            for(r3 in 1:R){
              termt <- termt+A[r1,r2,r3]*U2[pp,r2]*newfmean1[index,r3]
            }
          }
          byterm[index,r1] <- termt
        }
      }
      model_mean <- lm(z0[data1$group==s] ~ byterm-1,weights=Omega1[data1$group==s])
      U1[s,] <- model_mean$coefficients
    }
    
    #U2
    for(p in 1:np){
      newdata1 <- data1[data1$position==p,]
      newfmean1 <- matrix(fmean[data1$position==p,],ncol=R)
      byterm <- matrix(0,nrow=dim(newdata1)[1],ncol=R)
      for(i in unique(newdata1$subj)){
        index <- newdata1$subj==as.character(i) 
        ss <- newdata1$group[index][1]
        for(r2 in 1:R){
          termt <- rep(0,sum(index))
          for(r1 in 1:R){
            for(r3 in 1:R){
              termt <- termt+A[r1,r2,r3]*U1[ss,r1]*newfmean1[index,r3]
            }
          }
          byterm[index,r2] <- termt
        }
      }
      model_mean <- lm(z0[data1$position==p] ~ byterm-1,weights=Omega1[data1$position==p])
      U2[p,] <- model_mean$coefficients
    }
    
    #A
    byterm <- matrix(0,nrow = N_obs,ncol=R*R*R)
    for(i in 1:n){
      index <- data1$subj==as.character(i) 
      ss <- data1$group[index][1]
      pp <- data1$position[index][1]
      for(r1 in 1:R){
        for(r2 in 1:R){
          for(r3 in 1:R){
            byterm[index,(r1-1)*R*R+(r2-1)*R+r3] <- U1[ss,r1]*U2[pp,r2]*fmean[index,r3]
          }
        }
      }
    }
    model_mean <- lm(z0 ~ byterm-1,weights=Omega1)
    model_mean_coef <- model_mean$coefficients
    for(r1 in 1:R){
      for(r2 in 1:R){
        for(r3 in 1:R){
          A[r1,r2,r3] <- model_mean_coef[(r1-1)*R*R+(r2-1)*R+r3]
        }
      }
    }
    
    
    #tucker decomposition: HOOI
    Mt1_tnsr <- ttl(as.tensor(A), list(U1=U1,U2=U2,U3=U3), 1:3)
    tucker_decomp1 <- tucker(Mt1_tnsr, ranks = c(R, R, R))
    A <- tucker_decomp1$Z
    A <- A@data
    U1 <- tucker_decomp1$U[[1]]
    U2 <- tucker_decomp1$U[[2]]
    U3 <- tucker_decomp1$U[[3]]
    fmean <- tBmatrix_theta %*% U3
    Mt1_tnsrnew <- ttl(tucker_decomp1$Z, tucker_decomp1$U, 1:3)
    Mt1new <- Mt1_tnsrnew@data
    #fnorm(Mt1_tnsr-Mt1_tnsrnew) / fnorm(Mt1_tnsr)
    #tucker_decomp1$norm_percent
    
    #eta0
    for(s in 1:ns){
      for(p in 1:np){
        index <- (data1$group==s & data1$position==p)
        eta0[index] <- tBmatrix_theta[index,] %*% Mt1new[s,p,]
      }
    }
    
    
    
    #theta
    df_theta <- rep(0,L)
    for(jj in 1:10){
      theta_pre <- theta
      for(l in 1:L){
        zn1 <- c()
        byterm <- c()
        for(i in 1:n){
          index <- data1$subj==as.character(i) 
          term1 <- (Exi2_list[[i]][l,l])^(-1/2)
          term2 <- (Y_star[index]-eta0[index])*Exi_list[[i]][l,]
          term3 <- 0
          for(l_prime in 1:L){
            if(l_prime != l){
              term3 <- term3 + Exi2_list[[i]][l_prime,l]*matrix(tBmatrix_theta[index,],ncol=nm+3) %*% theta[,l_prime]
            }
          }
          zn1 <- c(zn1, term1*(term2-term3) )
          byterm <- c(byterm, rep((Exi2_list[[i]][l,l])^(1/2),sum(index)))
        }
        newdata <- data1
        newdata$zn1 <- zn1
        newdata$Omega1 <- Omega1
        newdata$byterm <- byterm
        m1 <- mgcv::gam(zn1~s(time, by=byterm ,bs='ps', m=c(2,2), k=nm+3), drop.intercept = T,method="REML",data=newdata, weights=Omega1)
        if(length(m1$coefficients)==nm+2){
          theta[,l] <- G_half %*% c(m1$coefficients,0)
        }else{
          theta[,l] <- G_half %*% m1$coefficients
        }
        #theta[,l] <- G_half %*% m1$coefficients
        df_theta[l] <- pen.edf(m1)
      }
      if(max(abs(theta_pre - theta)) < 1e-6) break
    }
    
    
    
    #lambda
    sum <- rep(0,L)
    for(i in 1:n){
      index <- data1$subj==as.character(i) 
      sum <- sum + diag(Exi2_list[[i]])
    }
    lambda <- sum/n
    #lambda <- sapply(lambda,function(x){min(x,10)})
    
    #orthonormal
    temp <- tcrossprod(theta %*% diag(lambda,L,L), theta)
    temp <- 0.5*(temp+t(temp))
    eig1 <- eigen(temp)
    lambda <- eig1$values[1:L]
    for(k in 1:L){
      if(crossprod(theta[,k], eig1$vectors[,k]) < 0){
        theta[,k] <- - eig1$vectors[,k]
      }else{
        theta[,k] <- eig1$vectors[,k] 
      }
    } 
    phi <- tBmatrix_theta %*% theta
    Lambda <- diag(lambda,L,L)
    print("lambda")
    print(lambda)
    
    
    
    #calculate eta_i(Tij), Psi, Sigma_i, R_i 
    eta_list <- list()
    Sigma_list <- list()
    Sigma_list_inv <- list()
    tPSP_list <- list()
    LtPSP_list <- list()
    tPRP_list <- list()
    for(i in 1:n){
      index <- data1$subj==as.character(i)
      phii <- matrix(phi[index,],ncol=L)
      eta_list[[i]] <-  log(data1$Nij)[index]+eta0[index] + phii %*% matrix(xi[i,],ncol=1)
      eta_list[[i]] <-  sapply(eta_list[[i]],function(x){min(x,4*log(10))})
      eta_list[[i]] <-  sapply(eta_list[[i]],function(x){max(x,-1*log(10))})
      temp <- as.vector(exp(-eta_list[[i]]))
      #temp <- sapply(temp,function(x){max(x,1e-4)})
      #temp <- sapply(temp,function(x){min(x,1e2)})
      Sigma_list[[i]] <- diag(temp,nrow=length(as.vector(eta_list[[i]])),ncol=length(as.vector(eta_list[[i]])))
      Sigma_list_inv[[i]] <- diag(1/temp,nrow=length(as.vector(eta_list[[i]])),ncol=length(as.vector(eta_list[[i]])))
      tPSP_list[[i]] <- t(phii) %*% Sigma_list_inv[[i]] %*% phii
      LtPSP_list[[i]] <- solve(diag(1/lambda,L,L)+tPSP_list[[i]],tol=tols)
      tPRP_list[[i]] <- tPSP_list[[i]]-tPSP_list[[i]] %*% LtPSP_list[[i]] %*% tPSP_list[[i]]
    }
    
    #update xi
    xi <- matrix(0,n,L)    #n*L
    for(i in 1:n){
      index <- data1$subj==as.character(i)
      phii <- matrix(phi[index,],ncol=L)
      xi[i,] <-  as.vector(LtPSP_list[[i]] %*% t(phii) %*% Sigma_list_inv[[i]] %*% (Y_star[index]-eta0[index]))
    }
    
    keep_A[iter,,,] <- A
    keep_U1[iter,,] <- U1
    keep_U2[iter,,] <- U2
    keep_U3[iter,,] <- U3
    keep_theta[iter,,] <- theta
    keep_lambda[iter,] <- lambda
    keep_xi[iter,,] <- xi
    
    if(iter >1){
      Out <- c()
      for(i in 1:n){
        if(iter < 20){
          N=2000
        }else{
          N=10000
        }
        index <- data1$subj==as.character(i) 
        mi <- length(data1$y[index])
        xii <- rmvn(N, rep(0,L), diag(lambda,L,L)) #N*L
        vari <- xii %*% t(matrix(phi[index,],ncol=L)) #N*mi
        #etai <- matrix(rep(eta0[index], N), nrow=N, byrow=T) + vari + log(data1$Nij[index])  #N*mi
        etai <- matrix(rep(eta0[index], N), nrow=N, byrow=T) + vari + matrix(rep(log(data1$Nij[index]),N),nrow=N,byrow=TRUE) #N*mi
        mui <- exp(etai) #N*mi
        #fi <- matrix(mapply(dpois,matrix(rep(data1$y[index], N), nrow=N, byrow=T),mui),nrow=N,byrow=F)#N*mi
        ai <- -rowSums(mui) + etai %*% matrix(data1$y[index], ncol=1, byrow=T)  #N*1
        max_ai <- max(ai)
        log_fi <- max_ai + log(sum(exp(ai-max_ai)))
        #fprodi <- apply(fi, 1, prod)#N*1
        #out <- MC_Expec_l(eta0[index],Phi[index,],lambda,2000,data1$y[index])
        Out <- cbind(Out, log_fi)
      }
      #flike <- mean(apply(Out, 1, prod))
      #log_like <- log(flike)
      ll[iter] <- sum(Out)
      loss[iter,1] <- abs(ll[iter]-ll[iter-1]) / (abs(ll[iter-1]) + tol0)
      
      A_rel <- max(abs(keep_A[iter,,,]-keep_A[iter-1,,,])/(abs(keep_A[iter-1,,,])+tol0),na.rm = T)
      U1_rel <- max(abs(keep_U1[iter,,] %*% t(keep_U1[iter,,])-keep_U1[iter-1,,] %*% t(keep_U1[iter-1,,]))/(abs(keep_U1[iter-1,,] %*% t(keep_U1[iter-1,,]))+tol0),na.rm = T)
      U2_rel <- max(abs(keep_U2[iter,,] %*% t(keep_U2[iter,,])-keep_U2[iter-1,,] %*% t(keep_U2[iter-1,,]))/(abs(keep_U2[iter-1,,] %*% t(keep_U2[iter-1,,]))+tol0),na.rm = T)
      U3_rel <- max(abs(keep_U3[iter,,] %*% t(keep_U3[iter,,])-keep_U3[iter-1,,] %*% t(keep_U3[iter-1,,]))/(abs(keep_U3[iter-1,,] %*% t(keep_U3[iter-1,,]))+tol0),na.rm = T)
      theta_rel <- max(abs(keep_theta[iter,,]-keep_theta[iter-1,,])/(abs(keep_theta[iter-1,,])+tol0),na.rm = T)
      lambda_rel <- max(abs(keep_lambda[iter,]-keep_lambda[iter-1,])/(abs(keep_lambda[iter-1,])+tol0),na.rm = T)
      max_rel <- max(A_rel,U1_rel,U2_rel,U3_rel,theta_rel,lambda_rel,na.rm = T)
      
      A_abs <- max(abs(keep_A[iter,,,]-keep_A[iter-1,,,]),na.rm = T)
      U1_abs <- max(abs(keep_U1[iter,,] %*% t(keep_U1[iter,,])-keep_U1[iter-1,,] %*% t(keep_U1[iter-1,,])),na.rm = T)
      U2_abs <- max(abs(keep_U2[iter,,] %*% t(keep_U2[iter,,])-keep_U2[iter-1,,] %*% t(keep_U2[iter-1,,])),na.rm = T)
      U3_abs <- max(abs(keep_U3[iter,,] %*% t(keep_U3[iter,,])-keep_U3[iter-1,,] %*% t(keep_U3[iter-1,,])),na.rm = T)
      theta_abs <- max(abs(keep_theta[iter,,]-keep_theta[iter-1,,]),na.rm = T)
      lambda_abs <- max(abs(keep_lambda[iter,]-keep_lambda[iter-1,]),na.rm = T)
      max_abs <- max(A_abs,U1_abs,U2_abs,U3_abs,theta_abs,lambda_abs,na.rm = T)
      loss[iter,2] <- max_rel
      loss[iter,3] <- max_abs
      
      A_rel_d[iter] <- A_rel
      U1_rel_d[iter] <- U1_rel
      U2_rel_d[iter] <- U2_rel
      U3_rel_d[iter] <- U3_rel
      theta_rel_d[iter] <- theta_rel
      lambda_rel_d[iter] <- lambda_rel
      A_abs_d[iter] <- A_abs
      U1_abs_d[iter] <- U1_abs
      U2_abs_d[iter] <- U2_abs
      U3_abs_d[iter] <- U3_abs
      theta_abs_d[iter] <- theta_abs
      lambda_abs_d[iter] <- lambda_abs
    }
    
    #if(iter > 20){
    #  flag <- as.numeric(loss[iter,2] < tol1)+ as.numeric(loss[iter,3] < tol2)+ as.numeric(loss[iter,1] < tol3)
    #}else{
    #  flag <- as.numeric(loss[iter,2] < tol1)+ as.numeric(loss[iter,3] < tol2) 
    #}
    flag <- as.numeric(loss[iter,2] < tol1)+ as.numeric(loss[iter,3] < tol2) 
    print("loss")
    print(loss[iter,])
    
    if((iter > Miniter) && (flag >= 1)){
      break
    }
    
  } #end iter
  end <- iter
  
  #calculate marginal likelihood
  n_like <- 20
  log_like_v <- rep(0,n_like)
  for(i_like in 1:n_like){
    Out <- c()
    for(i in 1:n){
      N=10000
      index <- data1$subj==as.character(i) 
      mi <- length(data1$y[index])
      xii <- rmvn(N, rep(0,L), diag(lambda,L,L)) #N*L
      vari <- xii %*% t(matrix(phi[index,],ncol=L)) #N*mi
      #etai <- matrix(rep(eta0[index], N), nrow=N, byrow=T) + vari + log(data1$Nij[index]) #N*mi
      etai <- matrix(rep(eta0[index], N), nrow=N, byrow=T) + vari + matrix(rep(log(data1$Nij[index]),N),nrow=N,byrow=TRUE) #N*mi
      mui <- exp(etai) #N*mi
      #fi <- matrix(mapply(dpois,matrix(rep(data1$y[index], N), nrow=N, byrow=T),mui),nrow=N,byrow=F)#N*mi
      ai <- -rowSums(mui) + etai %*% matrix(data1$y[index], ncol=1, byrow=T)  #N*1
      max_ai <- max(ai)
      log_fi <- max_ai + log(sum(exp(ai-max_ai)))
      #fprodi <- apply(fi, 1, prod)#N*1
      #out <- MC_Expec_l(eta0[index],Phi[index,],lambda,2000,data1$y[index])
      Out <- cbind(Out, log_fi)
    }
    log_like_v[i_like] <- sum(Out)
  }
  log_like <- mean(log_like_v)
  
  #calculate dof
  df_A <- R*R*R
  df_U1 <- ns*R-R*(R+1)/2
  df_U2 <- np*R-R*(R+1)/2
  df <- df_A + df_U1 + df_U2 + sum(df_U3) - (R*(R+1))/2 + sum(df_theta) + L- (L*(L+1))/2 
  df_out <- c(df_U3,df_theta)
  
  AIC <- -2*log_like + 2*df
  BIC <- -2*log_like + log(N_obs)*df
  
  
  tock <- proc.time()
  out <- list("data"=data1,"keep_A"=keep_A ,"keep_U1"=keep_U1,"keep_U2"=keep_U2,"keep_U3"=keep_U3, "keep_theta"=keep_theta, "keep_lambda"=keep_lambda,
              "keep_xi"=keep_xi,"keep_loss"=loss,"ll"=ll,
              "tBmatrix"=tBmatrix, "tBmatrix_theta"=tBmatrix_theta,"converge"=iter, "tB"=tB, "tB_theta"=tB_theta,"time"=tock-tick,
              "AIC"=AIC, "BIC"=BIC,"N_obs"=N_obs,"df_out"=df_out,"n"=n, "df"=df, "log_like"=log_like,
              "xi"=xi,"A"=A,"U1"=U1,"U2"=U2,"U3"=U3, "theta"=theta,"lambda"=lambda,"eta0"=eta0,"phi"=phi,"fmean"=fmean,
              "X"=X,"R"=R,"L"=L,
              "A_rel_d"=A_rel_d,"U1_rel_d"=U1_rel_d,"U2_rel_d"=U2_rel_d,"U3_rel_d"=U3_rel_d,"theta_rel_d"=theta_rel_d, "lambda_rel_d"=lambda_rel_d,
              "A_abs_d"=A_abs_d,"U1_abs_d"=U1_abs_d,"U2_abs_d"=U2_abs_d,"U3_abs_d"=U3_abs_d,"theta_abs_d"=theta_abs_d, "lambda_abs_d"=lambda_abs_d)
  return(out)
}

# noisy sequence function
# internal to multiness_sim for gaussian model

# array of low rank matrices with common factors and
# symmetric iid normal noise

# rho controls the correlation between V and U
# there will also be correlation between say U_1 and U_2 induced by this
noisy_sequence <- function(n,m,d1,d2,sigma,rho,
                           dependence_type="all",hollow=identity,
                           gamma=1){
  # make this robust to scalar sigma
  if(length(sigma)==1){
    sigma_vec <- rep(sigma,m)
  }
  else{
    sigma_vec <- sigma
  }
  # generate common low rank structure
  if(d1 > 0){
    V <- matrix(stats::rnorm(n*d1,sd=gamma),n,d1)
  }
  else{
    V <- NULL
  }
  if(d2 > 0){
    if(dependence_type=="all"){
      # generate individual low rank structure
      W <- array(stats::rnorm(n*d2*m,sd=gamma),c(n,d2,m))
      # combine V and W to get correlated U
      U <- array(NA,dim(W))
      for(ii in 1:m){
        B <- t(rstiefel::rustiefel(d1,min(d1,d2)))
        w_scale <- rep(sqrt(1-rho^2),min(d1,d2))
        if(d2>d1){
          B <- rbind(B,matrix(0,d2-d1,d1))
          w_scale <- c(w_scale,rep(1,d2-d1))
        }
        U[,,ii] <- rho*(V %*% t(B)) + (W[,,ii] %*% diag(w_scale))
      }
    }
    if(dependence_type=="U_only"){ # introduce dependence through a single common factor in the U's
      U_all <- gamma*sqrt(n)*rstiefel::rustiefel(n,1+m*d2)
      U <- array(NA,c(n,d2,m))
      for(ii in 1:m){
        rho_bar <- sqrt(1-rho^2)
        col_ind <- 1+(ii-1)*d2+1:d2
        U[,,ii] <- rho*matrix(rep(U_all[,1],d2),n,d2) + sqrt(1-rho^2)*U_all[,col_ind]
      }
    }
  }
  else{
    U <- NULL
  }
  # allocate space
  P <- A <- array(NA,c(n,n,m))
  for(ii in 1:m){
    # calculate P and E
    P[,,ii] <- hollow(tcrossprod(cbind(V,U[,,ii])))
    E <- matrix(stats::rnorm(n^2,sd=sigma_vec[ii]),n,n)
    E[lower.tri(E)] <- t(E)[lower.tri(E)]
    A[,,ii] <- hollow(P[,,ii] + E)
  }
  return(list(V=V,U=U,P=P,A=A,d1=d1,d2=d2))
}

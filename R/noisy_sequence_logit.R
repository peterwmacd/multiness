# noisy_sequence_logit
# internal to multiness_sim in the logistic case

noisy_sequence_logit <- function(n,m,d1,d2,
                                 density_shift,gamma,rho,
                                 hollow=identity){
  # generate common low rank structure
  if(d1 > 0){
    V <- matrix(stats::rnorm(n*d1,sd=gamma),n,d1)
  }
  else{
    V <- NULL
  }
  if(d2 > 0){
    # generate individual low rank structure
    W <- array(stats::rnorm(n*d2*m,sd=gamma),c(n,d2,m))
    # combine V and W to get correlated U
    U <- array(NA,dim(W))
    for(ii in 1:m){
      if(d1 > 0){
        B <- t(rand_orth(d1,min(d1,d2)))
        w_scale <- rep(sqrt(1-rho^2),min(d1,d2))
        if(d2>d1){
          B <- rbind(B,matrix(0,d2-d1,d1))
          w_scale <- c(w_scale,rep(1,d2-d1))
        }
        U[,,ii] <- rho*(V %*% t(B)) + (W[,,ii] %*% diag(w_scale))
      }
      else{
        U[,,ii] <- W[,,ii]
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
    P.temp <- expit(tcrossprod(cbind(V,U[,,ii])) - density_shift)
    P[,,ii] <- hollow(P.temp)
    P.tri <- c(P.temp[lower.tri(P.temp,diag=T)])
    A.temp <- matrix(NA,n,n)
    A.temp[lower.tri(A.temp,diag=T)] <- as.numeric(stats::rbinom(n=length(P.tri),size=1,prob=P.tri))
    A.temp[upper.tri(A.temp)] <- t(A.temp)[upper.tri(A.temp)]
    A[,,ii] <- hollow(A.temp)
  }
  return(list(V=V,U=U,P=P,A=A,d1=d1,d2=d2))
}

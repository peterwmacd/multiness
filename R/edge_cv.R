# utility functions for edge cross-validation
# (1) generating edge missingness
# (2) evaluating holdout error

# generate a (symmetric) missingness pattern matrix
# (generate as a boolean array to save storage space)
# subset: matrix of same size with the subset of entries to
# remove
missing_mat_sym <- function(n,p,subset=matrix(TRUE,n,n)){
  # list of upper triangle indices
  ut <- which(subset&upper.tri(subset),arr.ind=TRUE)
  n_ut <- nrow(ut)
  # sample rows to remove
  ut_samp <- ut[sample(1:n_ut,floor(p*n_ut)),]
  samp <- rbind(ut_samp,ut_samp[,c(2,1)])
  # zero entries
  M <- matrix(TRUE,n,n)
  for(ii in 1:nrow(samp)){
    r <- samp[ii,]
    M[r[1],r[2]] <- FALSE
  }
  return(M)
}

# generate a (symmetric) missingness pattern array
# PWM: maybe make this visible to streamline layer holdout experiment
missing_array_sym <- function(n,m,p,subset=array(TRUE,c(n,n,m))){
  A <- array(NA,c(n,n,m))
  for(ii in 1:m){
    A[,,ii] <- missing_mat_sym(n,p,subset=subset[,,ii])
  }
  return(A)
}

# get the prediction error for a given fit, missingness
holdout_error <- function(A,fit,hollow,misspattern,
                          link=identity){
  m <- dim(A)[3]
  # initialize error vector
  obs <- NULL
  est <- NULL
  for(ii in 1:m){
    tri_ind <- upper.tri(A[,,ii],diag=!hollow)
    miss <- !c(misspattern[,,ii][tri_ind])
    obs <- c(obs,c(A[,,ii][tri_ind])[miss])
    P <- link(fit$F_hat + fit$G_hat[[ii]])
    est <- c(est,c(P[tri_ind])[miss])
  }
  return(sum((obs-est)^2))
}

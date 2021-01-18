# prox gd refit
#library(rARPACK)
#library(glmnet)
#source("aux_functions.R")

# refitting eigenvalues with glm as post-processing for 
# convex/soft-thresholding approach

# input: fit object from nuc-norm, glm family for glmnet
# output: new fit object, factor and matrix reps.

# coded for matrix-valued estimators version

prox_gd_refit <- function(A,fit,family="gaussian",
                          misspattern=NULL,hollow=F,
                          return_factors=F,ignore_neg=F){
  # get dimensions
  n <- dim(A)[1]
  if(hollow){
    nh <- n*(n-1)*0.5
  }
  else{
    nh <- n*(n-1)*0.5 + n
  }
  m <- dim(A)[3]
  if(ignore_neg){
    d1 <- sum(RSpectra::eigs_sym(fit$F_hat,fit$F_rank)$values>0)
    d2 <- sapply(1:m,function(ii){
      sum(RSpectra::eigs_sym(fit$G_hat[[ii]],fit$G_rank[ii])$values>0)
    })
  }
  else{
    d1 <- fit$F_rank
    d2 <- fit$G_rank
  }
  m2 <- sum(d2>0)
  # don't run if all ranks are zero
  if((d1+sum(d2))==0){
    if(return_factors){
      fit$V_hat <- matrix(0,n,1)
      fit$U_hat <- lapply(1:m,function(ii){matrix(0,n,1)})
    }
    return(fit)
  }
  # get response
  y <- c(apply(A,3,function(M){M[upper.tri(M,diag=!hollow)]}))
  # get predictors
  # V predictors
  if(d1>0){
    evec_v <- RSpectra::eigs_sym(fit$F_hat,d1)$vectors
    X_v <- NULL
    for(ii in 1:d1){
      temp <- outer(evec_v[,ii],evec_v[,ii])
      X_v <- cbind(X_v,rep(c(temp[upper.tri(temp,diag=!hollow)]),m))
    }
  }
  else{
    X_v <- NULL
  }
  # U predictors
  if(m2==0){
    all_X_u <- NULL
  }
  else{
    # initialize objects
    X_u <- evec_u <- list()
    #all_X_u <- Matrix(data=0,nrow=m*nh,ncol=sum(d2))
    length(X_u) <- m
    length(evec_u) <-  m2
    kk <- 0
    for(ii in 1:m){
      if(d2[ii]>0){
        kk <- kk+1
        evec_u[[kk]] <- RSpectra::eigs_sym(fit$G_hat[[ii]],d2[ii])$vectors
        for(jj in 1:d2[ii]){
          temp <- outer(evec_u[[kk]][,jj],evec_u[[kk]][,jj])
          X_u[[ii]] <- cbind(X_u[[ii]],c(temp[upper.tri(temp,diag=!hollow)]))
        }
      }
      else{
        X_u[[ii]] <- matrix(0,nh,1)
      }
    }
    all_X_u <- Matrix::bdiag(X_u)
    all_X_u <- all_X_u[,Matrix::colSums(all_X_u)!=0]
  }
  # complete predictors
  X <- cbind(X_v,all_X_u)
  # remove missing entries if needed
  if(!is.null(misspattern)){
    nonmissing <- c(apply(misspattern,3,function(M){M[upper.tri(M)]}))
    y <- y[nonmissing]
    X <- X[nonmissing,,drop=F]
  }
  # refit with glmnet, unpenalized
  # glmnet can handle sparse matrices as produced by bdiag
  if(dim(X)[2]>1){
    refit <- glmnet::glmnet(X,y,family=family,
                    alpha=0,lambda=0,
                    standardize=F,intercept=F)
  }
  else{
    temp <- stats::glm(y~X-1,family=family)
    refit <- list()
    refit$beta <- temp$coefficients
  }
  # return updated estimates
  # note: if new eigenvalues are not positive, can't return PSD rep.
  u_ind <- 1+d1+cumsum(c(0,d2))
  if(return_factors & all(refit$beta>0)){
    # store common factors
    if(d1>1){
      V_hathat <- evec_v %*% diag(sqrt(refit$beta[1:d1]))
    }
    if(d1==1){
      V_hathat <- matrix(evec_v*sqrt(refit$beta[1]),n,1)
    }
    if(d1==0){
      V_hathat <- matrix(0,n,1)
    }
    # store individual factors
    U_hathat <- list()
    kk <- 0
    for(ii in 1:m){
      if(d2[ii]>0){
        kk <- kk+1
        if(d2[ii]==1){
          U_hathat[[ii]] <- matrix(evec_u[[kk]] * sqrt(refit$beta[u_ind[ii]]),n,1)
        }
        else{
          U_hathat[[ii]] <- evec_u[[kk]] %*% diag(sqrt(refit$beta[u_ind[ii]:(u_ind[ii+1]-1)]))
        }
      }
      else{
        U_hathat[[ii]] <- matrix(0,n,1)
      }
    }
    # store common matrix
    F_hathat <- tcrossprod(V_hathat)
    # store individual matrices
    G_hathat <- lapply(U_hathat,tcrossprod)
  }
  else{
    # add a warning message if the factors can't be returned
    if(return_factors){
      warning("Can't factorize result, returning matrices instead")
    }
    # store common factors
    V_hathat <- NULL
    # store individual factors
    U_hathat <- NULL
    # store common matrix
    if(d1>0){
      F_hathat <- einfo_to_mat(list(vecs=evec_v,vals=refit$beta[1:d1]))
    }
    else{
      F_hathat <- matrix(0,n,n)
    }
    # store indiviual matrices
    G_hathat <- list()
    kk <- 0
    for(ii in 1:m){
      if(d2[ii]>0){
        kk <- kk+1
        G_hathat[[ii]] <- einfo_to_mat(list(vecs=evec_u[[kk]],vals=refit$beta[u_ind[ii]:(u_ind[ii+1]-1)]))
      }
      else{
        G_hathat[[ii]] <- matrix(0,n,n)
      }
    }
  }
  return(list(V_hat=V_hathat,U_hat=U_hathat,
              F_hat=F_hathat,G_hat=G_hathat,
              F_rank=d1,G_rank=d2,
              K=fit$K,convergence=fit$convergence))
}

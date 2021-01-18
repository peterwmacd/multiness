# utility/helper functions for proximal gradient descent algorithm

# update steps for F and G
# F update for step size is prescaled by m
V_update_f <- function(V_old,U_old,A_bar,link,
                       eta,lambda,
                       max_rank,
                       soft=T,eig_maxitr,
                       hollow=T,misspattern=NULL){
  m <- length(U_old)
  VV <- einfo_to_mat(V_old)
  if(hollow){ 
    thresh_mat <- VV + eta*hollowize(A_bar - (1/m)*sum_einfo_list(U_old,link[[1]],VV,misspattern))
  }
  else{
    thresh_mat <- VV + eta*(A_bar - (1/m)*sum_einfo_list(U_old,link[[1]],VV,misspattern))
  }
  V_new <- sv_thresh_f(thresh_mat,
                       thresh = (eta/m)*lambda,
                       max_rank = max_rank,
                       soft=soft,
                       eig_maxitr = eig_maxitr)
  return(V_new)
}

U_update_f <- function(U_old,V_old,A,link,
                       eta,lambda,
                       max_rank,soft=T,eig_maxitr,
                       hollow=T,misspattern=NULL){
  m <- dim(A)[3]
  VV <- einfo_to_mat(V_old)
  U_new <- lapply(1:m,function(ii){
    UU <- einfo_to_mat(U_old[[ii]])
    if(hollow){
      if(is.null(misspattern)){
        thresh_mat <- UU + eta*hollowize(A[,,ii] - link[[1]](VV+UU))
      }
      else{
        thresh_mat <- UU + eta*(hollowize(A[,,ii] - link[[1]](VV+UU))*misspattern[,,ii])
      }
    }
    else{
      if(is.null(misspattern)){
        thresh_mat <- UU + eta*(A[,,ii] - link[[1]](VV+UU))
      }
      else{
        thresh_mat <- UU + eta*((A[,,ii] - link[[1]](VV+UU))*misspattern[,,ii])
      }
    }
    return(sv_thresh_f(thresh_mat,
                       thresh = eta*lambda[ii],
                       max_rank=max_rank,
                       soft=soft,
                       eig_maxitr=eig_maxitr))
  }
  )
  return(U_new)
}

# evaluate the (least squares) objective function at a given F and G
objective <- function(F_0,G_0,A,la,lambda,link,hollow=F,misspattern=NULL,
                      max_rank=Inf,nuc_F=NULL,nuc_G=NULL){
  m <- dim(A)[3]
  if(link[[1]](1)==1){
    if(hollow){
      if(is.null(misspattern)){
        f1 <- .5*sum(sapply(1:m,function(ii){norm(hollowize(A[,,ii] - F_0 - G_0[[ii]]),type="F")^2}))
      }
      else{
        f1 <- .5*sum(sapply(1:m,function(ii){norm(hollowize(A[,,ii] - F_0 - G_0[[ii]])*misspattern[,,ii],type="F")^2}))
      }
    }
    else{
      if(is.null(misspattern)){
        f1 <- .5*sum(sapply(1:m,function(ii){norm(A[,,ii] - F_0 - G_0[[ii]],type="F")^2}))
      }
      else{
        f1 <- .5*sum(sapply(1:m,function(ii){norm((A[,,ii] - F_0 - G_0[[ii]])*misspattern[,,ii],type="F")^2}))
      }
    }
  }
  else{
    if(hollow){
      if(is.null(misspattern)){
        f1 <- sum(sapply(1:m,function(ii){sum(hollowize(-A[,,ii]*(F_0 + G_0[[ii]]) + log(1 + exp(F_0 + G_0[[ii]]))))}))
      }
      else{
        f1 <- sum(sapply(1:m,function(ii){sum(hollowize(-A[,,ii]*(F_0 + G_0[[ii]]) + log(1 + exp(F_0 + G_0[[ii]])))*misspattern[,,ii])}))
      }
    }
    else{
      if(is.null(misspattern)){
        f1 <- sum(sapply(1:m,function(ii){sum(-A[,,ii]*(F_0 + G_0[[ii]]) + log(1 + exp(F_0 + G_0[[ii]])))}))
      }
      else{
        f1 <- sum(sapply(1:m,function(ii){sum((-A[,,ii]*(F_0 + G_0[[ii]]) + log(1 + exp(F_0 + G_0[[ii]])))*misspattern[,,ii])}))
      }
    }
  }
  if(is.null(nuc_F)){
    f2 <- la*nuclear(F_0,max_rank)
  }
  else{
    f2 <- la*nuc_F
  }
  if(is.null(nuc_G)){
    f3 <- sum(lambda*sapply(G_0,nuclear,max_rank=max_rank))
  }
  else{
    f3 <- sum(lambda*nuc_G)
  }
  return(f1+f2+f3)
}

objective_factor <- function(V,U,A,la,lambda,link,hollow=F,misspattern=NULL,max_rank=Inf){
  # calculate fitted matrices
  F_0 <- einfo_to_mat(V)
  G_0 <- lapply(U,einfo_to_mat)
  # calculate nuclear norms
  nuc_F <- sum(abs(V$vals))
  nuc_G <- sapply(U,function(x){sum(abs(x$vals))})
  return(objective(F_0,G_0,A,la,lambda,link,hollow,misspattern,
                   max_rank,nuc_F,nuc_G))
}
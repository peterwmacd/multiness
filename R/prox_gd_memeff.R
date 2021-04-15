#setwd("~/nuclear_norm")
#source("aux_functions.R")
#source("data_functions.R")
#source("alt_svd.R")

# Method 1 - proximal gradient descent

prox_gd_memeff <- function(A,lambda,alpha,
                             link=list(identity,identity),
                             eta,init='svd',V_init=NULL,U_init=NULL,
                             pos=F,block=T,soft=T,hollow=T,misspattern=NULL,
                             eps=1e-4,max_rank=10,K_max=50,eig_maxitr=1000,init_rank=max_rank,
                             verbose=FALSE,check_obj=FALSE,eig_prec=1e-2){
  # dimensions
  n <- dim(A)[1]
  m <- dim(A)[3]
  # common penalty parameter
  la <- alpha*sqrt(mean(lambda^2))
  # precalculate mean of A
  if(is.null(misspattern)){
    A_bar <- apply(A,c(1,2),mean)
  }
  else{
    A_bar <- apply(A*misspattern,c(1,2),mean)
  }
  # initialization method
  # initialize with a low-rank approximation
  if(init=='svd'){
    V_old <- rank_thresh_f(link[[2]](rank_thresh(A_bar,max_rank,pos,eig_maxitr)),init_rank,pos,eig_maxitr)
    U_old <- lapply(1:m,function(ii){rank_thresh_f(link[[2]](rank_thresh(A[,,ii] - einfo_to_mat(V_old),max_rank,pos,eig_maxitr)),init_rank,pos,eig_maxitr)})
  }
  # initialize completely at random
  # if(init=='random'){
  #   F_old <- list(tcrossprod(matrix(rnorm(n*max_rank),n,max_rank)),max_rank)
  #   G_old <- lapply(1:m,function(ii){list(tcrossprod(matrix(rnorm(n*max_rank),n,max_rank)),max_rank)})
  # }
  # initialize at zero
  if(init=='zero'){
    V_old <- list(vals=rep(0,1),vecs=matrix(0,n,1))
    U_old <- lapply(1:m,function(ii){list(vals=rep(0,1),vecs=matrix(0,n,1))})
  }
  # initialize with previous fixed output
  if(init=="fix"){
    V_old <- V_init
    U_old <- U_init
  }
  # initialize variables
  K_count <- 0
  obj <- Inf
  eta_curr <- eta
  for(kk in 1:K_max){
    K_count <- K_count+1
    # update V
    V_new <- V_update_f(V_old,U_old,A_bar,link,
                      eta_curr,la,
                      max_rank,pos,soft,eig_maxitr,hollow,misspattern)
    # use V_old or V_new depending on block updating
    # update U
    if(block){
      U_new <- U_update_f(U_old,V_new,A,link,
                          eta_curr,lambda,
                          max_rank,pos,soft,eig_maxitr,hollow,misspattern)
    }
    else{
      U_new <- U_update_f(U_old,V_old,A,link,
                          eta_curr,lambda,
                          max_rank,pos,soft,eig_maxitr,hollow,misspattern)
    }
    # check convergence
    if(check_obj){
      obj_new <- objective_factor(V_new,U_new,
                                  A,la,lambda,link,
                                  hollow,misspattern,max_rank)
      keep <- (obj_new < obj)
      crit <- (obj - obj_new) / obj_new
      # also converge if eta gets very small
      conv <- (keep & (crit < eps))
      # update step sizes
      if(!keep){
        eta_curr <- .5*eta_curr
      }
    }
    else{
      keep <- TRUE
      obj_new <- Inf
      crit <- mean(abs(c(einfo_to_mat(V_new))-c(einfo_to_mat(V_old))))
      conv <- (crit < eps)
    }
    # PWM: make a more clear printout using cat(), only print if the step was taken
    if(verbose){
      if(keep){
        if(check_obj){
          cat("\nCompleted (attempted) iteration",K_count,
              "\nstep size eta =",eta_curr,
              "\nstopping criterion (relative decrease in obj.)",crit,
              "\n")
        }
        else{
          cat("\nCompleted (attempted) iteration",K_count,
              "\nstep size eta =",eta,
              "\nstopping criterion (av. entrywise change in est. of F)",crit,
              "\n")
        }
      }
    }
    # stop if converged, otherwise update estimates and step size
    if(conv){
      break
    }
    else{
      if(keep){
        # update
        V_old <- V_new
        U_old <- U_new
        obj <- obj_new
        eta_curr <- min(eta,2*eta_curr)
      }
    }
  }
  # set a convergence flag like in stats::optim
  # did the loop exit with convergence
  convergence <- as.integer(!conv)

  # trim very small eigenvalues (numerical precision)
  pos_ev <- c(min(V_new$vals) < -eig_prec,
              sapply(U_new,function(x){min(x$vals) < -eig_prec}))
  # check if ASE can be run on the results
  ase_ok <- !any(pos_ev)
  # *ase_ok is no longer used, from before ASE could
  # handle negative eigenvalues

  G_hat <- lapply(U_new,einfo_to_mat,eig_prec=eig_prec)
  return(list(F_hat=einfo_to_mat(V_new,eig_prec=eig_prec),G_hat=G_hat,
              F_rank=sum(abs(V_new$vals)>eig_prec),G_rank=sapply(U_new,function(x){sum(abs(x$vals)>eig_prec)}),
              K=K_count,convergence=convergence,ase_ok=ase_ok))
}

# functions to perform (hard/soft/rank constrained)
# singular value thresholding

# soft singular value thresholding for symmetric matrices
sv_thresh_f <- function(M,thresh,max_rank=dim(M)[1]-1,
                        soft=TRUE,pos=FALSE,
                        eig_maxitr=1000){
  # dimensions
  n <- dim(M)[1]
  # truncated eigendecomposition
  e <- RSpectra::eigs_sym(M,min(max_rank,n-1),
                which="LM",opts=list(maxitr=eig_maxitr))
  # thresholded rank
  if(pos){
    which_r <- which(e$values > thresh)
  }
  else{
    which_r <- which(abs(e$values) > thresh)
  }
  r <- length(which_r)
  if(r==0){
    return(list(vals=rep(0,1),vecs=matrix(0,n,1)))
  }
  else{
    if(soft){
      vals <- e$values[which_r] + (-sign(e$values[which_r]))*thresh
    }
    else{
      vals <- e$values[which_r]
    }
    vecs <- e$vectors[,which_r]
    return(list(vals=vals,vecs=vecs))
  }
}
# return a val/vec (einfo) object

# hard singular value thresholding for symmetric PSD matrices
# for initialization use the psd version (i.e. ignore negative eigenvalues)
rank_thresh <- function(M,max_rank,pos=FALSE,eig_maxitr=1000){
  # dimensions
  n <- dim(M)[1]
  # PSD version
  if(pos){
    estr <- "LA"
  }
  else{
    estr <- "LM"
  }
  # truncated eigendecomposition
  e <- RSpectra::eigs_sym(M,min(max_rank,n-1),
                which=estr,opts=list(maxitr=eig_maxitr))
  return(einfo_to_mat(list(vals=e$values,vecs=e$vectors)))
}

# rank_thresholding which returns an einfo object
rank_thresh_f <- function(M,max_rank,pos=FALSE,eig_maxitr=1000){
  # dimensions
  n <- dim(M)[1]
  # PSD version
  if(pos){
    estr <- "LA"
  }
  else{
    estr <- "LM"
  }
  # truncated eigendecomposition
  e <- RSpectra::eigs_sym(M,min(max_rank,n-1),
                which=estr,opts=list(maxitr=eig_maxitr))
  return(list(vals=e$values,vecs=e$vectors))
}

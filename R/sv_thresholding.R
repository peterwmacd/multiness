# functions to perform (hard/soft/rank constrained)
# singular value thresholding

# soft singular value thresholding for symmetric matrices
sv_thresh_f <- function(M,thresh,max_rank=dim(M)[1]-1,
                        soft=T,eig_maxitr=1000){
  # dimensions
  n <- dim(M)[1]
  # truncated eigendecomposition
  e <- RSpectra::eigs_sym(M,min(max_rank,n-1),
                which="LM",opts=list(maxitr=eig_maxitr))
  # thresholded rank
  r <- sum(abs(e$values) > thresh)
  if(r==0){
    return(list(vals=rep(0,1),vecs=matrix(0,n,1)))
  }
  else{
    if(soft){
      vals <- e$values[1:r] + (-sign(e$values[1:r]))*thresh
    }
    else{
      vals <- e$values[1:r]
    }
    vecs <- e$vectors[,1:r]
    return(list(vals=vals,vecs=vecs))
  }
}
# return a val/vec (einfo) object

# hard singular value thresholding for symmetric PSD matrices
# for initialization use the psd version (i.e. ignore negative eigenvalues)
rank_thresh <- function(M,max_rank,eig_maxitr=1000){
  # dimensions
  n <- dim(M)[1]
  # truncated eigendecomposition
  e <- RSpectra::eigs_sym(M,min(max_rank,n-1),
                which="LM",opts=list(maxitr=eig_maxitr))
  return(einfo_to_mat(list(vals=e$values,vecs=e$vectors)))
}

# rank_thresholding which returns an einfo object
rank_thresh_f <- function(M,max_rank,eig_maxitr=1000){
  # dimensions
  n <- dim(M)[1]
  # truncated eigendecomposition
  e <- RSpectra::eigs_sym(M,min(max_rank,n-1),
                which="LA",opts=list(maxitr=eig_maxitr))
  return(list(vals=e$values,vecs=e$vectors))
}

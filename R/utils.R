# utility functions

# logit and expit link functions

#' Inverse logistic link function
#'
#' \code{expit} applies the inverse logistic link function \eqn{f(x) = e^x / (1+e^x)}.
#'
#' @usage
#' expit(x)
#'
#' @param x A numeric vector.
#'
#' @export
expit <- function(x){
  1 - 1/(1+exp(x))
}

#' Logistic link function
#'
#' \code{logit} applies the logistic link function \eqn{f(x) = log(x / (1-x))}.
#'
#' @usage
#' logit(x,tol=1e-6)
#'
#' @param x A numeric vector with values in the interval \[0,1\].
#' @param tol A positive scalar which bounds the entries of x away from 0 and 1
#' for numerical stability.
#' Defaults to \code{tol=1e-6}
#'
#' @export
logit <- function(x,tol=1e-6){ # tolerance parameter to avoid NaN's
  y <- pmin(pmax(x,tol),1-tol)
  log(y/(1-y))
}

# restore a matrix from eig information
einfo_to_mat <- function(einfo,eig_prec=0){
  eval <- einfo$vals*(abs(einfo$vals) > eig_prec)
  return((einfo$vecs %*% (eval *t(einfo$vecs))))
}

# make an nxn matrix 'hollow' (zero diagonal)
hollowize <- function(M){
  n <- dim(M)[1]
  M*(1-diag(n))
}

# sum a list of matrices stored as eigenstructure (plus a common intercept/offset)
sum_einfo_list <- function(lof,link,offset,misspattern=NULL){
  m <- length(lof)
  n <- dim(offset)[1]
  if(m==1){
    return(einfo_to_mat(lof[[1]])*misspattern[,,1])
  }
  else{
    if(is.null(misspattern)){
      sum_einfo <- function(base,einfo){
        base + link(offset + einfo_to_mat(einfo))
      }
      return(Reduce(sum_einfo,lof,init=matrix(0,n,n)))
    }
    else{
      base <- matrix(0,n,n)
      for(ii in 1:m){
        base <- base + link(offset + einfo_to_mat(lof[[ii]]))*misspattern[,,ii]
      }
      return(base)
    }
  }
}

#' Adjacency Spectral Embedding (ASE)
#'
#' \code{ase} calculates the \eqn{d}-dimensional adjacency spectral embedding of a symmetric,
#' positive semi-definite \eqn{n \times n} matrix \eqn{M}.
#'
#' @usage
#' ase(M,d)
#'
#' @param M A symmetric and positive semi-definite matrix.
#' @param d A non-negative integer embedding dimension.
#'
#' @return An \eqn{n \times d} matrix \eqn{X} which factorizes
#' \eqn{M} as \eqn{M = XX^{\intercal}}. If \eqn{d=0}, \code{ase} returns
#' an \eqn{n \times 1} matrix of zeros.
#'
#' @export
ase <- function(M,d){
  if(d==0){
    X.hat <- matrix(0,nrow(M),1)
  }
  else{
    temp <- RSpectra::eigs_sym(M,d,which="LA")
    if(d>1){
      X.hat <- temp$vectors %*% diag(sqrt(temp$values))
    }
    else{
      X.hat <- temp$vectors*sqrt(temp$values)
    }
  }
  return(X.hat)
}

# evaluate the nuclear norm of a matrix (with a max rank)
nuclear <- function(M,max_rank){
  suppressWarnings(sum(abs(RSpectra::svds(M,min(max_rank,nrow(M)))$d)))
}

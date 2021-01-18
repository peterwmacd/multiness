#' Simulate from the MultiNeSS model
#'
#' \code{multiness_sim} simulates a realization of the Gaussian
#' or logistic MultiNeSS model with Gaussian latent positions.
#'
#' The common and individual latent positions, \eqn{V} and \eqn{U_k}
#' respectively, are generated as
#' Gaussian random variables with standard deviation \code{lp_opts$gamma}, and
#' dependence controlled by the optional
#' arguments \code{lp_opts$dependence_type} and \code{lp_opts$rho}.
#'
#' Under the Gaussian model, the \eqn{n \times n} adjacency matrix for layer \eqn{k=1,...,m}
#' has independent Gaussian entries with standard deviation \code{sigma} and
#' mean given by
#' \deqn{E(A_k) = VV^{\intercal} + U_kU_k^{\intercal}.}
#'
#' Under the logistic model, the \eqn{n \times n} adjacency matrix for layer \eqn{k=1,...,m}
#' has independent Bernoulli entries with mean given by
#' \deqn{E(A_k) = g(VV^{\intercal} + U_kU_k^{\intercal}),}
#' where \eqn{g} denotes the element-wise application of the inverse logistic
#' link (\code{\link{expit}}) function. Under both models, \code{self_loops} provides
#' an option to set the diagonal entries of the adjacency matrices to zero.
#'
#' @usage
#' multiness_sim(n,m,d1,d2,model,sigma,self_loops,lp_opts)
#'
#' @param n A positive integer, the number of nodes.
#' @param m A positive integer, the number of layers.
#' @param d1 A non-negative integer, the number of common latent dimensions.
#' @param d2 A non-negative integer, the number of individual latent dimensions.
#' @param model A string which provides choice of model,
#' either \code{'gaussian'} or \code{'logistic'}. Defaults to \code{'gaussian'}.
#' @param sigma A positive scalar or numeric vector of length \code{m},
#' the entry-wise standard deviation for the Gaussian noise for all layers
#' (if a scalar) or for each layer (if a vector). Ignored under the logistic
#' model. Defaults to \code{1}.
#' @param self_loops A Boolean, if \code{FALSE}, all diagonal entries are set
#' to zero. Defaults to \code{TRUE}.
#' @param lp_opts A list, containing additional optional arguments controlling the
#' properties of the latent positions:
#' \describe{
#'     \item{dependence_type}{A string, valid choices are \code{'all'} or
#'     \code{'U_only'} for the Gaussian model; \code{'all'} for the logistic
#'     model. If \code{'all'}, \eqn{V} and \eqn{U_k}; and \eqn{U_k} and \eqn{U_l}
#'     (for \eqn{k \neq l}) have expected canonical correlation approximately equal to
#'     |\eqn{rho}| (see \code{rho}).
#'     If \code{'U_only'}, \eqn{U_k} and \eqn{U_l} (for \eqn{k \neq l}) have expected
#'     canonical correlation approximately equal to |\eqn{rho}| (see \code{rho}).
#'     Defaults to \code{'all'}.}
#'     \item{gamma}{A positive scalar, the standard deviation of the entries of
#'     the latent position matrices \eqn{V} and \eqn{U_k}. Defaults to \code{1}.}
#'     \item{rho}{A positive scalar in the interval (-1,1), controls the expected canonical
#'     correlation between latent position matrices (see \code{dependence_type}).
#'     Defaults to \code{0}.}
#' }
#'
#'
#'
#' @return A list is returned with the realizations of the latent dimensions
#' and the multiplex network:
#' \item{A}{An array of dimension \eqn{n \times n \times m}, the realized
#' multiplex network.}
#' \item{V}{A matrix of dimension \eqn{n \times d1}, the realized common
#' latent positions. If \code{d1=0}, returns \code{NULL}.}
#' \item{U}{An array of dimension \eqn{n \times d2 \times m}, the realized
#' individual latent positions. If \code{d2=0}, returns \code{NULL}.}
#'
#' @examples
#' # gaussian model, uncorrelated latent positions
#' data1 <- multiness_sim(n=100,m=4,d1=2,d2=2,
#'                       model="gaussian")
#'
#' # logistic model, correlated latent positions
#' data2 <- multiness_sim(n=100,m=4,d1=2,d2=2,
#'                        model="logistic",
#'                        self_loops=FALSE,
#'                        lp_opts=list(dependence_type="all",rho=.3))
#'
#' @export
multiness_sim <- function(
  n,
  m,
  d1,
  d2, # just an integer for now
  model='gaussian', # either gaussian or logistic
  sigma=1, # either a scalar or a length m vector, unused in logit model
  self_loops=TRUE, # boolean (not a function!)
  lp_opts=list() # a list including gamma (entry sd), dependence info
){
  # PWM: fix this later
  # check all the dimensions are integer (can be coerced to integer?)
  # especially d2
  # if(!(is.integer(n) & is.integer(m) & is.integer(d1) & is.integer(d2))){
  #   stop('All dimensions must be integer')
  # }

  # set hollow function
  if(self_loops){
    hollow_fun <- identity
  }
  else{
    hollow_fun <- hollowize
  }

  # check form of sigma, convert to a vector if it isn't a scalar
  if(!(length(sigma)==m)){
    if(length(sigma)==1){
      sigma_vec <- rep(sigma,m)
    }
    else{
      stop('standard deviation must be either a scalar or a vector of length m')
    }
  }
  else{
    sigma_vec <- sigma
  }

  # set default latent position options
  # gamma
  if(is.null(lp_opts$gamma)){
    lp_opts$gamma <- 1
  }
  # rho
  if(is.null(lp_opts$rho)){
    lp_opts$rho <- 0
  }
  # dependence_type
  if(is.null(lp_opts$dependence_type)){
    # no dependence if the list entry is missing
    lp_opts$dependence_type <- 'all'
  }
  else{
    if(!(lp_opts$dependence_type=='all' || lp_opts$dependence_type=='U_only')){
      stop('Not a valid dependence type')
    }
  }

  # choose model and generate data from internal function
  if(model=='gaussian' || model=='Gaussian'){
    temp <- noisy_sequence(n=n,
                           m=m,
                           d1=d1,
                           d2=d2,
                           sigma=sigma_vec,
                           hollow=hollow_fun,
                           gamma=lp_opts$gamma,
                           rho=lp_opts$rho,
                           dependence_type=lp_opts$dependence_type)
  }
  else{
    if(model=='logistic'){
      if(lp_opts$dependence_type == 'U_only'){
        stop('Not a valid dependence type')
      }
      temp <- noisy_sequence_logit(n=n,
                                   m=m,
                                   d1=d1,
                                   d2=d2,
                                   gamma=lp_opts$gamma,
                                   rho=lp_opts$rho,
                                   hollow=hollow_fun)
    }
    else{
      stop('Not a valid model')
    }
  }

  # collect output and return
  return(list(A=temp$A,V=temp$V,U=temp$U))
}

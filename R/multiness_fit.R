#' Fit the MultiNeSS model
#'
#' \code{multiness_fit} fits the Gaussian or logistic MultiNeSS model
#' with various options for parameter tuning.
#'
#' A MultiNeSS model is fit to an \eqn{n \times n \times m} array \eqn{A} of
#' symmetric adjacency matrices on a common set of nodes. Fitting
#' proceeds by convex proximal gradient descent on the entries of
#' \eqn{F = VV^{\intercal}} and \eqn{G_k = U_kU_k^{\intercal}}, see
#' \href{https://arxiv.org/abs/2012.14409}{MacDonald et al., (2020)},
#' Section 3.2. Additional optional arguments for
#' the gradient descent routine can be provided in \code{optim_opts}.
#' \code{refit} provides an option
#' to perform an additional refitting step to debias the eigenvalues
#' of the estimates, see
#' \href{https://arxiv.org/abs/2012.14409}{MacDonald et al., (2020)}, Section 3.3.
#'
#' By default, \code{multiness_fit} will return estimates of the matrices
#' \eqn{F} and \eqn{G_k}. \code{optim_opts$return_posns} provides an option
#' to instead return estimates of latent positions \eqn{V} and \eqn{U_k}
#' based on the adjacency spectral embedding (if such a factorization exists).
#'
#' Tuning parameters \eqn{\alpha} and \eqn{\lambda_k} in the nuclear norm penalty
#' \deqn{\lambda ||F||_* + \sum_k \lambda \alpha_k ||G_k||_*}
#' are either set by the
#' user (\code{tuning='fixed'}), selected adaptively using the
#' \code{\link[denoiseR:estim_sigma]{denoiseR}}
#' package to estimate entry-wise variance (\code{tuning='adaptive'}), or
#' selected using edge cross-validation (\code{tuning='cv'}). For more details
#' see \href{https://arxiv.org/abs/2012.14409}{MacDonald et al., (2020)},
#' Section 3.4. Additional optional arguments for parameter tuning
#' can be provided in \code{tuning_opts}.
#'
#' @usage
#' multiness_fit(A,model,self_loops,refit,tuning,tuning_opts,optim_opts)
#'
#' @param A An \eqn{n \times n \times m} array containing edge entries for
#' an undirected multiplex network on \eqn{n} nodes and \eqn{m} layers.
#' @param model A string which provides choice of model,
#' either \code{'gaussian'} or \code{'logistic'}. Defaults to \code{'gaussian'}.
#' @param self_loops A Boolean, if \code{FALSE}, all diagonal entries are ignored in
#' optimization. Defaults to \code{TRUE}.
#' @param refit A Boolean, if \code{TRUE}, a refitting step is performed to
#' debias the eigenvalues of the estimates. Defaults to \code{TRUE}.
#' @param tuning A string which provides the tuning method, valid options are
#' \code{'fixed'}, \code{'adaptive'}, or \code{'cv'}. Defaults to \code{'adaptive'}.
#' @param tuning_opts A list, containing additional optional arguments controlling
#' parameter tuning. The arguments used depends on the choice of tuning method.
#' If \code{tuning='fixed'}, \code{multiness_fit} will utilize the following
#' arguments:
#' \describe{
#'     \item{lambda}{A positive scalar,
#'     the \eqn{\lambda} parameter in the nuclear norm penalty, see Details.
#'     Defaults to \code{2.309 * sqrt(n*m)}.}
#'     \item{alpha}{A positive scalar or numeric vector of length \code{m}, the parameters \eqn{\alpha_k} in the
#'     nuclear norm penalty, see Details. If a scalar is provided all \eqn{\lambda_k} parameters are set to that
#'     value. Defaults to \code{1/sqrt(m)}}
#' }
#' If \code{tuning='adaptive'}, \code{multiness_fit} will utilize the following
#' arguments:
#' \describe{
#'     \item{layer_wise}{A Boolean, if \code{TRUE}, the entry-wise variance
#'     is estimated individually for each layer. Otherwise the estimates are
#'     pooled. Defaults to \code{TRUE}.}
#'     \item{penalty_const}{A positive scalar \eqn{C} which scales the
#'     penalty parameters (see Details).
#'     Defaults to \code{2.309}.}
#'     \item{penalty_const_lambda}{A positive scalar \eqn{c} which scales only the \eqn{\lambda}
#'     penalty parameter (see Details).
#'     Defaults to \code{1}.}
#' }
#' If \code{tuning='cv'}, \code{multiness_fit} will utilize the following
#' arguments:
#' \describe{
#'    \item{layer_wise}{A Boolean, if \code{TRUE}, the entry-wise variance
#'     is estimated individually for each layer. Otherwise the estimates are
#'     pooled. Defaults to \code{TRUE}.}
#'     \item{N_cv}{A positive integer, the number of repetitions of edge
#'     cross-validation performed for each parameter setting. Defaults to \code{3}.}
#'     \item{p_cv}{A positive scalar in the interval (0,1), the proportion
#'     of edge entries held out in edge cross-validation. Defaults to \eqn{0.1}.}
#'     \item{penalty_const_lambda}{A positive scalar \eqn{c} which scales only the \eqn{\lambda}
#'     penalty parameter (see Details).
#'     Defaults to \code{1}.}
#'     \item{penalty_const_vec}{A numeric vector with positive entries, the candidate
#'     values of constant \eqn{C} to scale the penalty parameters (see Details).
#'     An optimal constant is chosen by edge cross-validation. Defaults to
#'     \code{c(1,1.5,...,3.5,4)}.}
#'     \item{refit_cv}{A Boolean, if \code{TRUE}, a refitting step is
#'     performed when fitting the model for edge cross-validation. Defaults
#'     to \code{TRUE}}
#'     \item{verbose_cv}{A Boolean, if \code{TRUE}, console output will
#'     provide updates on the progress of edge cross-validation. Defaults
#'     to \code{FALSE}.}
#' }
#' @param optim_opts A list, containing additional optional arguments controlling
#' the proximal gradient descent algorithm.
#' \describe{
#'     \item{check_obj}{A Boolean, if \code{TRUE}, convergence is determined
#'     by checking the decrease in the objective. Otherwise it is determined by
#'     checking the average entry-wise difference in consecutive values of \eqn{F}.
#'     Defaults to \code{TRUE}.}
#'     \item{eig_maxitr}{A positive integer, maximum iterations for internal
#'     eigenvalue solver. Defaults to \code{1000}.}
#'     \item{eig_prec}{A positive scalar, estimated eigenvalues below this
#'     threshold are set to zero. Defaults to \code{1e-2}.}
#'     \item{eps}{A positive scalar, convergence threshold for proximal gradient
#'     descent. Defaults to \code{1e-6}.}
#'     \item{eta}{A positive scalar, step size for proximal gradient descent.
#'     Defaults to \code{1} for the Gaussian model, \code{5} for the logistic
#'     model.}
#'     \item{init}{A string, initialization method. Valid options are
#'     \code{'fix'} (using initializers \code{optim_opts$V_init} and
#'     \code{optim_opts$U_init}), \code{'zero'} (initialize all parameters at zero),
#'      or \code{'svd'} (initialize with a truncated SVD with rank \code{optim_opts$init_rank}).
#'      Defaults to \code{'zero'}.}
#'     \item{K_max}{A positive integer, maximum iterations for proximal gradient
#'     descent. Defaults to \code{100}.}
#'     \item{max_rank}{A positive integer, maximum rank for internal eigenvalue
#'     solver. Defaults to \code{sqrt(n)}.}
#'     \item{missing_pattern}{An \eqn{n \times n \times m} Boolean array with \code{TRUE}
#'     for each observed entry and \code{FALSE} for missing entries. If unspecified, it
#'     is set to \code{!is.na(A)}.}
#'     \item{positive}{A Boolean, if \code{TRUE}, singular value thresholding only retains
#'     positive eigenvalues. Defaults to \code{FALSE}.}
#'     \item{return_posns}{A Boolean, if \code{TRUE}, returns estimates
#'     of the latent positions based on ASE. Defaults to \code{FALSE}.}
#'     \item{verbose}{A Boolean, if \code{TRUE}, console output will provide
#'     updates on the progress of proximal gradient descent. Defaults to
#'     \code{FALSE}.}
#' }
#'
#' @return A list is returned with the MultiNeSS model estimates, dimensions of
#' the common and individual latent spaces, and some additional optimization
#' output:
#' \item{F_hat}{An \eqn{n \times n} matrix estimating the common part of the expected
#' adjacency matrix, \eqn{F = VV^{\intercal}}. If \code{optim_opts$return_posns}
#' is \code{TRUE}, this is not returned.}
#' \item{G_hat}{A list of length \eqn{m}, the collection of \eqn{n \times n} matrices
#' estimating the individual part of each adjacency matrix, \eqn{G_k = U_kU_k^{\intercal}}.
#' If \code{optim_opts$return_posns}
#' is \code{TRUE}, this is not returned.}
#' \item{V_hat}{A matrix estimating the common latent positions.
#' Returned if \code{optim_opts$return_posns} is \code{TRUE}.}
#' \item{U_hat}{A list of length \eqn{m}, the collection of matrices
#' estimating the individual latent positions.
#' Returned if \code{optim_opts$return_posns} is \code{TRUE}.}
#' \item{d1}{A non-negative integer, the estimated common dimension of the
#' latent space.}
#' \item{d2}{An integer vector of length \eqn{m}, the estimated individual
#' dimension of the latent space for each layer.}
#' \item{K}{A positive integer, the number of iterations run in proximal
#' gradient descent.}
#' \item{convergence}{An integer convergence code, \code{0} if proximal
#' gradient descent converged in fewer than \code{optim_opts$K_max} iterations,
#' \code{1} otherwise.}
#' \item{lambda}{A positive scalar, the tuned \eqn{\lambda}
#' penalty parameter (see Details).}
#' \item{alpha}{A numeric vector of length \eqn{m}, the tuned \eqn{\alpha} penalty parameters
#' (see Details).}
#'
#' @examples
#' # gaussian model data
#' data1 <- multiness_sim(n=100,m=4,d1=2,d2=2,
#'                      model="gaussian")
#'
#' # multiness_fit with fixed tuning
#' fit1 <- multiness_fit(A=data1$A,
#'                       model="gaussian",
#'                       self_loops=TRUE,
#'                       refit=FALSE,
#'                       tuning="fixed",
#'                       tuning_opts=list(lambda=40,alpha=1/2),
#'                       optim_opts=list(max_rank=20,verbose=TRUE))
#'
#' # multiness_fit with adaptive tuning
#' fit2 <- multiness_fit(A=data1$A,
#'                       refit=TRUE,
#'                       tuning="adaptive",
#'                       tuning_opts=list(layer_wise=FALSE),
#'                       optim_opts=list(return_posns=TRUE))
#'
#' # logistic model data
#' data2 <- multiness_sim(n=100,m=4,d1=2,d2=2,
#'                        model="logistic",
#'                        self_loops=FALSE)
#'
#' # multiness_fit with cv tuning
#' fit3 <- multiness_fit(A=data2$A,
#'                       model="logistic",
#'                       self_loops=FALSE,
#'                       tuning="cv",
#'                       tuning_opts=list(N_cv=2,
#'                                        penalty_const_vec=c(1,2,2.309,3),
#'                                        verbose_cv=TRUE))
#'
#' @export
multiness_fit <- function(
  A, # adjacency array (n x n x m)
  model="gaussian", # model (gaussian or logistic for now)
  self_loops=TRUE, # does the data contain self_loops
  refit=TRUE,
  tuning="adaptive", # type of tuning (adaptive, fixed or cv)
  # additional options
  tuning_opts=list(),
  optim_opts=list()
){
  # get dimensions (PWM: fill these below later)
  n <- dim(A)[1]
  m <- dim(A)[3]

  # check dimensions of the adjacency array:
  if(!(length(dim(A))==3 & n==dim(A)[2])){
    stop("Invalid dimensions")
  }
  # check model:
  if(!(model=="gaussian" || model=="Gaussian" || model=="logistic")){
    stop("Not a valid model")
  }
  # check tuning:
  if(!(tuning=="fixed" || tuning=="adaptive" || tuning=="cv")){
    stop("Not a valid tuning method")
  }

  # set default optimization options:

  # check_obj
  if(is.null(optim_opts$check_obj)){
    optim_opts$check_obj <- TRUE
  }
  # eig_maxitr
  if(is.null(optim_opts$eig_maxitr)){
    optim_opts$eig_maxitr <- 1000
  }
  # eig_prec
  if(is.null(optim_opts$eig_prec)){
    optim_opts$eig_prec <- 1e-2
  }
  # eps
  if(is.null(optim_opts$eps)){
    optim_opts$eps <- 1e-6
  }
  # eta
  if(is.null(optim_opts$eta)){
    if(model == "gaussian"){
      optim_opts$eta <- 1
    }
    if(model == "logistic"){
      optim_opts$eta <- 5
    }
  }
  # init (and V_init/U_init for fix init, init_rank for svd)
  if(is.null(optim_opts$init)){
    optim_opts$init <- "zero"
  }
  else{
    if(optim_opts$init=="svd"){
      if(is.null(optim_opts$init_rank)){
        optim_opts$init_rank <- round(sqrt(n))
      }
    }
    if(optim_opts$init=="fix"){
      if(is.null(optim_opts$U_init) || is.null(optim_opts$V_init)){
        optim_opts$init=="zero"
      }
    }
  }
  # K_max
  if(is.null(optim_opts$K_max)){
    optim_opts$K_max <- 100
  }
  # max_rank
  if(is.null(optim_opts$max_rank)){
    optim_opts$max_rank <- round(sqrt(n))
  }
  # missing_pattern
  if(is.null(optim_opts$missing_pattern)){
    # set to NAs if no pattern is specified
    if(any(is.na(A))){
      optim_opts$missing_pattern <- !is.na(A)
    }
  }
  # positive
  if(is.null(optim_opts$positive)){
    optim_opts$positive <- FALSE
  }
  # return_posns
  if(is.null(optim_opts$return_posns)){
    optim_opts$return_posns <- FALSE
  }
  # verbose
  if(is.null(optim_opts$verbose)){
    optim_opts$verbose <- FALSE
  }

  # set link functions and family variable for glm
  if(model=="gaussian" || model=="Gaussian"){
    linkl <- list(identity,identity)
    linkf <- identity
    family_str <- "gaussian"
  }
  if(model=="logistic"){
    linkl <- list(expit,logit)
    linkf <- expit
    family_str <- "binomial"
  }

  # set default tuning options and tune the different methods
  # each of these chunks define tuned parameters:
  # 'lambda_vec_tuned'
  # 'alpha_tuned'

  # FIXED tuning
  if(tuning=="fixed"){
    # set default options
    # lambda PWM: note that this depends on the model choice
    if(is.null(tuning_opts$lambda)){
      if(model=="gaussian"){
        tuning_opts$lambda <- (4/sqrt(3))*sqrt(n)*sqrt(m)
      }
      else{
        tuning_opts$lambda <- .5*(4/sqrt(3))*sqrt(n)*sqrt(m)
      }
    }
    # alpha
    if(is.null(tuning_opts$alpha)){
      tuning_opts$alpha <- 1/sqrt(m)
    }

    # set parameters
    if(!(length(tuning_opts$alpha)==m)){
      if(length(tuning_opts$alpha)==1){
        alpha_vec_tuned <- rep(tuning_opts$alpha,m)
      }
      else{
        stop('alpha must be either a scalar or a vector of length m')
      }
    }
    else{
      alpha_vec_tuned <- tuning_opts$alpha
    }
    lambda_tuned <- tuning_opts$lambda
  }

  # ADAPTIVE tuning
  if(tuning=="adaptive"){
    # set defaults
    # layer_wise
    if(is.null(tuning_opts$layer_wise)){
      tuning_opts$layer_wise <- TRUE
    }
    # penalty_const
    if(is.null(tuning_opts$penalty_const)){
      tuning_opts$penalty_const <- 4/sqrt(3)
    }
    # penalty_const_alpha
    if(is.null(tuning_opts$penalty_const_lambda)){
      tuning_opts$penalty_const_lambda <- 1
    }

    # set parameters
    sigma_hat_vec <- apply(A,3,denoiseR::estim_sigma,method="MAD")
    if(!tuning_opts$layer_wise){
      sigma_hat_vec <- rep(mean(sigma_hat_vec),m)
    }
    # euclidean norm of sigma_hat_vec
    sigma_hat_en <- sqrt(sum(sigma_hat_vec^2))
    # final tuning parameters
    lambda_tuned <- tuning_opts$penalty_const*tuning_opts$penalty_const_lambda*sqrt(n)*sigma_hat_en
    alpha_vec_tuned <- sigma_hat_vec / (tuning_opts$penalty_const_lambda*sigma_hat_en)
  }

  # CROSS-VALIDATION tuning
  if(tuning=="cv"){
    # set default options
    # layer_wise
    if(is.null(tuning_opts$layer_wise)){
      tuning_opts$layer_wise <- TRUE
    }
    # N_cv (CV stability reps)
    if(is.null(tuning_opts$N_cv)){
      tuning_opts$N_cv <- 3
    }
    # p_cv (missing proportion)
    if(is.null(tuning_opts$p_cv)){
      tuning_opts$p_cv <- .1
    }
    # penalty_const_alpha
    if(is.null(tuning_opts$penalty_const_lambda)){
      tuning_opts$penalty_const_lambda <- 1
    }
    # penalty_const_vec (range of tuning constants)
    if(is.null(tuning_opts$penalty_const_vec)){
      tuning_opts$penalty_const_vec <- seq(1,4,.5)
    }
    # refit_cv (refit for CV prediction)
    if(is.null(tuning_opts$refit_cv)){
      tuning_opts$refit_cv <- TRUE
    }
    # verbose_cv
    if(is.null(tuning_opts$verbose_cv)){
      tuning_opts$verbose_cv <- FALSE
    }

    # set parameters
    # calibrate CV (up to proportionality)
    sigma_hat_vec <- apply(A,3,denoiseR::estim_sigma,method="MAD")
    if(!tuning_opts$layer_wise){
      sigma_hat_vec <- rep(mean(sigma_hat_vec),m)
    }
    # euclidean norm of sigma_hat_vec
    sigma_hat_en <- sqrt(sum(sigma_hat_vec^2))
    # set lambda up to a constant
    lambda_propto <- tuning_opts$penalty_const_lambda*sqrt(n)*sigma_hat_en
    # set alpha without tuning
    alpha_vec_tuned <- sigma_hat_vec / (tuning_opts$penalty_const_lambda*sigma_hat_en)

    # edge cross-validation:
    # find the subset of non-zero entries for holdout (if applicable)
    # allocate CV result
    LL <- length(tuning_opts$penalty_const_vec)
    cv_results <- matrix(Inf,tuning_opts$N_cv,LL)
    # loop through CV repeats
    for(ii in 1:tuning_opts$N_cv){
      # generate missingness pattern
      # augment with existing missingness if applicable
      if(is.null(optim_opts$missing_pattern)){
        MP_cv <- missing_array_sym(n,m,tuning_opts$p_cv)
      }
      else{
        MP_cv <- missing_array_sym(n,m,tuning_opts$p_cv,
                                   subset=optim_opts$missing_pattern)
      }
      # augment with existing missingness if applicable
      # loop through lambdas
      for(jj in 1:LL){
        # set lambda
        lambda_cv_current <- tuning_opts$penalty_const_vec[jj]*lambda_propto
        # fit model
        fit_cv <- prox_gd_memeff(A=A,
                                 lambda=lambda_cv_current,
                                 alpha=alpha_vec_tuned,
                                 link=linkl,
                                 eta=optim_opts$eta,
                                 init=optim_opts$init,
                                 V_init=optim_opts$V_init,U_init=optim_opts$U_init,
                                 pos=optim_opts$positive,
                                 block=T,soft=T,
                                 hollow=!self_loops,
                                 misspattern = MP_cv,
                                 eps=optim_opts$eps,
                                 max_rank=optim_opts$max_rank,init_rank=optim_opts$init_rank,
                                 eig_maxitr=optim_opts$eig_maxitr,
                                 verbose=FALSE,
                                 check_obj=optim_opts$check_obj,
                                 eig_prec=optim_opts$eig_prec)
        # refit if applicable and evaluate error
        if(tuning_opts$refit_cv){
          # refit
          refit_cv <- prox_gd_refit(A=A,
                                    fit=fit_cv,
                                    family=family_str,
                                    misspattern=MP_cv,
                                    hollow=!self_loops,
                                    return_factors=FALSE,
                                    ignore_neg=FALSE)
          # cv error
          cv_results[ii,jj] <- holdout_error(A=A,
                                             fit=refit_cv,
                                             hollow=!self_loops,
                                             misspattern=MP_cv,
                                             link=linkf)

        }
        else{
          # cv error
          cv_results[ii,jj] <- holdout_error(A=A,
                                             fit=fit_cv,
                                             hollow=!self_loops,
                                             misspattern=MP_cv,
                                             link=linkf)
        }
        if(tuning_opts$verbose_cv){
          cat("\nCompleted CV rep",ii,"for constant",jj,"of",LL,"\n")
        }
      }
    }
    # evaluated CV errors
    # select constant
    lambda_const_select <- tuning_opts$penalty_const_vec[which.min(Matrix::colMeans(cv_results))]
    # set final tuning parameter(s)
    lambda_tuned <- lambda_const_select*lambda_propto
  }

  # fit the final model (first stage)
  fit <- prox_gd_memeff(A=A,
                        lambda=lambda_tuned,
                        alpha=alpha_vec_tuned,
                        link=linkl,
                        eta=optim_opts$eta,
                        init=optim_opts$init,
                        V_init=optim_opts$V_init,U_init=optim_opts$U_init,
                        pos=optim_opts$positive,
                        block=T,soft=T,
                        hollow=!self_loops,
                        misspattern = optim_opts$missing_pattern,
                        eps=optim_opts$eps,
                        max_rank=optim_opts$max_rank,init_rank=optim_opts$init_rank,
                        K_max=optim_opts$K_max,
                        eig_maxitr=optim_opts$eig_maxitr,
                        verbose=optim_opts$verbose,
                        check_obj=optim_opts$check_obj,
                        eig_prec=optim_opts$eig_prec)
  # if required, refit the model (second stage)
  if(refit){
    refit <- prox_gd_refit(A=A,
                           fit=fit,
                           family=family_str,
                           misspattern=optim_opts$missing_pattern,
                           hollow=!self_loops,
                           return_factors=FALSE,
                           ignore_neg=FALSE)
  }
  else{
    refit <- fit
  }
  # return results, factorize if required
  if(!optim_opts$return_posns){
    return(list(F_hat=refit$F_hat,
                G_hat=refit$G_hat,
                d1=refit$F_rank,
                d2=refit$G_rank,
                K=refit$K,
                convergence=refit$convergence,
                lambda=lambda_tuned,
                alpha=alpha_vec_tuned))
  }
  else{
    V_hat <- ase(refit$F_hat,refit$F_rank)
    U_hat <- lapply(1:m,
                    function(kk){
                      ase(refit$G_hat[[kk]],refit$G_rank[kk])
                    })
    return(list(V_hat=V_hat,
                U_hat=U_hat,
                d1=refit$F_rank,
                d2=refit$G_rank,
                K=refit$K,
                convergence=refit$convergence,
                lambda=lambda_tuned,
                alpha=alpha_vec_tuned))
  }
}


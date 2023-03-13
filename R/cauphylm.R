#' @title Maximum Likelihood estimator for a Cauchy model
#'
#' @description
#' Find the maximum likelihood, using numerical optimization with \code{\link{optim}}.
#'
#' @param formula a model formula.
#' @param data a data frame containing variables in the model.
#' If not found in data, the variables are taken from current environment.
#' @param phy a phylogenetic tree of class \code{\link[ape]{phylo}}.
#' @param model a model for the trait evolution. Only \code{"cauchy"} is allowed.
#' @param lower.bound named list with lower bound values for the parameters.
#' @param upper.bound named list with upper bound values for the parameters. See Details.
#' @param starting.value starting value for the parameters.
#' This should be a named list, with \code{disp} the dispersion parameter.
#' The default initial values are computed from standard statistics used on (independent) Cauchy variables, see details.
#' @param hessian if \code{TRUE}, then the numerical hessian is computed, for confidence interval computations. See \code{\link{compute_vcov}}.
# @param optim.algo optimization algorithms to be used in \code{\link[nloptr]{nloptr}}.
#' 
#' @details 
#' Unless specified by the user, the initial values for the parameters are taken according to the following heuristics:
#' \itemize{
#'  \item{\code{mu}}{ is the trimmed mean of the trait, keeping only 24\% of the observations, as advocated in Rothenberg et al. 1964 (for models with the intercept only)}
#'  \item{\code{coef}}{ are obtained from a robust regression using \code{\link[robustbase]{lmrob.S}} (for linear models with several predictors)}
#'  \item{\code{disp}}{ is initialized from the trait centered and normalized by tip heights, with one of the following statistics, taken from Rousseeuw & Croux 1993:}
#'  \itemize{
#'  \item{\code{IQR}}{ half of the inter-quartile range (see \code{\link{IQR}});}
#'  \item{\code{MAD}}{ median absolute deviation with constant equal to 1 (see \code{\link{mad}});}
#'  \item{\code{Sn}}{ Sn statistics with constant 0.7071 (see \code{\link[robustbase]{Sn}});}
#'  \item{\code{Qn}}{ Qn statistics with constant 1.2071 (see \code{\link[robustbase]{Qn}});}
#' }
#' }
#' 
#' For \code{model=lambda}, the maximum default value is computed using \code{phytools:::maxLambda},
#' and is the ratio between the maximum height of a tip node over the maximum height of an internal node.
#' This can be larger than 1.
#'
#' @return
#' \item{coefficients}{the named vector of coefficients.}
#' \item{disp}{the maximum likelihood estimate of the dispersion parameter.}
#' \item{logLik}{the maximum of the log likelihood.}
#' \item{p}{the number of all parameters of the model.}
#' \item{aic}{AIC value of the model.}
#' \item{fitted.values}{fitted values}
#' \item{residuals}{raw residuals}
#' \item{y}{response}
#' \item{X}{design matrix}
#' \item{n}{number of observations (tips in the tree)}
#' \item{d}{number of dependent variables}
#' \item{formula}{the model formula}
#' \item{call}{the original call to the function}
#' \item{model}{the phylogenetic model for the covariance}
#' 
#' @seealso \code{\link{fitCauchy}}, \code{\link[phylolm]{phylolm}}
#' 
#' @references
#' Rothenberg T. J., Fisher F. M., Tilanus C. B. 1964. A Note on Estimation from a Cauchy Sample. Journal of the American Statistical Association. 59:460–463.
#' 
#' @export
#'
cauphylm <- function(formula, data = list(), phy,
                     model = c("cauchy", "lambda"),
                     lower.bound = list(disp = 0, lambda = 0), 
                     upper.bound = list(disp = Inf, lambda = NULL),
                     starting.value = list(disp = NULL,
                                           lambda = NULL),
                     hessian = FALSE) {
  ## Model matrix
  mf <- model.frame(formula = formula, data = data)
  mf <- checkTraitTree(mf, phy)
  X = model.matrix(attr(mf, "terms"), data = mf)
  y = model.response(mf)
  d = ncol(X)
  n <- length(phy$tip.label)
  
  ## Check that response is a vector
  if (is.matrix(y) && ncol(y) > 1) {
    stop(paste0("The response vector y in the formula is multivariate (it has several columns).\n",
                "  Please fit each column one by one: 'cauphylm' can only handle a simple (univariate) response vector y."))
  }
  
  model <- match.arg(model)
  
  res <- fitCauchy.internal(phy, X, y, 
                            model = model,
                            method = "fixed.root",
                            starting.value = starting.value,
                            lower.bound = lower.bound, 
                            upper.bound = upper.bound)
  
  sol <- res$param
  coefs <- res$param[grepl("coef", names(res$param))]
  names(coefs) <- colnames(X)
  
  res <- list(coefficients = coefs,
              disp = safe_get(res$param, "disp"),
              logLik = res$logLikelihood,
              p = length(res$param),
              aic = 2 * (length(res$param)) - 2 * res$logLikelihood,
              fitted.values = drop(X %*% coefs),
              residuals = y - drop(X %*% coefs),
              y = y,
              X = X,
              n = n,
              d = d,
              formula = formula,
              call = match.call(),
              model = model,
              sigma2_error = 0,
              phy = phy,
              lambda = safe_get(res$param, "lambda"),
              method = "fixed.root",
              root_tip_reml = res$rootTip)
  
  ## vcov
  if (hessian) {
    res <- compute_vcov.cauphylm(res)
  }
  
  class(res) <- "cauphylm"
  return(res)
}

#' @title Compute Approximated Variance Covariance Matrix
#'
#' @description
#' Find the approximated vcov matrix of the parameters.
#'
#' @param obj a fitted object, either with \code{\link{cauphylm}} or \code{\link{fitCauchy}}
#' 
#' @details 
#' This function uses the numerical inverse of the Hessian of the likelihood at the optimal value
#' to approximate the variance covariance matrix.
#' It can be used to compute confidence intervals with functions \code{\link{confint.cauphylm}}
#' or \code{\link{confint.cauphyfit}}
#'
#' @return
#' The same object, with added vcov entry.
#' 
#' @export
#'
compute_vcov <- function(obj) {
  UseMethod("compute_vcov")
}

##
#' @export
#' @method compute_vcov cauphylm
##
compute_vcov.cauphylm <- function(obj) {
  param_names <- getParamNames(obj$model, obj$X)
  all_params_names <- c(names(obj$coefficients), param_names[!grepl("coef", param_names)])
  all_params <- c(obj$coefficients, obj$disp)
  if (obj$model == "lambda") all_params <- c(all_params, obj$lambda)
  names(all_params) <- all_params_names
  obj$all_params <- all_params
  # approxHessian <- numDeriv::hessian(minus_like, -obj$loglikelihood, tree = tree, trait = trait)
  # approxHessian3 <- pracma::hessian(minus_like, -obj$loglikelihood, tree = tree, trait = trait)
  minus_like <- function(param) {
    names(param) <- param_names
    centralTips <- drop(obj$X %*% param[grepl("coef", names(param))])
    phy_trans <- transformBranchLengths(obj$phy, obj$model, param)
    return(-logDensityTipsCauchy(phy_trans, obj$y - centralTips, 0, param["disp"], method = "fixed.root"))
  }
  # approxHessian <- nlme::fdHess(obj$all_params, minus_like)
  # obj$vcov <- solve(approxHessian$Hessian)
  approxHessian <- pracma::hessian(f = minus_like, x0 = obj$all_params)
  obj$vcov <- solve(approxHessian)
  colnames(obj$vcov) <- rownames(obj$vcov) <- all_params_names
  return(obj)
}

##
#' @export
#' @method print cauphylm
##
print.cauphylm <- function(x, digits = max(3, getOption("digits") - 3), ...){
  # Call
  cat("Call:\n")
  print(x$call)
  cat("\n")
  ## AIC
  aiclogLik = c(x$aic,x$logLik)
  names(aiclogLik) = c("AIC","logLik")
  print(aiclogLik, digits = digits)
  ## Parameters
  cat("\nParameter estimate(s) using ML:\n")
  cat("dispersion:",x$disp,"\n")
  if (!is.null(x$lambda)) cat(x$model,":",x$lambda, "\n")
  ## Coefficients
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

##
#' Generic Methods for cauphylm
#' 
#' @export
#' @inheritParams phylolm::vcov.phylolm
#' @method vcov cauphylm
vcov.cauphylm <- function(object, ...) {
  if (is.null(object$vcov)) {
    object <- compute_vcov.cauphylm(object)
  }
  return(object$vcov)
}
#' @export
#' @method logLik cauphylm
#' @rdname vcov.cauphylm
logLik.cauphylm <- function(object, ...){
  res <- phylolm::logLik.phylolm(object, ...)
  class(res) <- "logLik.cauphylm"
  res
}
#' @export
#' @inheritParams phylolm::print.logLik.phylolm
#' @method print logLik.cauphylm
#' @rdname vcov.cauphylm
print.logLik.cauphylm <- phylolm::print.logLik.phylolm
##
#' @export
#' @inheritParams phylolm::AIC.logLik.phylolm
#' @method AIC logLik.cauphylm
#' @rdname vcov.cauphylm
AIC.logLik.cauphylm <- phylolm::AIC.logLik.phylolm
#' @export
#' @method AIC cauphylm
#' @rdname vcov.cauphylm
AIC.cauphylm <- phylolm::AIC.phylolm
# #' @export
# #' @inheritParams phylolm::extractAIC.phylolm
# #' @rdname vcov.cauphylm
# extractAIC.cauphylm <- phylolm::extractAIC.phylolm
# #' @export
# #' @rdname vcov.cauphylm
# nobs.cauphylm <- phylolm::nobs.phylolm
#' @export
#' @inheritParams phylolm::predict.phylolm
#' @method predict cauphylm
#' @rdname vcov.cauphylm
predict.cauphylm <- phylolm::predict.phylolm
#' @export
#' @inheritParams stats::confint
#' @method confint cauphylm
#' @rdname vcov.cauphylm
confint.cauphylm <- function(object, parm, level = 0.95, ...){
  message("Approximated asymptotic confidence interval using the Hessian.")
  cf <- object$all_params
  ses <- sqrt(diag(vcov.cauphylm(object)))
  pnames <- names(ses)
  if (missing(parm)) 
    parm <- pnames
  else if (is.numeric(parm)) 
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- stats::qnorm(a)
  pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
  ci <- array(NA_real_, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ci[] <- cf[parm] + ses[parm] %o% fac
  ci
}
#' @title Phylogenetic Regression using a Cauchy Process
#'
#' @description
#' Perform a phylogenetic regression using the Cauchy Process, by numerical optimization.
#'
#' @param formula a model formula.
#' @param data a data frame containing variables in the model.
#' If not found in data, the variables are taken from current environment.
#' @param starting.value named list initial values for the parameters.
#' See Details for the default values.
#' @inheritParams fitCauchy
#' 
#' @details 
#' This function fits a Cauchy Process on the phylogeny, using maximum likelihood
#' and the \code{"fixed.root"} method (see \code{\link{fitCauchy}}).
#' It further assumes that the root value \code{x0} is a linear combination of the
#' covariables in \code{formula}.
#' The corresponding regression model is:
#' \deqn{Y = X \beta + E,}
#' with:
#' \itemize{
#' \item{\eqn{Y}}{ the vector of traits at the tips of the tree;}
#' \item{\eqn{X}}{ the regression matrix of covariables in \code{formula};}
#' \item{\eqn{\beta}}{ the vector of coefficients;}
#' \item{\eqn{E}}{ a centered error vector that is Cauchy distributed,
#' and can be seen as the result of a Cauchy process starting at 0 at the root,
#' and with a dispersion \code{disp} (see \code{\link{fitCauchy}}).}
#' }
#' 
#' Unless specified by the user, the initial values for the parameters 
#' are taken according to the following heuristics:
#' \itemize{
#'  \item{\code{coefficients}:}{ \eqn{\beta} are obtained from a robust regression using \code{\link[robustbase]{lmrob.S}};}
#'  \item{\code{disp}:}{ is initialized from the trait centered and normalized 
#'  by tip heights, with one of the following statistics, taken from Rousseeuw & Croux 1993:}
#'  \itemize{
#'  \item{\code{IQR}:}{ half of the inter-quartile range (see \code{\link{IQR}});}
#'  \item{\code{MAD}:}{ median absolute deviation with constant equal to 1 (see \code{\link{mad}});}
#'  \item{\code{Sn}:}{ Sn statistics with constant 0.7071 (see \code{\link[robustbase]{Sn}});}
#'  \item{\code{Qn}:}{ Qn statistics with constant 1.2071 (see \code{\link[robustbase]{Qn}}).}
#' }
#' }
#' 
#' Unless specified by the user, \code{disp} is taken positive unbounded.
#' 
#' The function uses \code{\link{nloptr}} for the numerical optimization of the 
#' (restricted) likelihood, computed with function \code{\link{logDensityTipsCauchy}}.
#' It uses algorithms \href{https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#bobyqa}{\code{BOBYQA}}
#' and \href{https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#mlsl-multi-level-single-linkage}{\code{MLSL_LDS}}
#' for local and global optimization.
#' 
#' If \code{model="lambda"}, the CP is fit on a tree with branch lengths re-scaled 
#' using the Pagel's lambda transform (see \code{\link[phylolm]{transf.branch.lengths}}),
#' and the \code{lambda} value is estimated using numerical optimization.
#' The default initial value for the \code{lambda} parameter is computed using adequate robust moments.
#' The default maximum value is computed using \code{phytools:::maxLambda},
#' and is the ratio between the maximum height of a tip node over the maximum height of an internal node.
#' This can be larger than 1.
#' The default minimum value is 0.
#'
#' @return
#' \item{coefficients}{the named vector of estimated coefficients.}
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
#' \item{phy}{the phylogenetic tree}
#' \item{lambda}{the ml estimate of the lambda parameter (for \code{model="lambda"})}
#' 
#' @seealso \code{\link{fitCauchy}}, \code{\link{confint.cauphylm}}, \code{\link{ancestral}}, \code{\link{increment}}, \code{\link{logDensityTipsCauchy}}, \code{\link[phylolm]{phylolm}}
#' 
#' @examples
#' # Simulate tree and data
#' set.seed(1289)
#' phy <- ape::rphylo(20, 0.1, 0)
#' error <- rTraitCauchy(n = 1, phy = phy, model = "cauchy",
#'                       parameters = list(root.value = 0, disp = 0.1))
#' x1 <- ape::rTraitCont(phy, model = "BM", sigma = 0.1, root.value = 0)
#' trait <- 3 + 2*x1 + error
#' # Fit the data
#' fit <- cauphylm(trait ~ x1, phy = phy)
#' fit
#' # Approximate confidence intervals
#' confint(fit)
#' 
#' @references
#' Bastide, P. and Didier, G. 2023. The Cauchy Process on Phylogenies: a Tractable Model for Pulsed Evolution. Systematic Biology. doi:10.1093/sysbio/syad053.
#' 
#' Rothenberg T. J., Fisher F. M., Tilanus C. B. 1964. A Note on Estimation from a Cauchy Sample. Journal of the American Statistical Association. 59:460–463.
#' 
#' Rousseeuw P.J., Croux C. 1993. Alternatives to the Median Absolute Deviation. Journal of the American Statistical Association. 88:1273–1283.
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
  # Checks
  check_binary_tree(phy)
  
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
#' Find the approximated variance covariance matrix of the parameters.
#'
#' @param obj a fitted object, either with \code{\link{fitCauchy}} or \code{\link{cauphylm}}.
#' 
#' @details 
#' This function computes the numerical Hessian of the likelihood at the optimal value
#' using function \code{\link[pracma]{hessian}}, and then uses its inverse 
#' to approximate the variance covariance matrix.
#' It can be used to compute confidence intervals with functions \code{\link{confint.cauphylm}}
#' or \code{\link{confint.cauphyfit}}.
#' 
#' \code{\link{confint.cauphylm}} and \code{\link{confint.cauphyfit}}
#' internally call \code{compute_vcov}, but do not save the result.
#' This function can be used to save the vcov matrix.
#'
#' @return
#' The same object, with added vcov entry.
#' 
#' @seealso \code{\link{fitCauchy}}, \code{\link{cauphylm}}, 
#' \code{\link{confint.cauphylm}}, \code{\link{confint.cauphyfit}},
#' \code{\link{vcov.cauphylm}}, \code{\link{vcov.cauphyfit}}
#' 
#' @examples
#' # Simulate tree and data
#' set.seed(1289)
#' phy <- ape::rphylo(20, 0.1, 0)
#' dat <- rTraitCauchy(n = 1, phy = phy, model = "cauchy",
#'                     parameters = list(root.value = 10, disp = 0.1))
#' # Fit the data, without computing the Hessian at the estimated parameters.
#' fit <- fitCauchy(phy, dat, model = "cauchy", method = "reml", hessian = FALSE)
#' # Precompute the vcov matrix
#' fit <- compute_vcov(fit)
#' # Approximate confidence intervals
#' confint(fit)
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
    return(-logDensityTipsCauchy(phy_trans, obj$y - centralTips, 0, param["disp"], method = "fixed.root", do_checks = FALSE))
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
#' @inheritParams phylolm::print.phylolm
#' @rdname vcov.cauphylm
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
#' Generic Methods for S3 class \code{cauphylm}.
#' 
#' @export
#' @param object an object of class \code{cauphylm}.
#' 
#' @return
#' Same value as the associated methods from the \code{stats} package:
#' \itemize{
#' \item{\code{\link[stats]{vcov}}}{ an estimated covariance matrix, see \code{\link{compute_vcov}};}
#' \item{\code{\link[stats]{logLik}}}{ an object of class \code{\link[stats]{logLik}};}
#' \item{\code{\link[stats]{AIC}}}{ a numeric value;}
#' \item{\code{\link[stats]{confint}}}{ a matrix (or vector) with columns giving lower and upper confidence limits for each parameter;}
#' \item{\code{\link[stats]{coef}}}{ coefficients extracted from the model;}
#' \item{\code{\link[stats]{predict}}}{ a vector of predicted values.}
#' }
#' 
#' @examples
#' # Simulate tree and data
#' set.seed(1289)
#' phy <- ape::rphylo(20, 0.1, 0)
#' error <- rTraitCauchy(n = 1, phy = phy, model = "cauchy",
#'                       parameters = list(root.value = 0, disp = 0.1))
#' x1 <- ape::rTraitCont(phy, model = "BM", sigma = 0.1, root.value = 0)
#' trait <- 3 + 2*x1 + error
#' # Fit the data
#' fit <- cauphylm(trait ~ x1, phy = phy)
#' fit
#' # vcov matrix
#' vcov(fit)
#' # Approximate confidence intervals
#' confint(fit)
#' # log likelihood of the fitted object
#' logLik(fit)
#' # AIC of the fitted object
#' AIC(fit)
#' # predicted values
#' predict(fit)
#' # coefficients
#' coef(fit)
#' 
#' @seealso \code{\link{cauphylm}}, \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}
#' \code{\link[stats]{AIC}}, \code{\link[stats]{confint}}, \code{\link[stats]{coef}},
#' \code{\link[stats]{predict}}, \code{\link[phylolm]{predict.phylolm}}
#' 
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
#' @method print logLik.cauphylm
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
  if (is.null(object$vcov)) {
    object <- compute_vcov.cauphylm(object)
  }
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
#' @export
#' @method coef cauphylm
#' @rdname vcov.cauphylm
coef.cauphylm <- function(object, ...){
  ## Regression matrix for fixed root
  param_names <- getParamNames(object$model, object$X)
  all_params_names <- c(names(object$coefficients), param_names[!grepl("coef", param_names)])
  all_params <- c(object$coefficients, object$disp)
  if (object$model == "lambda") all_params <- c(all_params, object$lambda)
  names(all_params) <- all_params_names
  return(all_params)
}
#' @title Maximum Likelihood estimator for a Cauchy model
#'
#' @description
#' Find the maximum likelihood, using numerical optimization with \code{\link{optim}}.
#'
#' @param phy a phylogenetic tree of class \code{\link[ape]{phylo}}.
#' @param trait named vector of traits at the tips.
#' @param model a model for the trait evolution. Only \code{"cauchy"} is allowed.
#' @param starting.value starting value for the parameters of the Cauchy.
#' This should be a named list, with \code{mu} and \code{disp} the ancestral central and dispersion parameters.
#' The default initial values are computed from standard statistics used on (independent) Cauchy variables, see details.
#' @param lower.bound named list with lower bound values for the parameters.
#' @param upper.bound named list with upper bound values for the parameters. See Details.
#' @param method the method used to fit the process. One of \code{reml} (the default), \code{random.root} or \code{fixed.root}.
#' See details.
#' @param root.edge multiplicative factor for the root dispersion, equal to the length of the root edge. Ignored if \code{method!=random.root}.
#' @param hessian if \code{TRUE}, then the numerical hessian is computed, for confidence interval computations. See \code{\link{compute_vcov}}.
#' @param optim if "local", only a local optimization around the initial parameter values is performed (the default). If "global", a global maximization is attempted using the "MLSL" approach (see \code{\link{nloptr}}).
#' @param method.init.disp the initialization method for the dispersion. One of "Qn", "Sn", "MAD", "IQR". Default to the "Qn" statistics. See Details.
#' 
#' @details 
#' Unless specified by the user, the initial values for the parameters are taken according to the following heuristics:
#' \itemize{
#'  \item{\code{mu}}{ is the trimmed mean of the trait, keeping only 24\% of the observations, as advocated in Rothenberg et al. 1964 (for \code{method="fixed.root"})}
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
#' The \code{method} argument specifies the method used for the fit:
#' \itemize{
#'   \item{\code{method=reml}}{ 
#'   the dispersion parameter is fitted using the REML criterion,
#'   obtained by re-rooting the tree to one of the tips. 
#See \code{\link{logREMLTipsCauchy}}.
#'   }
#'   \item{\code{method=random.root}}{ 
#'   the root value is assumed to be a random Cauchy variable, centered at \code{mu=0},
#'   and with a dispersion \code{disp_root = disp * root.edge}.
#'   }
#'   \item{\code{method=fixed.root}}{ 
#'   the model is fitted conditionally on the root value mu,
#'   i.e. with a model where the root value is fixed and inferred from the data.
#'   }
#' }
#' In the first two cases, the optimization is done on the dispersion only,
#' while in the last case the optimization is on the root value and the dispersion.
#'
#' @return A list, with the maximum likelihood rate parameter, and the likelihood value.
#' 
#' @examples
#' # Simulate tree and data
#' phy <- ape::rphylo(5, 0.1, 0)
#' dat <- rTraitCauchy(n = 1, phy = phy, model = "cauchy", parameters = list(root.value = 0, disp = 1))
#' # Fit the data
#' fit <- fitCauchy(phy, dat, model = "cauchy", method = "reml")
#' fit
#' 
#' @references
#' Rothenberg T. J., Fisher F. M., Tilanus C. B. 1964. A Note on Estimation from a Cauchy Sample. Journal of the American Statistical Association. 59:460–463.
#' Rousseeuw P.J., Croux C. 1993. Alternatives to the Median Absolute Deviation. Journal of the American Statistical Association. 88:1273–1283.
#' 
#' @seealso \code{\link{cauphylm}}, \code{geiger::fitContinuous}
#' 
#' @export
#'
fitCauchy <- function(phy, trait, 
                      model = c("cauchy", "lambda"),
                      method = c("reml", "random.root", "fixed.root"),
                      starting.value = list(mu = NULL, disp = NULL, lambda = NULL),
                      lower.bound = list(disp = 0, lambda = 0), 
                      upper.bound = list(disp = Inf, lambda = NULL),
                      root.edge = 100,
                      hessian = FALSE,
                      optim = c("local", "global"),
                      method.init.disp = c("Qn", "Sn", "MAD", "IQR")) {
  
  model <- match.arg(model)
  method <- match.arg(method)
  optim <- match.arg(optim)
  method.init.disp <- match.arg(method.init.disp)
  
  res <- fitCauchy.internal(phy, X = NULL, trait, 
                            model = model,
                            method = method,
                            starting.value = starting.value,
                            lower.bound = lower.bound, 
                            upper.bound = upper.bound,
                            root.edge = root.edge,
                            optim = optim,
                            method.init.disp = method.init.disp)
  
  res <- list(x0 = safe_get(res$param, "coef1"),
              disp = safe_get(res$param, "disp"),
              lambda = safe_get(res$param, "lambda"),
              logLik = res$logLikelihood,
              p = length(res$param),
              aic = 2 * length(res$param) - 2 * res$logLikelihood,
              trait = trait,
              y = trait,
              n = Ntip(phy),
              d = 1,
              call = match.call(),
              model = model,
              sigma2_error = 0,
              phy = phy,
              method = method,
              random.root = (method == "random.root"),
              reml = (method == "reml"),
              root_tip_reml = res$rootTip)
  if (method == "random.root") {
    res$x0 <- 0.0
    res$phy$root.edge <- root.edge
  }
  
  ## vcov
  if (hessian) {
    res <- compute_vcov.cauphyfit(res)
  }
  
  class(res) <- "cauphyfit"
  return(res)
}

#' @export
#' @method compute_vcov cauphyfit
##
compute_vcov.cauphyfit <- function(obj) {
  ## Regression matrix for fixed root
  X <- NULL
  if (is.null(X) && obj$method == "fixed.root") {
    X <- matrix(rep(1, length(obj$y)), nrow = length(obj$y))
    colnames(X) <- "coef1"
  }
  # parameters
  param_names <- getParamNames(obj$model, X)
  all_params_names <- sub("coef1", "x0", param_names)
  all_params <- obj$disp
  if (obj$method == "fixed.root") all_params <- c(obj$x0, all_params)
  if (obj$model == "lambda") all_params <- c(all_params, obj$lambda)
  names(all_params) <- all_params_names
  obj$all_params <- all_params
  # likelihood
  minus_like <- switch(obj$method,
                       reml = minusLikelihoodREML,
                       fixed.root = minusLikelihoodFixedRoot(X),
                       random.root = minusLikelihoodRandomRoot)
  minus_like_untransformed <- function(param, param_names, ...) {
    names(param) <- param_names
    param <- transform_values(param)
    return(minus_like(param, param_names, ...))
  }
  # approxHessian <- nlme::fdHess(pars = obj$all_params, fun = minus_like_untransformed,
                                # param_names = param_names, tree = obj$phy, trait = obj$trait, Xdesign = X, model = obj$model, rootTip = obj$root_tip_reml)
  # approxHessian <- transform_hessian(approxHessian, obj$all_params)
  # obj$vcov <- solve(approxHessian$Hessian)
  approxHessian <- pracma::hessian(f = minus_like_untransformed, x0 = obj$all_params,
                                   param_names = param_names, tree = obj$phy, trait = obj$trait, Xdesign = X, model = obj$model, rootTip = obj$root_tip_reml)
  obj$vcov <- solve(approxHessian)
  colnames(obj$vcov) <- rownames(obj$vcov) <- all_params_names
  return(obj)
}

##
#' @export
#' @method print cauphyfit
##
print.cauphyfit <- function(x, digits = max(3, getOption("digits") - 3), ...){
  # Call
  cat("Call:\n")
  print(x$call)
  cat("\n")
  ## AIC
  aiclogLik = c(x$aic, x$logLik)
  if (x$reml) {
    names(aiclogLik) = c("AIC", "logLik (restricted)")
  } else {
    names(aiclogLik) = c("AIC", "logLik")
  }
  print(aiclogLik, digits = digits)
  ## Parameters
  if (x$reml) {
    cat("\nParameter estimate(s) using REML:\n")
  } else {
    cat("\nParameter estimate(s) using ML:\n")
  }
  cat("dispersion:", x$disp, "\n")
  if (!is.null(x$x0)) cat("root value:", x$x0, "\n")
  if (!is.null(x$lambda)) cat(x$model,":",x$lambda, "\n")
}

##
#' Generic Methods for cauphyfit
#' 
#' @export
#' @inheritParams phylolm::vcov.phylolm
#' @method vcov cauphyfit
vcov.cauphyfit <- function(object, ...) {
  if (is.null(object$vcov)) {
    object <- compute_vcov.cauphyfit(object)
  }
  return(object$vcov)
}
#' @export
#' @method logLik cauphyfit
#' @rdname vcov.cauphyfit
logLik.cauphyfit <- function(object, ...){
  res <- phylolm::logLik.phylolm(object, ...)
  class(res) <- "logLik.cauphyfit"
  res
}
#' @export
#' @inheritParams phylolm::print.logLik.phylolm
#' @method print logLik.cauphyfit
#' @rdname vcov.cauphyfit
print.logLik.cauphyfit <- phylolm::print.logLik.phylolm
##
#' @export
#' @inheritParams phylolm::AIC.logLik.phylolm
#' @method AIC logLik.cauphyfit
#' @rdname vcov.cauphyfit
AIC.logLik.cauphyfit <- phylolm::AIC.logLik.phylolm
#' @export
#' @method AIC cauphyfit
#' @rdname vcov.cauphyfit
AIC.cauphyfit <- phylolm::AIC.phylolm
# #' @export
# #' @inheritParams phylolm::extractAIC.phylolm
# #' @rdname vcov.cauphyfit
# extractAIC.cauphyfit <- phylolm::extractAIC.phylolm
# #' @export
# #' @rdname vcov.cauphyfit
# nobs.cauphyfit <- phylolm::nobs.phylolm
# #' @export
# #' @inheritParams phylolm::predict.phylolm
# #' @rdname vcov.cauphyfit
# predict.cauphyfit <- phylolm::predict.phylolm
#' @export
#' @inheritParams stats::confint
#' @method confint cauphyfit
#' @rdname vcov.cauphyfit
confint.cauphyfit <- function(object, parm, level = 0.95, ...){
  message("Approximated asymptotic confidence interval using the Hessian.")
  if (is.null(object$vcov)) {
    object <- compute_vcov.cauphyfit(object)
  }
  ses <- sqrt(diag(vcov.cauphyfit(object)))
  pnames <- names(ses)
  cf <- object$all_params
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
#' @importFrom stats reorder
#' @export
#' @method coef cauphyfit
#' @rdname vcov.cauphyfit
coef.cauphyfit <- function(object, ...){
  ## Regression matrix for fixed root
  X <- NULL
  if (is.null(X) && object$method == "fixed.root") {
    X <- matrix(rep(1, length(object$y)), nrow = length(object$y))
    colnames(X) <- "coef1"
  }
  # parameters
  param_names <- getParamNames(object$model, X)
  all_params_names <- sub("coef1", "x0", param_names)
  all_params <- object$disp
  if (object$method == "fixed.root") all_params <- c(object$x0, all_params)
  if (object$model == "lambda") all_params <- c(all_params, object$lambda)
  names(all_params) <- all_params_names
  return(all_params)
}

#' @title Method for Profiling \code{cauphyfit} Objects
#'
#' @method profile cauphyfit
#' 
#' @description
#' Investigates the profile log-likelihood function for a fitted model of class \code{cauphyfit}.
#' 
#' @param fitted the \code{cauphyfit} fitted model object.
#' @param which the original model parameters which should be profiled. This can be a numeric or character vector. By default, all parameters are profiled.
#' @param level highest confidence level for parameters intervals, computed using the approximated Hessian (see \code{\link{compute_vcov}}).
#' @param npoints number of points to profile the likelihood for each parameter.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details 
#' This function computes a confidence interval for the parameters using \code{\link{confint.cauphyfit}},
#' and then compute the likelihood function between the bounds of the interval, for each parameter,
#' all other parameters being fixed.
#'
#' @return An object of class \code{profile.cauphyfit}, which is a list with an element for each parameter being profiled.
#' The elements are data-frames with two variables
#'  * par.vals a matrix of parameter values for each fitted model.
#'  * profLogLik	the profile log likelihood.
#' 
#' @examples
#' phy <- ape::rphylo(5, 0.1, 0)
#' dat <- rTraitCauchy(n = 1, phy = phy, model = "cauchy", parameters = list(root.value = 0, disp = 1))
#' fit <- fitCauchy(phy, dat, model = "cauchy", method = "reml")
#' pr <- profile(fit)
#' plot(pr)
#' 
#' @seealso \code{\link{fitCauchy}}, \code{\link{profile}}, \code{\link{plot.profile.cauphyfit}}
#' 
#' @export
#'
profile.cauphyfit <- function(fitted, which = 1:npar, level = 0.80, npoints = 100, ...){
  # Fitted object
  fitted <- compute_vcov(fitted)
  ci <- suppressMessages(confint.cauphyfit(fitted, level = level))
  ci["disp", 1] <- max(0, ci["disp", 1]) # make sure lower bound is larger than 0
  estim <- fitted$all_params
  names_params <- names(estim)

  # likelihood
  X <- NULL
  if (is.null(X) && fitted$method == "fixed.root") {
    X <- matrix(rep(1, length(fitted$y)), nrow = length(fitted$y))
    colnames(X) <- "coef1"
  }
  minus_like <- switch(fitted$method,
                       reml = minusLikelihoodREML,
                       fixed.root = minusLikelihoodFixedRoot(X),
                       random.root = minusLikelihoodRandomRoot)
  like_untransformed <- function(param, param_names, ...) {
    names(param) <- param_names
    param <- transform_values(param)
    return(- minus_like(param, param_names, ...))
  }

  # select params
  npar <- length(names_params)
  if (is.character(which)) {
    whichind <- match(which, names_params)
    if (anyNA(whichind)) message(paste0("Parameters ",
                                        paste0(which[is.na(whichind)], collapse = ", "),
                                        " are not in the fitted object, and will be ignored."))
    which <- whichind[!is.na(whichind)]
  }
  if (any(which > npar)) {
    message(paste0("Parameters with indexes ",
                   paste0(which[which > npar], collapse = ", "),
                   " are not in the fitted object, and will be ignored."))
    which <- which[which <= npar]
  }
  
  # result
  par.vals.estim <- rep(1, npoints + 1) %*% t(estim)
  par.vals.estim <- as.data.frame(par.vals.estim)
  res <- list()
  for (param in which) {
    grid <- seq(ci[param, 1], ci[param, 2], length.out = npoints)
    grid <- c(grid, estim[param])
    grid <- sort(grid)
    par.vals <- par.vals.estim
    par.vals[, param] <- grid
    profLogLik <- apply(par.vals, 1, like_untransformed,
                        param_names = names_params,
                        tree = fitted$phy, trait = fitted$trait, Xdesign = X, model = fitted$model, rootTip = fitted$root_tip_reml)
    res[[names_params[param]]] <- list(par.vals = par.vals,
                                       profLogLik = profLogLik)
  }
  class(res) <- "profile.cauphyfit"
  return(res)
}

#' @importFrom graphics abline
NULL

##
#' @title Plot for class \code{profile.cauphyfit}
#'
#' @description
#' This function takes an object of class code{\link{profile.cauphyfit}},
#' and plots the profile likelihood for each parameter.
#'
#' @param x an object of class \code{profile.cauphyfit}
#' @param n.col the number of columns on which to display the plot. Can be left blank.
#' @param ... further arguments to be passed to \code{\link{plot}}.
#' 
#' @return
#' NULL
#' 
#' @examples
#' phy <- ape::rphylo(5, 0.1, 0)
#' dat <- rTraitCauchy(n = 1, phy = phy, model = "cauchy", parameters = list(root.value = 0, disp = 1))
#' fit <- fitCauchy(phy, dat, model = "cauchy", method = "fixed.root")
#' pr <- profile(fit)
#' plot(pr)
#' 
#' @seealso \code{\link{profile.cauphyfit}}, \code{\link{fitCauchy}}
#' 
#' @export
#'
plot.profile.cauphyfit <- function(x, n.col, ...){
  nparams <- length(x)
  if (missing(n.col)) {
    n.col <- ifelse(nparams <= 3, nparams, 3)
  }
  n.lines <- (nparams %/% n.col) + ifelse(nparams %% n.col == 0, 0, 1)
  y_lab <- "Profile Log Likelihood"
  scr <- split.screen(c(n.lines, n.col))
  on.exit(close.screen(all.screens = TRUE))
  for (i in seq_len(nparams)) {
    screen(scr[i])
    plot(x[[i]]$par.vals[, names(x)[i]],
         x[[i]]$profLogLik,
         xlab = names(x)[i],
         ylab = y_lab, type = 'l', ...)
    abline(v = x[[i]]$par.vals[which.max(x[[i]]$profLogLik), names(x)[i]], lty = 2)
  }
}

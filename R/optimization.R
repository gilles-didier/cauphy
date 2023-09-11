#' @importFrom stats optim model.frame model.matrix model.response IQR median quantile mad
NULL

#' @importFrom robustbase lmrob lmrob.control
NULL

#' @importFrom grDevices dev.off hcl.colors png rainbow
NULL

#' @importFrom utils setTxtProgressBar txtProgressBar
NULL

#' @title Maximum Likelihood estimator for a Cauchy model
#'
#' @description
#' Find the maximum likelihood, using numerical optimization with \code{\link{optim}}.
#'
#' @inheritParams fitCauchy
#' @inheritParams cauphylm
#'
#' @return A list, with the maximum likelihood rate parameter, and the likelihood value.
#' 
#' @references
#' Rothenberg T. J., Fisher F. M., Tilanus C. B. 1964. A Note on Estimation from a Cauchy Sample. Journal of the American Statistical Association. 59:460–463.
#' 
#' @seealso cauphylm
#' 
#' @keywords internal
#'
fitCauchy.internal <- function(phy, X, y, 
                               model = c("cauchy", "lambda"),
                               method = c("reml", "random.root", "fixed.root"),
                               number.params,
                               starting.value = list(x0 = NULL, disp = NULL, lambda = NULL),
                               lower.bound = list(disp = 0, lambda = 0), 
                               upper.bound = list(disp = Inf, lambda = 1),
                               root.edge = 100,
                               optim = c("local", "global"),
                               method.init.disp = "Qn", ...) {
  # Checks
  check_binary_tree(phy)
  
  if (is.null(dim(y))) y <- as.matrix(y) # y is always a matrix with as many rows as the number of species
  y <- checkTraitTree(y, phy)
  
  if (ncol(y) == 1) checkDuplicates(y, phy) # If univarite, check replicates
  
  model <- match.arg(model)
  method <- match.arg(method)
  
  ## Root edge
  if (method == "random.root") {
    if (!is.null(phy$root.edge)) warning("The root edge of the tree is specified, but will be overridden by the `root.edge` argument of the function.")
    phy$root.edge <- root.edge
  }
  
  ## Regression matrix for fixed root
  if (is.null(X) && method == "fixed.root") {
    X <- matrix(rep(1, ncol(y)), nrow = ncol(y))
    colnames(X) <- "coef"
  }
  
  if (!(is.null(X) || ncol(X) == 1) && ncol(y) > 1) stop("Regression is only available for univariate traits.")
  
  if (model == "lambda" && ncol(y) > 1) stop("The Cauchy lambda model is only available for a univariate trait.")
  
  # Number of parameters
  if (missing(number.params)) {
    nmean <- max(ifelse(is.null(X), 0, ncol(X)),
                 ifelse(method == "fixed.root", ncol(y), 0))
    number.params <- list(mean = nmean,
                          angle = ncol(y)*(ncol(y)-1),
                          disp = ncol(y),
                          lambda = ifelse(model=="lambda", 1, 0))
  } else {
    # TODO: checks of number of parameters
  }
  
  # starting values
  start.values <- getStartingValues(model, phy, X, y, number.params, starting.value, method.init.disp, method)
  start.values <- transform_values(start.values)
  # lower
  lower.default <- list(coef = -Inf,
                        angle = 0,
                        disp = 0,
                        lambda = 0)
  lower.values <- getBounds(model, phy, X, y, number.params, lower.bound, lower.default)
  lower.values <- transform_values(lower.values)
  # upper
  upper.default <- list(coef = Inf,
                        angle = pi,
                        disp = Inf,
                        lambda = maxLambda(phy))
  upper.values <- getBounds(model, phy, X, y, number.params, upper.bound, upper.default)
  upper.values <- transform_values(upper.values)
  # param names
  param_names <- getParamNames(number.params)
  
  # root tip for reml
  rootTip <- NULL
  if (method == "reml") {
    rootTip <- which.min(colSums(cophenetic.phylo(phy)))
  }

  ## Actual fit
  minus_like <- switch(method,
                       reml = minusLikelihoodREML(ncol(y)),
                       fixed.root = minusLikelihoodFixedRoot(X, ncol(y)),
                       random.root = minusLikelihoodRandomRoot(ncol(y)))
  
  res <- fit_function(minus_like, phy, y, X,  model, start.values, lower.values, upper.values, optim = optim, rootTip = rootTip)
  res$rootTip <- rootTip
  
  return(res)
}

transform_values <- function(param) {
  # disp
  disp_params <- grepl("disp", names(param))
  param[disp_params] <- log_sorted_transform(param[disp_params])
  # lambda
  lambda_params <- grepl("lambda", names(param))
  param[lambda_params] <- log(param[lambda_params])
  # angle
  angle_params <- grepl("angle", names(param))
  param[angle_params] <- logit_transform(param[angle_params], 0, pi)
  return(param)
}

back_transform_values <- function(param) {
  # disp
  disp_params <- grepl("disp", names(param))
  param[disp_params] <- log_sorted_back_transform(param[disp_params])
  # lambda
  lambda_params <- grepl("lambda", names(param))
  param[lambda_params] <- exp(param[lambda_params])
  # angle
  angle_params <- grepl("angle", names(param))
  param[angle_params] <- logit_back_transform(param[angle_params], 0, pi)
  return(param)
}

logit_transform <- function(x, lower_bound, upper_bound) {
  x <- (x - lower_bound)/(upper_bound - lower_bound)
  return(log(x/(1-x)))
}

logit_back_transform <- function(x, lower_bound, upper_bound) {
  x <- 1/(1 + exp(-x))
  return(lower_bound + (upper_bound - lower_bound) * x)
}

log_sorted_transform <- function(x) {
  if (all(x == 0)) return(rep(-Inf, length(x)))
  if (all(x == Inf)) return(rep(Inf, length(x)))
  y <- log(x)
  for (i in seq_along(y)[-1]) {
    y[i] <- log(y[i-1] - y[i])
  }
  return(y)
}

log_sorted_back_transform <- function(x) {
  if (all(x == -Inf)) return(rep(0, length(x)))
  if (all(x == Inf)) return(rep(Inf, length(x)))
  y <- x
  for (i in seq_along(y)[-1]) {
    y[i] <- y[i-1] - exp(x[i])
  }
  return(exp(y))
}

# transform_hessian <- function(hessian, param) {
#   param["disp"] <- log(param["disp"])
#   if (!is.na(param["lambda"])) param["lambda"] <- log(param["lambda"])
#   return(param)
# }

#' @title Maximum Likelihood estimator for a Cauchy model
#'
#' @description
#' Find the maximum likelihood.
#'
#' @inheritParams fitCauchy
#' @param minus_like a function giving the minus likelihood of the model.
#' @param rootTip root tip for the reml.
#' 
#' @return A list, with the maximum likelihood rate parameter, and the likelihood value.
#' 
#' @keywords internal
#'
fit_function <- function(minus_like, tree, trait, X, model, start.values, lower.values, upper.values, optim, rootTip) {
  param_names <- names(start.values)
  # opt
  opt <- do_optim(minus_like, start.values, lower.values, upper.values, optim = optim,
                  param_names = param_names, tree = tree, trait = trait, Xdesign = X, model = model, rootTip = rootTip)
  # result
  sol <- opt$solution
  names(sol) <- param_names
  sol <- back_transform_values(sol)
  # parameters
  return(list(param = sol,
              logLikelihood = -opt$objective))
}

#' @title Minus Likelihood function for a Cauchy model
#'
#' @description
#' Gives the minus likelihood function, with fixed root.
#'
#' @param Xdesign the design matrix
#' @param pdim dimension of the process
#' 
#' @return A function
#' 
#' @keywords internal
#'
minusLikelihoodFixedRoot <- function(Xdesign, pdim) {
  if (pdim == 1) {
    if (ncol(Xdesign) == 1 && sum(Xdesign) == nrow(Xdesign)) {
      return(minusLikelihoodFixedRoot_mu)
    }
    return(minusLikelihoodFixedRoot_lm)
  } else if (pdim == 2) {
    return(minusLikelihoodFixedRoot_bi)
  } else {
    stop("Dimension above 2 not yet implemented.")
  }
}

#' @title Minus Likelihood function for a Cauchy model
#'
#' @description
#' Gives the minus likelihood function, with fixed root.
#'
#' @inheritParams fitCauchy.internal
#' @param param the parameters where to evaluate the function
#' @param param_names the parameter names
#' @param Xdesign the design matrix
#' 
#' @return The value of minus the log-likelihood.
#' 
#' @keywords internal
#'
minusLikelihoodFixedRoot_mu <- function(param, param_names, tree, trait, Xdesign, model, rootTip) {
  names(param) <- param_names
  param <- back_transform_values(param)
  phy_trans <- transformBranchLengths(tree, model, param)
  # tree_height <- max(node.depth.edgelength(phy_trans))
  # return(length(phy_trans$tip.label) * log(param["disp"] * 100 / tree_height) - logDensityTipsCauchy(phy_trans, trait / param["disp"] / 100 * tree_height, param[1], 0.01 * tree_height, method = "fixed.root"))  
  return(-logDensityTipsCauchy(phy_trans, trait, param[1], param["disp"], method = "fixed.root", do_checks = FALSE))
}

#' @title Minus Likelihood function for a Cauchy model
#'
#' @description
#' Gives the minus likelihood function, with fixed root.
#'
#' @inheritParams minusLikelihoodFixedRoot_mu
#' 
#' @return The value of minus the log-likelihood.
#' 
#' @keywords internal
#'
minusLikelihoodFixedRoot_bi <- function(param, param_names, tree, trait, Xdesign, model, rootTip) {
  names(param) <- param_names
  param <- back_transform_values(param)
  phy_trans <- transformBranchLengths(tree, model, param)
  # tree_height <- max(node.depth.edgelength(phy_trans))
  # return(length(phy_trans$tip.label) * log(param["disp"] * 100 / tree_height) - logDensityTipsCauchy(phy_trans, trait / param["disp"] / 100 * tree_height, param[1], 0.01 * tree_height, method = "fixed.root"))  
  return(-logDensityTipsCauchyBi(phy_trans, trait,
                                 root.value = param[grepl("coef", names(param))],
                                 disp = param[grepl("disp", names(param))],
                                 angle = param[grepl("angle", names(param))],
                                 method = "fixed.root"))
}

#' @title Minus Likelihood function for a Cauchy model
#'
#' @description
#' Gives the minus likelihood function, with fixed root, lm model
#'
#' @inheritParams minusLikelihoodFixedRoot_mu
#' 
#' @return The value of minus the log-likelihood.
#' 
#' @keywords internal
#'
minusLikelihoodFixedRoot_lm <- function(param, param_names, tree, trait, Xdesign, model, rootTip) {
  names(param) <- param_names
  param <- back_transform_values(param)
  centralTips <- drop(Xdesign %*% param[grepl("coef", names(param))])
  phy_trans <- transformBranchLengths(tree, model, param)
  return(-logDensityTipsCauchy(phy_trans, trait - centralTips, 0, param["disp"], method = "fixed.root", do_checks = FALSE))
}

#' @title Minus Likelihood function for a Cauchy model
#'
#' @description
#' Gives the minus likelihood function, with random root.
#'
#' @inheritParams minusLikelihoodFixedRoot
#' 
#' @return The value of minus the log-likelihood.
#' 
#' @keywords internal
#'
minusLikelihoodRandomRoot <- function(pdim) {
  if (pdim == 1) {
    return(minusLikelihoodRandomRoot_uni)
  } else if (pdim == 2) {
    return(minusLikelihoodRandomRoot_bi)
  } else {
    stop("Dimension above 2 not yet implemented.")
  }
}

#' @title Minus Likelihood function for a Cauchy model
#'
#' @description
#' Gives the minus likelihood function, with random root.
#'
#' @inheritParams minusLikelihoodFixedRoot_mu
#' 
#' @return The value of minus the log-likelihood.
#' 
#' @keywords internal
#'
minusLikelihoodRandomRoot_uni <- function(param, param_names, tree, trait, Xdesign, model, rootTip) {
  names(param) <- param_names
  param <- back_transform_values(param)
  phy_trans <- transformBranchLengths(tree, model, param)
  return(-logDensityTipsCauchy(phy_trans, trait, root.value = 0.0, disp = param["disp"], method = "random.root", do_checks = FALSE))
}

#' @title Minus Likelihood function for a Cauchy model
#'
#' @description
#' Gives the minus likelihood function, with random root.
#'
#' @inheritParams minusLikelihoodFixedRoot_mu
#' 
#' @return The value of minus the log-likelihood.
#' 
#' @keywords internal
#'
minusLikelihoodRandomRoot_bi <- function(param, param_names, tree, trait, Xdesign, model, rootTip) {
  names(param) <- param_names
  param <- back_transform_values(param)
  phy_trans <- transformBranchLengths(tree, model, param)
  return(-logDensityTipsCauchyBi(phy_trans, trait,
                                 root.value = c(0.0, 0.0),
                                 disp = param[grepl("disp", names(param))],
                                 angle = param[grepl("angle", names(param))],
                                 method = "random.root"))
}

#' @title Minus REML function for a Cauchy model
#'
#' @description
#' Gives the minus REML function.
#'
#' @inheritParams minusLikelihoodFixedRoot
#' 
#' @return The value of minus the log-REML.
#' 
#' @keywords internal
#'
minusLikelihoodREML <- function(pdim) {
  if (pdim == 1) {
    return(minusLikelihoodREML_uni)
  } else if (pdim == 2) {
    return(minusLikelihoodREML_bi)
  } else {
    stop("Dimension above 2 not yet implemented.")
  }
}

#' @title Minus REML function for a Cauchy model
#'
#' @description
#' Gives the minus REML function.
#'
#' @inheritParams minusLikelihoodFixedRoot_mu
#' 
#' @return The value of minus the log-REML.
#' 
#' @keywords internal
#'
minusLikelihoodREML_uni <- function(param, param_names, tree, trait, Xdesign, model, rootTip) {
  names(param) <- param_names
  param <- back_transform_values(param)
  phy_trans <- transformBranchLengths(tree, model, param)
  # tree_height <- max(node.depth.edgelength(phy_trans))
  # return((length(phy_trans$tip.label) - 1) * log(param["disp"] * 1000 / tree_height) - logDensityTipsCauchy(tree = phy_trans, tipTrait = trait / param["disp"] / 1000 * tree_height, root.value = NULL, disp = 0.001 * tree_height, method = "reml", rootTip = rootTip))  
  return(-logDensityTipsCauchy(tree = phy_trans, tipTrait = trait, root.value = NULL, disp = param["disp"], method = "reml", rootTip = rootTip, do_checks = FALSE))
}

#' @title Minus REML function for a Cauchy model
#'
#' @description
#' Gives the minus REML function.
#'
#' @inheritParams minusLikelihoodFixedRoot_bi
#' 
#' @return The value of minus the log-REML.
#' 
#' @keywords internal
#'
minusLikelihoodREML_bi <- function(param, param_names, tree, trait, Xdesign, model, rootTip) {
  names(param) <- param_names
  param <- back_transform_values(param)
  phy_trans <- transformBranchLengths(tree, model, param)
  return(-logDensityTipsCauchyBi(phy_trans, trait,
                                 root.value = NULL,
                                 disp = param[grepl("disp", names(param))],
                                 angle = param[grepl("angle", names(param))],
                                 method = "reml",
                                 rootTip = rootTip))
}

#' @title Initialization of the position parameter.
#'
#' @description
#' Initialize using the trimmed mean of the trait.
#'
#' @inheritParams fitCauchy
#'
#' @return The initial position parameter.
#' 
#' @references
#' Rothenberg T. J., Fisher F. M., Tilanus C. B. 1964. A Note on Estimation from a Cauchy Sample. Journal of the American Statistical Association. 59:460–463.
#' Rousseeuw P.J., Croux C. 1993. Alternatives to the Median Absolute Deviation. Journal of the American Statistical Association. 88:1273–1283.
#' 
#' @keywords internal
#'
initPositionParameter <- function(trait) {
  return(apply(trait, 2, function(tt) mean(tt, trim = 0.76)))
}

#' @title Initialization of the dispersion parameter.
#'
#' @description
#' Initialize the dispersion parameter.
#' 
#' @param center the center parameter of the distribution
#' @param method.init.disp the robust estimator method
#' 
#' @details 
#' Constants are taken from Rousseeuw & Croux, 1993.
#' 
#' @references
#' Rousseeuw P.J., Croux C. 1993. Alternatives to the Median Absolute Deviation. Journal of the American Statistical Association. 88:1273–1283.
#' 
#' @inheritParams fitCauchy
#'
#' @return The initial dispersion parameter.
#' 
#' @keywords internal
#'
initDispersionParameter <- function(tree, trait, center, method.init.disp = c("Qn", "Sn", "MAD", "IQR"), method) {
  if (is.null(center)) center <- initPositionParameter(trait)
  # root_edge <- tree$root.edge
  # if (is.null(root_edge)) root_edge <- 0.0
  # root_edge <- 0.0
  norm_trait <- (trait - center) / diag(ape::vcv(tree))
  method.init.disp <- match.arg(method.init.disp)
  res <- switch (method.init.disp,
                 IQR = 0.5 * IQR(norm_trait),
                 MAD = mad(norm_trait, center = 0, constant = 1),
                 Sn = robustbase::Sn(norm_trait, constant = 0.7071),
                 Qn = robustbase::Qn(norm_trait, constant = 1.2071)
  )
  return(res)
}

#' @title Initialization of the lambda parameter.
#'
#' @description
#' Initialize the lambda parameter.
#' 
#' @param disp_hat the previously estimated dispersion
#' @param tol the percentage of tip pairs to keep for the computation. Default to 0.1.
#' @inheritParams fitCauchy
#' 
#' @details 
#' Use linear combinations of cherries to get a first estimate of lambda.
#' 
#' @return The initial lambda parameter.
#' 
#' @keywords internal
#'
initLambdaParameter <- function(tree, trait, disp_hat, tol = 0.1) {
  ## Tree distance matrices
  dist_phylo <- cophenetic.phylo(tree)
  vcv_tree <- vcv(tree)
  tree_heights <- diag(vcv_tree)
  ## Keep only tol % closest tip paris
  thershold <- quantile(dist_phylo[upper.tri(dist_phylo, diag = FALSE)], tol)
  tip_pairs <- which(upper.tri(dist_phylo) & dist_phylo <= thershold, arr.ind = TRUE)
  ## Normalize tip pairs differences
  diff_tips <- abs(trait[tip_pairs[, 1]] - trait[tip_pairs[, 2]])
  diff_tips <- diff_tips - disp_hat * (tree_heights[tip_pairs[, 1]] + tree_heights[tip_pairs[, 2]])
  diff_tips <- diff_tips / (- 2 * disp_hat * apply(tip_pairs, 1, function(tt) vcv_tree[tt[1], tt[2]]))
  ## Estimation
  lambda_estim <- median(diff_tips)
  lambda_estim <- abs(lambda_estim)
  ## Bounds
  max_lambda <- maxLambda(tree)
  # if (lambda_estim < 0) return(0)
  if (lambda_estim > max_lambda) return(max_lambda)
  return(lambda_estim)
}

#' @importFrom phylolm transf.branch.lengths
NULL

#' @title Transform branch lengths
#'
#' @description
#' Transform branch lengths for pagel lambda model
#'
#' @param model the model
#' @param phy the phylogenetic tree
#' @param param the parameters
#'
#' @keywords internal
#'
transformBranchLengths <- function(phy, model, param) {
  phy_trans <- phy
  phy_trans$root.edge <- NULL
  if (model == "lambda") phy_trans <- phylolm::transf.branch.lengths(phy_trans, model = model, parameters = list(lambda = param["lambda"]))$tree
  phy_trans$root.edge <- phy$root.edge
  return(phy_trans)
}

#' @title Maximum lambda value
#'
#' @description
#' Find maximum lambda value.
#' Function taken from \code{phytools:::maxLambda}.
#'
#' @param phy the phylogenetic tree
#'
#' @keywords internal
#'
maxLambda <- function (phy) {
  if (!inherits(phy, "phylo")) stop("tree should be an object of class \"phylo\".")
  if (is.ultrametric(phy)) {
    h <- ape::node.depth.edgelength(phy)
    return(max(h[phy$edge[, 2]]) / max(h[phy$edge[, 1]]))
  }
  else return(1)
}

#' @title Get starting values
#'
#' @description
#' Get starting values given the model.
#'
#' @param model the model
#' @param phy the phylogenetic tree
#' @param X model matrix
#' @param y the response vector
#' @param starting.value the input starting values
#'
#' @keywords internal
#'
getStartingValues <- function(model, phy, X, y, number.params, starting.value, method.init.disp, method) {
  tmp_fun <- switch(model,
                    cauchy = getStartingValuesCauchy,
                    lambda = getStartingValuesLambda
  )
  ss <- tmp_fun(phy, X, y, number.params, starting.value, method.init.disp, method)
  names(ss) <- getParamNames(number.params)
  return(ss)
}

#' @title Get starting values for a Cauchy
#'
#' @description
#' Get starting values for a Cauchy process
#'
#' @inheritParams getStartingValues
#'
#' @keywords internal
#'
getStartingValuesCauchy <- function(phy, X, y, number.params, starting.value, method.init.disp, method) {
  # starting values
  if (!is.null(X)) {
    # intercept only
    if (ncol(X) == 1 && all.equal(sum(X), nrow(X))) {
      if (is.null(starting.value$x0)) {
        start.coef <- initPositionParameter(y)
      } else {
        start.coef <- starting.value$x0
        if (!check_numeric_vector(start.coef, y)) stop("Starting value for x0 should be a vector of real numbers with dimension equal to the number of traits.")
      }
    } else {
      start.coef <- robustbase::lmrob.S(X, y, control = robustbase::lmrob.control())$coefficients
    }
  } else {
    start.coef <- NULL
  }
  center <- start.coef
  if (is.null(center)) center <- initPositionParameter(y)
  # angles
  start.angle <- NULL
  if (ncol(y) == 2) {
    if (is.null(starting.value$disp)) {
      svd_init <- svd(t(t(y) - center))
      theta_init_1 <- atan2(svd_init$v[2, 1], svd_init$v[1, 1])
      if (theta_init_1 < 0) theta_init_1 <- theta_init_1 + pi
      theta_init_2 <- (theta_init_1 + pi / 2) %% pi
      start.angle <- c(theta_init_1, theta_init_2)
    } else {
      start.angle <- starting.value$angle
      if (!(check_numeric_vector(start.angle, y) && all(start.angle >= 0) && all(start.angle < pi) && all(diff(sort(start.angle)) > 0.001))) stop("Starting value for the directions should be a vector of angles between 0 and pi, all differrent, with dimension equal to the number of traits.")
    }
  } 
  # disp
  if (is.null(starting.value$disp)) {
    if (ncol(y) == 1) {
      start.disp <- initDispersionParameter(phy, y, start.coef, method.init.disp = method.init.disp, method = method)
    }
    if (ncol(y) == 2) {
      rot_mat <- sapply(start.angle, function(x) c(cos(x), sin(x)))
      dat_trans <- y %*% rot_mat
      mu_trans <- t(rot_mat) %*% center
      start.disp <- c(initDispersionParameter(phy, dat_trans[, 1], mu_trans[1,1], method.init.disp = method.init.disp, method = method),
                      initDispersionParameter(phy, dat_trans[, 2], mu_trans[2,1], method.init.disp = method.init.disp, method = method))
      oo <- order(start.disp, decreasing = TRUE)
      start.disp <- start.disp[oo]
      start.angle <- start.angle[oo]
    }
  } else {
    start.disp <- starting.value$disp
    if (!(check_numeric_vector(start.disp, y) && all(start.disp > 0) && all(diff(start.disp) < -0.001))) stop("Starting value for the dispersion should be a vector of positive real numbers, all different and sorted in a decreasing order, with dimension equal to the number of traits.")
  }

  start.values <- c(start.coef, start.angle, start.disp)
  
  return(start.values)
}

check_numeric_vector <- function(v, y) {
  return(is.null(dim(v)) & length(v) == ncol(y) & is.numeric(v))
}

#' @title Get starting values for a Cauchy Lambda
#'
#' @description
#' Get starting values for a Cauchy Lambda process
#'
#' @inheritParams getStartingValues
#'
#' @keywords internal
#'
getStartingValuesLambda <- function(phy, X, y, number.params, starting.value, method.init.disp, method) {
  ## Cauchy parameters
  start.cauchy <- getStartingValuesCauchy(phy, X, y, number.params, starting.value, method.init.disp, method)
  disp_hat <- ifelse(length(start.cauchy) == 2, start.cauchy[2], start.cauchy[1])
  ## Lambda 
  if (is.null(starting.value$lambda)) {
    start.lambda <- initLambdaParameter(phy, y, disp_hat)
  } else {
    start.lambda <- starting.value$lambda
    if (!(is.null(dim(start.lambda)) && length(start.lambda) == 1 && is.numeric(start.lambda) && start.lambda >= 0 && start.lambda <= 1)) stop("Starting value for the lambda parameter should be between 0 and 1.")
  }
  return(c(start.cauchy, start.lambda))
}

#' @title Get bounds
#'
#' @description
#' Get bounds given the model.
#'
#' @param model the model
#' @param phy the phylogenetic tree
#' @param X model matrix
#' @param y the response vector
#' @param values the input values for the bound
#' @param default.values the default values for the bound
#'
#' @keywords internal
#'
getBounds <- function(model, phy, X, y, number.params, values, default.values) {
  tmp_fun <- switch(model,
                    cauchy = getBoundsCauchy,
                    lambda = getBoundsLambda
  )
  ss <- tmp_fun(phy, X, y, number.params, values, default.values)
  names(ss) <- getParamNames(number.params)
  return(ss)
}

#' @title Get bound for a Cauchy process
#'
#' @description
#' Get bounds for a Cauchy process
#'
#' @inheritParams getBounds
#'
#' @keywords internal
#'
getBoundsCauchy <- function(phy, X, y, number.params, values, default.values) {
  if (is.null(values$disp)) {
    bound.disp <- rep(default.values$disp, number.params$disp)
  } else {
    bound.disp <- values$disp
    if (!(check_numeric_vector(bound.disp, y) && all(bound.disp >= 0))) stop("Upper and lower values for the dispersion should be positive real numbers.")
  }
  if (is.null(values$angle)) {
    bound.angle <- rep(default.values$angle, number.params$angle)
  } else {
    bound.angle <- values$angle
    if (!(check_numeric_vector(bound.angle, y) && all(bound.angle >= 0) && all(bound.angle < pi))) stop("Upper and lower values for the angles should be between 0 and pi.")
  }
  if (number.params$mean > 0) {
    bound.coef <- rep(default.values$coef, number.params$mean)
  } else {
    bound.coef <- NULL
  }
  bound.values <- c(bound.coef, bound.angle, bound.disp)
  return(bound.values)
}

#' @title Get bounds for a Cauchy Lambda
#'
#' @description
#' Get bounds for a Cauchy Lambda process
#'
#' @inheritParams getBounds
#'
#' @keywords internal
#'
getBoundsLambda <- function(phy, X, y, number.params, values, default.values) {
  if (is.null(values$lambda)) {
    bound.lambda <- default.values$lambda
  } else {
    bound.lambda <- values$lambda
    if (!(is.null(dim(bound.lambda)) && length(bound.lambda) == 1 && is.numeric(bound.lambda) && bound.lambda >= 0 && bound.lambda <= 1)) stop("Bounds for the lambda parameter should be between 0 and 1.")
  }
  return(c(getBoundsCauchy(phy, X, y, number.params, values, default.values), bound.lambda))
}

#' @title Get parameter names
#'
#' @description
#' Get the names of the parameters depending on the model, for a univariate trait.
#'
#' @param model model
#' @param X model matrix
#'
#' @keywords internal
#'
getParamNamesUni <- function(model, X) {
  pnames <- list(mean = ifelse(is.null(X), 0, ncol(X)),
                 angle = 0,
                 disp = 1,
                 lambda = ifelse(model == "lambda", 1, 0))
  return(getParamNames(pnames))
}

#' @title Get parameter names
#'
#' @description
#' Get the names of the parameters depending on the model.
#' @param number.params named list with number of parameters.
# param coef number of regression or mean coefficients
# param angle number of direction coefficients
# param disp number of dispersion coefficients
# param lambda number of lambda coefficients
#'
#' @keywords internal
#'
getParamNames <- function(number.params) {
  # Assumes that number.params$coef = 1 or number.params$mean = 1
  return(c(getStandardName("coef", number.params$mean),
           getStandardName("angle", number.params$angle),
           getStandardName("disp", number.params$disp),
           getStandardName("lambda", number.params$lambda)))
}

getStandardName <- function(name, n) {
  if (n == 0) return(NULL)
  if (n == 1) return(name)
  if (n > 1) return(paste0(name, 1:n))
}

#' @title Check Matrix Parameter
#'
#' @description
#' Check that the parameters are compatible with the tree. Throws an error if not.
#'
#' @param trait (named) vector or matrix of traits being tested.
#' @param tree A phylogenetic tree.
#' @param name name of the trait. Default to 'trait'.
#'
#' @keywords internal
#'
checkTraitTree <- function(trait, tree, name = "trait") {
  N <- length(tree$tip.label)
  if (is.null(dim(trait))) { # trait is a vector
    if (length(trait) != N) {
      stop(paste0("`", name, "` should have the same length as the number of taxa in the tree."))
    }
    if ((is.null(tree$tip.label) || is.null(names(trait)))){
      stop(paste0("`", name, "` and/or the tips of the phylogeny are not named. I could not check for consistency. Please give names consistent names to the tree tip labels and the row names of matrix `", name, "` to avoid any ambiguity."))
    } else {
      if (!all(tree$tip.label == names(trait))){
        # Match
        tree_data_cor <- match(tree$tip.label, names(trait))
        data_tree_cor <- match(names(trait), tree$tip.label)
        if (anyNA(tree_data_cor)) {
          # Species in the tree NOT in data
          stopMessage <- paste0("Species '", paste(tree$tip.label[is.na(tree_data_cor)], collapse = "', '"), "' are in the tree but not in ", name, ".")
          if (anyNA(data_tree_cor)) {
            # Species in data NOT in the tree
            stop(stopMessage, "\n  ", "Species '", paste(names(trait)[is.na(data_tree_cor)], collapse = "', '"), "' are in ", name, " but not in the tree.")
            
          }
        }
        if (length(unique(tree_data_cor)) != length(tree$tip.label)){
          stop(paste0("`", name, "` names do not match the tip labels."))
        }
        warning(paste0("`", name, "` was not sorted in the correct order, when compared with the tips label. I am re-ordering it."))
        trait <- trait[tree_data_cor]
      }
    }
  } else { # trait is a matrix
    if (nrow(trait) != N) {
      stop(paste0("`", name, "` should have as many rows as the number of taxa in the tree."))
    }
    if ((is.null(tree$tip.label) || is.null(rownames(trait)))){
      stop(paste0("`", name, "` and/or the tips of the phylogeny are not named. I could not check for consistency. Please give names consistent names to the tree tip labels and the row names of matrix `", name, "` to avoid any ambiguity."))
    } else {
      if (!all(tree$tip.label == rownames(trait))){
        # Match
        tree_data_cor <- match(tree$tip.label, rownames(trait))
        data_tree_cor <- match(rownames(trait), tree$tip.label)
        if (anyNA(tree_data_cor)) {
          # Species in the tree NOT in data
          stopMessage <- paste0("Species '", paste(tree$tip.label[is.na(tree_data_cor)], collapse = "', '"), "' are in the tree but not in ", name, ".")
          if (anyNA(data_tree_cor)) {
            # Species in data NOT in the tree
            stop(stopMessage, "\n  ", "Species '", paste(rownames(trait)[is.na(data_tree_cor)], collapse = "', '"), "' are in ", name, " but not in the tree.")
            
          }
        }
        if (length(unique(tree_data_cor)) != length(tree$tip.label)){
          stop(paste0("`", name, "` names do not match the tip labels."))
        }
        warning(paste0("`", name, "` was not sorted in the correct order, when compared with the tips label. I am re-ordering it."))
        trait <- trait[tree_data_cor, , drop = FALSE]
      }
    }
  }
  return(trait)
}

#' @title NLOPT optimization
#'
#' @description
#' Perform the optimization
#'
#' @param minus_like the function to be minimized
#' @param start.values vector of starting values
#' @param lower.values vector of lower values
#' @param upper.values vector of upper values
#' @param options_nloptr options to be passed to \code{nloptr}
#' @param optim should a global optimization be performed first ?
#' @param ... further arguments to be passed to minus_like
#'
#' @keywords internal
#'
do_optim <- function(minus_like, start.values, lower.values, upper.values,
                     options_nloptr = list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1e-8, "maxeval" = 1000),
                     optim = c("local", "global"), ...) {
  init_local <- start.values
  optim <- match.arg(optim)
  global <- optim == "global"
  if (global) {
    local_opts <- list(algorithm = "NLOPT_LN_BOBYQA", xtol_rel = 1e-04, ftol_rel = 1e-5)
    opts_nlopt <- list("algorithm" = "NLOPT_GN_MLSL_LDS", "xtol_rel" = 1e-02, "maxeval" = 100)
    opts_nlopt[["local_opts"]] <- local_opts
    box_values <- abs(start.values) * 5
    opt_global <- nloptr::nloptr(x0 = start.values, eval_f = minus_like,
                                 lb = pmax(start.values - box_values, lower.values),
                                 ub = pmin(start.values + box_values, upper.values),
                                 opts = opts_nlopt, ...)
    init_local <- opt_global$solution
  }
  opt <- nloptr::nloptr(x0 = init_local, eval_f = minus_like,
                        lb = lower.values, ub = upper.values,
                        opts = options_nloptr, ...)
  # local around init
  # opt_fun <- function(algo, xtol_rel, x0) {
  #   opt <- nloptr::nloptr(x0 = x0, eval_f = minus_like,
  #                         lb = lower.values, ub = upper.values,
  #                         opts = list("algorithm" = algo, "xtol_rel" = xtol_rel, "maxeval" = 1000))
  #   return(opt)
  # }
  # all_opts <- lapply(optim.algo, opt_fun, xtol_rel = 1e-3, x0 = start.values)
  # best_algo <- which.min(sapply(all_opts, function(opt) opt$objective))
  # sol <- all_opts[[best_algo]]$solution
  # names(sol) <- param_names
  # opt <- opt_fun(optim.algo[[best_algo]], xtol_rel = 1e-8, x0 = sol)
  # opt <- opt_fun("NLOPT_LN_BOBYQA", xtol_rel = 1e-08, x0 = opt_global$solution)
}

#' @title Check For Duplicated Entries
#'
#' @description
#' Check that the trait are compatible with the tree. Throws an error if not.
#' Assumes that the trait and tree tips are in the same order, using function 
#' \code{\link{checkTraitTree}}.
#'
#' @param trait (named) vector or matrix of traits being tested.
#' @param tree phylogenetic tree.
#'
#' @keywords internal
#'
checkDuplicates <- function(trait, tree) {
  if (anyDuplicated(trait + diag(vcv(tree)))) stop("The trait vector has duplicated entries on tips that are equidistant from the root. The algorithm cannot currently handle this case. Please consider adding some noise to the tip trait values.")
}

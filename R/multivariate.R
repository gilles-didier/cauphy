#' @title Log Density of a Multivariate Cauchy Process
#'
#' @description
#' Compute the log density of the vector of trait at the tips of the phylogenetic tree, 
#' assuming a Cauchy process.
#' 
#' @param tree a phylogenetic tree of class \code{\link[ape]{phylo}}.
#' @param tipTrait a named matrix of tip trait values, with row names matching the tree labels.
#' @param rootTip the tip used to re-root the tree, when the REML method is used.
#' If \code{NULL}, the tip with the smallest distance to all the others is used (see Details).
#' Ignored in \code{method != "reml"}.
#' @param root.value the root starting value of the process.
#' @param transMat the transformation matrix.
#' @param method the method used to compute the likelihood.
#' One of \code{reml} (the default), \code{fixed.root} or \code{random.root}.
#' 
#' @return the log density value.
#' 
#' 
#' @keywords internal
#' 
#' 
logDensityTipsCauchyMulti <- function(tree, tipTrait, root.value = NULL, transMat, method = c("reml", "random.root", "fixed.root"), rootTip = NULL) {
  # dim
  tipTrait <- checkTraitTree(tipTrait, tree)
  p <- ncol(tipTrait)
  n <- nrow(tipTrait)
  if (method == "reml") n <- n - 1
  if (!all.equal(dim(transMat), c(p, p))) stop("'transMat' must have the same dimension as the trait.")
  # transform trait
  M <- rep(0, p)
  if (!is.null(root.value)) M <- root.value
  transTrait <- solve(transMat, t(tipTrait) - M)
  root.trans <- NULL
  if (!is.null(root.value)) root.trans <- 0.0
  # Independent likelihoods
  ind_ll <- apply(transTrait, 1, function(tt) logDensityTipsCauchy(tree = tree, tipTrait = tt,
                                                                   root.value = root.trans, disp = 1.0,
                                                                   method = method, rootTip = rootTip))
  # likelihood
  ll <- - n * determinant(transMat, logarithm = TRUE)$modulus + sum(ind_ll)
  return(as.vector(ll))
}

#' @title Log Density of a Bivariate Cauchy Process
#'
#' @description
#' Compute the log density of the vector of trait at the tips of the phylogenetic tree, 
#' assuming a Cauchy process.
#' 
#' @param tree a phylogenetic tree of class \code{\link[ape]{phylo}}.
#' @param tipTrait a named matrix of tip trait values, with row names matching the tree labels.
#' @param rootTip the tip used to re-root the tree, when the REML method is used.
#' If \code{NULL}, the tip with the smallest distance to all the others is used (see Details).
#' Ignored in \code{method != "reml"}.
#' @param root.value the root starting value of the process.
#' @param angle vector of the 2 angles of the 2 axes.
#' @param disp vector of the 2 dispersion values along the two axes.
#' @param method the method used to compute the likelihood.
#' One of \code{reml} (the default), \code{fixed.root} or \code{random.root}.
#' 
#' 
#' @return the log density value.
#' 
#' @seealso \code{\link{fitCauchy}}
#' 
#' @keywords internal
#' 
#' 
logDensityTipsCauchyBi <- function(tree, tipTrait,
                                   root.value = NULL, disp, angle,
                                   method = c("reml", "random.root", "fixed.root"), rootTip = NULL) {
  # dim
  tipTrait <- checkTraitTree(tipTrait, tree)
  p <- ncol(tipTrait)
  n <- nrow(tipTrait)
  if (method == "reml") n <- n - 1
  if (p != 2) stop("For the bivariate Cauchy, 'tipTrait' must be a matrix with 2 columns.")
  if (!is.null(root.value) && length(root.value) != 2) stop("For the bivariate Cauchy, if specified 'root.value' must a vector of length 2.")
  if (length(angle) != 2) stop("For the bivariate Cauchy, 'angle' must be a vector of length 2.")
  if (length(disp) != 2) stop("For the bivariate Cauchy, 'disp' must be a vector of length 2.")
  # transform trait
  basisMatr <- sapply(angle, function(x) c(cos(x), sin(x)))
  transTrait <- solve(basisMatr, t(tipTrait))
  if (!is.null(root.value)) root.value <- solve(basisMatr, root.value)
  # Independent likelihoods
  ind_ll <- sapply(1:2, function(k) logDensityTipsCauchy(tree = tree, tipTrait = transTrait[k, ],
                                                         root.value = root.value[k], disp = disp[k],
                                                         method = method, rootTip = rootTip))
  # likelihood
  ll <- - n * determinant(basisMatr, logarithm = TRUE)$modulus + sum(ind_ll)
  return(as.vector(ll))
}

#' @title Log Density of a Cauchy Process
#'
#' @description
#' Compute the log density of the vector of trait at the tips of the phylogenetic tree, 
#' assuming a Cauchy process.
#' 
#' @param tree a phylogenetic tree of class \code{\link[ape]{phylo}}.
#' @param tipTrait a names vector of tip trait values, with names matching the tree labels.
#' @param rootTip the tip used to re-root the tree, when the REML method is used.
#' If \code{NULL}, the tip with the smallest distance to all the others is used (see Details).
#' Ignored in \code{method != "reml"}.
#' @param root.value the root starting value of the process.
#' @param disp the dispersion value.
#' @param method the method used to compute the likelihood.
#' One of \code{reml} (the default), \code{fixed.root} or \code{random.root}.
#' See Details.
#' 
#' @details
#' The parameters of the Cauchy Process (CP)
#' are \code{disp}, the dispersion of the process,
#' and \code{root.value}, the starting value of the process at the root (for \code{method="fixed.root"}).
#' 
#' The model assumes that each increment of the trait \eqn{X} on a branch going from node \eqn{k} to \eqn{l} 
#' follows a Cauchy distribution, with a dispersion proportional to the length \eqn{t_l} of the branch:
#' \deqn{X_l - X_k \sim \mathcal{C}(0, \text{disp} \times t_l).}
#' 
#' The \code{method} argument specifies the type of likelihood that is computed:
#' \itemize{
#'   \item{\code{method="reml"}:}{ 
#'   the dispersion parameter is fitted using the REML criterion,
#'   obtained by re-rooting the tree to one of the tips.
#'   The default tip used to reroot the tree is:
#'   \code{rootTip = which.min(colSums(cophenetic.phylo(tree)))}.
#'   Any tip can be used, but this default empirically proved to be the most robust numerically;
#'   }
#'   \item{\code{method="random.root"}:}{ 
#'   the root value is assumed to be a random Cauchy variable, centered at \code{root.value=0},
#'   and with a dispersion \code{disp_root = disp * root.edge};
#'   }
#'   \item{\code{method="fixed.root"}:}{ 
#'   the model is fitted conditionally on the root value \code{root.value},
#'   i.e. with a model where the root value is fixed and inferred from the data.
#'   }
#' }
#' 
#' @return the log density value.
#' 
#' @seealso \code{\link{fitCauchy}}
#' 
#' @examples
#' phy <- ape::rphylo(5, 0.1, 0)
#' dat <- rTraitCauchy(n = 1, phy = phy, model = "cauchy", parameters = list(root.value = 0, disp = 1))
#' logDensityTipsCauchy(phy, dat, 0, 1, method = "fixed.root")
#' 
#' @author Paul Bastide \email{paul.bastide@m4x.org} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @export
#' 
#' 
logDensityTipsCauchy <- function(tree, tipTrait, root.value = NULL, disp, method = c("reml", "random.root", "fixed.root"), rootTip = NULL) {
  # type
  method <- match.arg(method)
  type <- switch(method,
                 reml = 2,
                 random.root = 0,
                 fixed.root = 1)
  # checks
  if (!inherits(tree, "phylo")) stop("object \"tree\" is not of class \"phylo\".")
  if (is.null(tree$edge.length)) stop("the tree has no branch lengths.")
  if (!is.binary(tree)) stop("The tree must be binary. Please use `ape::multi2di` before proceeding.")
  stopifnot(all.equal(matrix(as.integer(tree$edge), ncol = 2), tree$edge))
  if (method == "random.root" && is.null(root.value)) stop("Starting value must be specified for root node in the `random.root` method.")
  if (method == "random.root" && (is.null(tree$root.edge) || tree$root.edge == 0)) stop("In the random root model, the `root.edge` must be non NULL and non zero.")
  if ((method == "fixed.root") && (is.null(root.value))) stop ("Starting value must be specified for root node in the `fixed.root` method.")
  if (method == "reml" && !is.null(root.value)) stop("In the reml model, `root.value` cannot be specified.")
  tree$edge <- matrix(as.integer(tree$edge), ncol = 2)
  # rootTip
  if (!is.null(rootTip)) {
    rootTip <- rootTip - 1
  } else {
    rootTip <- which.min(colSums(cophenetic.phylo(tree))) - 1 #find_longest_tip_branch(tree) - 1
  }
  # likelihood
  res <-.Call("getLogDensityTipsCauchy", tree, tipTrait, names(tipTrait), root.value, disp, type, rootTip)
  return(res)
}

#' @importFrom foreach %dopar% %do%
NULL

#' @title Posterior density of a node
#'
#' @description
#' Compute the posterior density of a set of node values under a Cauchy process on a phylogenetic tree.
#' 
#' @param node the node for which to compute the posterior density.
#' @param vals the table of values where the density should be computed.
#' @inheritParams logDensityTipsCauchy
#' 
#' @return the posterior density value.
#' 
#' @seealso \code{\link{ancestral}}, \code{\link{fitCauchy}}
#' 
#' @examples
#' phy <- ape::rphylo(5, 0.1, 0)
#' dat <- rTraitCauchy(n = 1, phy = phy, model = "cauchy", parameters = list(root.value = 0, disp = 1))
#' posteriorDensityAncestral(7, 0.1, phy, dat, disp = 1)
#' 
#' @author Paul Bastide \email{paul.bastide@m4x.org} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @export
#' 
#' 
posteriorDensityAncestral <- function(node, vals, tree, tipTrait, root.value = NULL, disp, method = c("reml", "random.root", "fixed.root")) {
  # type
  method <- match.arg(method)
  type <- switch(method,
                 reml = 2,
                 random.root = 0,
                 fixed.root = 1)
  vals <- as.double(vals)
  # checks
  if (!inherits(tree, "phylo")) stop("object \"tree\" is not of class \"phylo\".")
  if (is.null(tree$edge.length)) stop("the tree has no branch lengths.")
  if (!is.binary(tree)) stop("The tree must be binary. Please use `ape::multi2di` before proceeding.")
  if (node <= length(tree$tip.label)) stop("Ancestral reconstruction is only allowed for ancestral nodes.")
  if (node > length(tree$tip.label) + Nnode(tree)) stop ("This node does not exist in the tree.")
  if ((method == "fixed.root") && (node == length(tree$tip.label) + 1)) stop ("Ancestral state reconstruction is not allowed for the root with the fixed root model.")
  if ((method == "fixed.root") && (is.null(root.value))) stop ("Starting value must be specified for root node in the `fixed.root` method.")
  if (method == "random.root" && is.null(root.value)) stop("Starting value must be specified for root node in the `random.root` method.")
  if (method == "random.root" && (is.null(tree$root.edge) || tree$root.edge == 0)) stop("In the random root model, the `root.edge` must be non NULL and non zero.")
  if (method == "reml" && !is.null(root.value)) stop("In the reml model, `root.value` cannot be specified.")
  res <-.Call("getPosteriorLogDensityAncestralCauchy", node - 1, vals, tree, tipTrait, names(tipTrait), root.value, disp, type)
  return(exp(res))
}


##
#' @title Posterior density of a node
#'
#' @description
#' Compute the posterior density of a node value under a fitted Cauchy process on a phylogenetic tree.
#' 
#' @param x an object of class \code{\link{fitCauchy}} or \code{\link{cauphylm}}.
#' @param node the vector of nodes for which to compute the posterior density. 
#' If not specified, the reconstruction is done on all the nodes.
#' @param values the vector of values where the density should be computed. 
#' If not specified, the reconstruction is done for a grid of \code{n_values} values 
#' between \code{1.5 * min(x$y)} and \code{1.5 * max(x$y)}.
#' @param n_values the number of point for the grid of values. 
#' Default to \code{100}. Ignored if \code{values} is provided.
#' @param n_cores number of cores for the parallelization. Default to 1.
# @param progress_bar if \code{TRUE}, show a progress bar for the computations.
#' @param ... other arguments to be passed to the method.
#' 
#' @details
#' This function assumes a Cauchy Process on the tree with fitted parameters 
#' (see \code{\link{fitCauchy}}),
#' and computes the posterior ancestral density of internal nodes, 
#' conditionally on the vector of tip values.
#' 
#' It computes the posterior density on all the points in \code{values},
#' that should be refined enough to get a good idea of the density curve.
#' 
#' @return an object of S3 class \code{ancestralCauchy},
#'  which is a matrix of posterior values, with nodes in rows and values in columns.
#'  
#' @examples
#' set.seed(1289)
#' # Simulate tree and data
#' phy <- ape::rphylo(10, 0.1, 0)
#' dat <- rTraitCauchy(n = 1, phy = phy, model = "cauchy",
#'                     parameters = list(root.value = 10, disp = 0.1))
#' # Fit the data
#' fit <- fitCauchy(phy, dat, model = "cauchy", method = "reml")
#' # Reconstruct the ancestral nodes
#' anc <- ancestral(fit)
#' plot_asr(fit, anc = anc)
#' plot(anc, type = "l")
#' # Refine grid for node 12 and 17
#' anc2 <- ancestral(fit, node = c(12, 17), n_values = 1000)
#' plot(anc2, type = "l")
#' 
#' @author Paul Bastide \email{paul.bastide@m4x.org} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @seealso \code{\link{fitCauchy}}, \code{\link{cauphylm}}, \code{\link{plot.ancestralCauchy}}, \code{\link{plot_asr}}, \code{\link{increment}}
#' 
#' @export
#'
##
ancestral <- function(x, ...) UseMethod("ancestral")

##
#' @describeIn ancestral \code{\link{cauphylm}} object
#' @export
##
ancestral.cauphylm <- function(x, node, values, n_values = 100, n_cores = 1, ...) {
  if (x$model != "cauchy") stop("Ancestral reconstruction is only available for the Cauchy process.")
  if (!(length(x$coefficients) == 1 && sum(x$X) == length(x$X))) stop("Ancestral reconstruction is only available for the Cauchy regression agains the intercept.")
  if (missing(node)) {
    if (x$method == "fixed.root") {
      node <- (Ntip(x$phy) + 2):(Ntip(x$phy) + Nnode(x$phy))
    } else {
      node <- (Ntip(x$phy) + 1):(Ntip(x$phy) + Nnode(x$phy))
    }
  } else {
    if (any(sapply(node, function(nn) !is.wholenumber(nn)))) stop("The 'node' must be whole numbers.")
  }
  if (missing(values)) {
    values <- seq(ifelse(min(x$y) < 0, 1.5 * min(x$y), 0.5 * min(x$y)),
                  ifelse(max(x$y) > 0, 1.5 * max(x$y), 0.5 * max(x$y)),
                  length.out = n_values)
  }
  anc_fun <- function(x, nn, values) posteriorDensityAncestral(node = nn, vals = values,
                                                               tree = x$phy, tipTrait = x$y,
                                                               root.value = x$coefficients, disp = x$disp,
                                                               method = "fixed.root", ...)
  anc <- parallel_construction(anc_fun, x, node, values, n_cores, progress_bar = FALSE)
  class(anc) <- "ancestralCauchy"
  attr(anc, "edge") <- FALSE
  return(anc)
}

##
#' @describeIn ancestral \code{\link{fitCauchy}} object
#' @export
##
ancestral.cauphyfit <- function(x, node, values, n_values = 100, n_cores = 1, ...) {
  if (x$model != "cauchy") stop("Ancestral reconstruction is only available for the Cauchy process.")
  if (missing(node)) {
    if (x$method == "fixed.root") {
      node <- (Ntip(x$phy) + 2):(Ntip(x$phy) + Nnode(x$phy))
    } else {
      node <- (Ntip(x$phy) + 1):(Ntip(x$phy) + Nnode(x$phy))
    }
  } else {
    if (any(sapply(node, function(nn) !is.wholenumber(nn)))) stop("The 'node' must be whole numbers.")
  }
  if (missing(values)) {
    values <- seq(ifelse(min(x$y) < 0, 1.5 * min(x$y), 0.5 * min(x$y)),
                  ifelse(max(x$y) > 0, 1.5 * max(x$y), 0.5 * max(x$y)),
                  length.out = n_values)
  }
  anc_fun <- function(x, nn, values) posteriorDensityAncestral(node = nn, vals = values,
                                                               tree = x$phy, tipTrait = x$trait,
                                                               root.value = x$x0, disp = x$disp, method = x$method, ...)
  anc <- parallel_construction(anc_fun, x, node, values, n_cores, progress_bar = FALSE)
  class(anc) <- "ancestralCauchy"
  attr(anc, "edge") <- FALSE
  return(anc)
}

parallel_construction <- function(anc_fun, x, node, values, n_cores, progress_bar) {
  
  if (length(node) == 1) {
    anc <- anc_fun(x, node, values)
    anc <- matrix(anc, ncol = length(values))
  } else if (n_cores > 1) {
    # combine_fun <- ifelse(progress_bar, rbindProgress(length(node)), rbind)
    combine_fun <- rbind
    
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    
    nn <- NULL
    anc <- foreach::foreach(nn = node, .packages = "cauphy", .combine = combine_fun) %dopar% {
      anc_fun(x, nn, values)
    }
    
  } else {
    combine_fun <- rbind
    
    nn <- NULL
    anc <- foreach::foreach(nn = node, .packages = "cauphy", .combine = combine_fun) %do% {
      anc_fun(x, nn, values)
    }
  }
  
  rownames(anc) <- node
  colnames(anc) <- values
  return(anc)
}

# # https://gist.github.com/kvasilopoulos/d49499ea854541924a8a4cc43a77fed0
# rbindProgress <- function(iterator){
#   pb <- txtProgressBar(min = 0, max = iterator - 1, style = 3)
#   count <- 0
#   return(
#     function(...) {
#       count <<- count + length(list(...)) - 1
#       setTxtProgressBar(pb, count)
#       flush.console()
#       rbind(...) 
#     }
#   )
# }


#' @title Posterior density of an increment
#'
#' @description
#' Compute the posterior density of a set of branch increments under a Cauchy process on a phylogenetic tree.
#' 
#' @param node the node ending the branch for which to compute the posterior density of the increment.
#' @param vals the table of values where the density should be computed.
#' @inheritParams logDensityTipsCauchy
#' 
#' @return the posterior density value.
#' 
#' @examples
#' set.seed(1289)
#' phy <- ape::rphylo(5, 0.1, 0)
#' dat <- rTraitCauchy(n = 1, phy = phy, model = "cauchy", parameters = list(root.value = 0, disp = 1))
#' posteriorDensityIncrement(2, 0.1, phy, dat, disp = 1)
#' 
#' @author Paul Bastide \email{paul.bastide@m4x.org} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @seealso \code{\link{increment}}, \code{\link{fitCauchy}}
#' 
#' @export
#' 
#' 
posteriorDensityIncrement <- function(node, vals, tree, tipTrait, root.value = NULL, disp, method = c("reml", "random.root", "fixed.root")) {
  ntaxa <- length(tree$tip.label)
  # type
  method <- match.arg(method)
  type <- switch(method,
                 reml = 2,
                 random.root = 0,
                 fixed.root = 1)
  vals <- as.double(vals)
  # checks
  if (!inherits(tree, "phylo")) stop("object \"tree\" is not of class \"phylo\".")
  if (is.null(tree$edge.length)) stop("the tree has no branch lengths.")
  if (!is.binary(tree)) stop("The tree must be binary. Please use `ape::multi2di` before proceeding.")
  if (node > length(tree$tip.label) + Nnode(tree)) stop ("This node does not exist in the tree.")
  if ((method == "fixed.root") && (node == length(tree$tip.label) + 1)) stop ("Ancestral increment reconstruction is not allowed for the root branch with the fixed root model.")
  if ((method == "reml") && (node == length(tree$tip.label) + 1)) stop ("Ancestral increment reconstruction is not allowed for the root branch with the reml model.")
  if ((method == "fixed.root") && (is.null(root.value))) stop ("Starting value must be specified for root node in the `fixed.root` method.")
  if (method == "random.root" && is.null(root.value)) stop("Starting value must be specified for root node in the `random.root` method.")
  if (method == "random.root" && (is.null(tree$root.edge) || tree$root.edge == 0)) stop("In the random root model, the `root.edge` must be non NULL and non zero.")
  if (method == "reml" && !is.null(root.value)) stop("In the reml model, `root.value` cannot be specified.")
  
  if ((method == "fixed.root") && 
      (node <= length(tree$tip.label)) && # tip
      (node %in% tree$edge[tree$edge[, 1] == (length(tree$tip.label) + 1), 2])) { # descending from root
    warning(paste0("This branch ends at a tip, and the root is fixed: the posterior increment density is a Dirac in ", tipTrait[tree$tip.label[node]] - root.value, "."))
  }
  
  res <-.Call("getPosteriorLogDensityIncrementCauchy", node - 1, vals, tree, tipTrait, names(tipTrait), root.value, disp, type)
  return(exp(res))
}


##
#' @title Posterior density of an increment
#'
#' @description
#' Compute the posterior density of a branch increment under a fitted Cauchy process on a phylogenetic tree.
#' 
#' @param x an object of class \code{\link{cauphylm}}.
#' @param node vector of nodes ending the branches for which to compute the posterior density of the increment. If not specified, the reconstruction is done on all the possible edges.
#' @param values the vector of values where the density should be computed. If not specified, the reconstruction is done for a grid of \code{n_values} values between \code{-1.5 * maxdiff} and \code{1.5 * maxdiff}, where \code{maxdiff} is the difference between the larger and smaller tip value.
#' @param n_values the number of point for the grid of values. Default to \code{100}. Ignored if \code{values} is provided.
#' @param n_cores number of cores for the parallelization. Default to 1.
# @param progress_bar if \code{TRUE}, show a progress bar for the computations.
#' @param ... other arguments to be passed to the method.
#' 
#' @details
#' This function assumes a Cauchy Process on the tree with fitted parameters 
#' (see \code{\link{fitCauchy}}),
#' and computes the posterior ancestral density of trait increments at branches
#' (ie, the difference between the traits value at the end and beginning of the branch),
#' conditionally on the vector of tip values.
#' 
#' It computes the posterior density on all the points in \code{values},
#' that should be refined enough to get a good idea of the density curve.
#'  
#' @return an object of S3 class \code{ancestralCauchy},
#'  which is a matrix of posterior increment values, with nodes in rows and values in columns.
#' 
#' @examples
#' set.seed(1289)
#' # Simulate tree and data
#' phy <- ape::rphylo(10, 0.1, 0)
#' dat <- rTraitCauchy(n = 1, phy = phy, model = "cauchy",
#'                     parameters = list(root.value = 10, disp = 0.1))
#' # Fit the data
#' fit <- fitCauchy(phy, dat, model = "cauchy", method = "reml")
#' # Reconstruct the ancestral increments
#' inc <- increment(fit)
#' plot_asr(fit, inc = inc)
#' plot(inc, node = c(2, 3, 4, 5, 8, 9), type = "l")
#' # Refine grid for edges ending at tips 3 and 8
#' inc2 <- increment(fit, node = c(3, 8), values = seq(-3, 3, 0.01))
#' plot(inc2, type = "l")
#' 
#' @author Paul Bastide \email{paul.bastide@m4x.org} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @seealso \code{\link{fitCauchy}}, \code{\link{cauphylm}}, \code{\link{plot.ancestralCauchy}}, \code{\link{plot_asr}}, \code{\link{ancestral}}
#' 
#' @export
#'
##
increment <- function(x, ...) UseMethod("increment")

##
#' @describeIn increment \code{\link{cauphylm}} object
#' @export
##
increment.cauphylm <- function(x, node, values, n_values = 100, n_cores = 1, ...) {
  if (x$model != "cauchy") stop("Ancestral reconstruction is only available for the Cauchy process.")
  if (!(length(x$coefficients) == 1 && sum(x$X) == length(x$X))) stop("Ancestral reconstruction is only available for the Cauchy regression agains the intercept.")
  if (missing(node)) {
    node <- 1:(Ntip(x$phy) + Nnode(x$phy))
    root_node <- Ntip(x$phy) + 1
    forbidden_nodes <- NULL
    if(x$method != "random.root") forbidden_nodes <- root_node
    node <- node[!(node %in% forbidden_nodes)]
  } else {
    if (any(sapply(node, function(nn) !is.wholenumber(nn)))) stop("The 'node' must be whole numbers.")
  }
  if (missing(values)) {
    maxdiff <- max(x$y) - min(x$y)
    values <- seq(-1.5 * maxdiff, 1.5 * maxdiff, length.out = n_values)
  }
  inc_fun <- function(x, nn, values) posteriorDensityIncrement(node = nn, vals = values,
                                                               tree = x$phy, tipTrait = x$y,
                                                               root.value = x$coefficients, disp = x$disp,
                                                               method = "fixed.root", ...)
  inc <- parallel_construction(inc_fun, x, node, values, n_cores, progress_bar = FALSE)
  class(inc) <- "ancestralCauchy"
  attr(inc, "edge") <- TRUE
  return(inc)
}

##
#' @describeIn increment \code{\link{fitCauchy}} object
#' @export
##
increment.cauphyfit <- function(x, node, values, n_values = 100, n_cores = 1, ...) {
  if (x$model != "cauchy") stop("Ancestral reconstruction is only available for the Cauchy process.")
  if (missing(node)) {
    node <- 1:(Ntip(x$phy) + Nnode(x$phy))
    root_node <- Ntip(x$phy) + 1
    forbidden_nodes <- NULL
    if(x$method != "random.root") forbidden_nodes <- root_node
    node <- node[!(node %in% forbidden_nodes)]
  } else {
    if (any(sapply(node, function(nn) !is.wholenumber(nn)))) stop("The 'node' must be whole numbers.")
  }
  if (missing(values)) {
    maxdiff <- max(x$y) - min(x$y)
    values <- seq(-1.5 * maxdiff, 1.5 * maxdiff, length.out = n_values)
  }
  inc_fun <- function(x, nn, values) posteriorDensityIncrement(node = nn, vals = values,
                                                               tree = x$phy, tipTrait = x$trait,
                                                               root.value = x$x0, disp = x$disp,
                                                               method = x$method, ...)
  inc <- parallel_construction(inc_fun, x, node, values, n_cores, progress_bar = FALSE)
  class(inc) <- "ancestralCauchy"
  attr(inc, "edge") <- TRUE
  return(inc)
}

##
#'@importFrom graphics close.screen screen split.screen title
NULL

##
#' @title Plot for class \code{ancestralCauchy}
#'
#' @description
#' This function takes an object of class \code{ancestralCauchy}, result of function
#' \code{\link{ancestral}} or \code{\link{increment}}, and plots the reconstructed states for given nodes.
#'
#' @param x an object of class \code{ancestralCauchy}, result of function
#' \code{\link{ancestral}} or \code{\link{increment}}.
#' @param node the vector of nodes where to plot the ancestral reconstruction.
#' @param n.col the number of columns on which to display the plot. Can be left blank.
#' @param ... further arguments to be passed to \code{\link{plot}}.
#' 
#' @return
#' NULL
#' 
#' @examples
#' set.seed(1289)
#' # Simulate tree and data
#' phy <- ape::rphylo(10, 0.1, 0)
#' dat <- rTraitCauchy(n = 1, phy = phy, model = "cauchy",
#'                     parameters = list(root.value = 10, disp = 0.1))
#' # Fit the data
#' fit <- fitCauchy(phy, dat, model = "cauchy", method = "reml")
#' # Reconstruct the ancestral values
#' inc <- increment(fit, node = c(3, 8), values = seq(-3, 3, 0.01))
#' plot(inc, type = "l")
#' anc <- ancestral(fit, node = c(12, 17), n_values = 1000)
#' plot(anc, type = "l")
#' 
#' @seealso  \code{\link{plot_asr}}, \code{\link{ancestral}}, \code{\link{increment}}, \code{\link{fitCauchy}}
#' 
#' @export
#'

plot.ancestralCauchy <- function(x, node, n.col, ...){
  values <- as.numeric(colnames(x))
  all_nodes <- as.numeric(rownames(x))
  if (missing(node)) {
    node <- all_nodes
  } else {
    if (any(sapply(node, function(nn) !is.wholenumber(nn)))) stop("The 'node' must be whole numbers.")
  }
  nn <- match(node, rownames(x))
  if (anyNA(nn)) {
    message(paste0("Nodes ", paste(node[is.na(nn)], collapse = ", "), " are not in the ancestralCauchy reconstruction object. They will not be plotted."))
    nn <- nn[!is.na(nn)]
  }
  if (length(nn) == 0) stop("There are no nodes left to plot.")
  if (missing(n.col)) {
    n.col <- ifelse(length(nn) <= 3, length(nn), 3)
  }
  n.lines <- (length(nn) %/% n.col) + ifelse(length(nn) %% n.col == 0, 0, 1)
  y_lab <- ifelse(attr(x, "edge"), "Posterior density of increment value", "Posterior density of trait value")
  title_plot <- ifelse(attr(x, "edge"), "Edge ending at node ", "Node ")
  scr <- split.screen(c(n.lines, n.col))
  on.exit(close.screen(all.screens = TRUE))
  for (i in seq_len(length(nn))) {
    screen(scr[i])
    dens <- x[nn[i], ]
    plot(values, dens, ylab = y_lab, ...)
    title(paste0(title_plot, rownames(x)[nn[i]]))
  }
}
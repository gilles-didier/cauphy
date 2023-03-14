#' @useDynLib cauphy, .registration=TRUE
#' @import ape
# @importFrom utils flush.console
#' @importFrom methods is
NULL

#library(ade4)
#library(ape)
#dyn.load("treeR.so")

#testTree <- function(tree, part) {
#phy <- read.tree(tree)
#.Call("drawPhylo", phy, "bof.tex")
#list <- scan(part, what = " ")
#new <- .Call("prunePhylo", phy, list)
#phy2 <- read.tree(text = new)
#}

#' @title Print a tree
#'
#' @description
#' \code{printRtree} prints a tree in Newick format
#' 
#' @param tree a phylogeny in \code{\link{ape}} \code{\link[ape]{phylo}} format.
#' 
#' @return No value returned.
#' 
#' @author Paul Bastide \email{paul.bastide@umontpellier.fr} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @keywords internal programming
#' 
printRTreeTest <- function(tree) {
	res <-.Call("printRTree", tree)
	res
}


#' @title Cauchy Trait Simulation
#'
#' @description
#' Simulate a continuous trait using the Cauchy Process
#' 
#' @param n number of independent replicates
#' @param phy a phylogeny in \code{\link{ape}} \code{\link[ape]{phylo}} format.
#' @param model a phylogenetic model. Default is "cauchy", for the Cauchy process. Alternative are "OU", "lambda", "kappa", and "delta".
#' @param parameters list of parameters for the model (see Details).
#' 
#' @return
#' If n=1, a numeric vector with names from the tip labels in the tree.
#' For more than 1 replicate, a matrix with the tip labels as row names, and one column per replicate.
#' 
#' @details
#' The default choice of parameters is as follow:
#' \itemize{
#'   \item{\code{model = cauchy}}{ 
#'   \code{root.value = 0}, \code{disp = 1}
#'   }
#'   \item{\code{model = lambda}}{ 
#'   \code{root.value = 0}, \code{disp = 1}, \code{lambda = 1}
#'   }
#'   \item{\code{model = kappa}}{ 
#'   \code{root.value = 0}, \code{disp = 1}, \code{kappa = 1}
#'   }
#'   \item{\code{model = delta}}{ 
#'   \code{root.value = 0}, \code{disp = 1}, \code{delta = 1}
#'   }
#'   \item{\code{model = OU}}{ 
#'   \code{root.value = 0}, \code{disp = 1}, \code{optimal.value = 0}, \code{selection.strength = 0}
#'   }
#' }
#' 
#' @seealso \code{\link[phylolm]{rTrait}}, \code{\link[ape]{rTraitCont}}
#' 
#' @examples
#' phy <- ape::rphylo(5, 0.1, 0)
#' y = rTraitCauchy(n = 1, phy = phy, model = "lambda", parameters = list(root.value = 0, disp = 0.1))
#' 
#' @author Paul Bastide \email{paul.bastide@umontpellier.fr} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @export
#' 
rTraitCauchy <- function(n = 1, phy,
                         model = c("cauchy", "OU", "lambda", "kappa", "delta"),
                         parameters = NULL) {
  ## Check parameters
  if (is.null(n) || length(n) > 1) stop("n needs to be an integer (number of replicates)")
  n = as.numeric(n)
  if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
  if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
  if (is.null(phy$tip.label)) stop("the tree has no tip labels.")
  phy <- reorder(phy, "pruningwise")
  ## Model and parameters
  model = match.arg(model)
  parameters.default = c(0, 1, 0, 0, 1, 1, 1)
  names(parameters.default) <- c("root.value", "disp", 
                                 "optimal.value", "selection.strength",
                                 "lambda", "kappa", "delta")
  if (is.null(parameters)) {
    parameters <- parameters.default
  } else {
    if (!is.list(parameters)) stop("please specify parameters as a list().")
    match_params <- match(names(parameters), names(parameters.default))
    if (any(is.na(match_params))) stop(paste0("Parameters ",
                                              paste(names(parameters)[is.na(match_params)], collapse = ", "),
                                              " are unknown."))
    specified <- !c(is.null(parameters$root.value), is.null(parameters$disp),
                    is.null(parameters$optimal.value), is.null(parameters$selection.strength),
                    is.null(parameters$lambda), is.null(parameters$kappa), is.null(parameters$delta))
    parameters.user <- c(parameters$root.value, parameters$disp,
                         parameters$optimal.value, parameters$selection.strength,
                         parameters$lambda, parameters$kappa, parameters$delta)
    parameters <- parameters.default
    parameters[specified] <- parameters.user
  }
  p <- as.list(parameters)
  if (model == "OU" & p$selection.strength == 0) model <- "cauchy"
  if (model == "lambda" & p$lambda == 1) model <- "cauchy"
  if (model == "kappa" & p$kappa == 1) model <- "cauchy"
  if (model == "delta" & p$delta == 1) model <- "cauchy"
  # Simulations
  if (model == "OU") {
    sim <- vapply(seq_len(n),
                  FUN = function(i) simulateTipsCauchyOU(tree = phy, start = p$root.value, disp = p$disp, str = p$selection.strength, mean = p$optimal.value),
                  FUN.VALUE = rep(0.0, length(phy$tip.label)))
  } else {
    if (model %in% c("lambda", "kappa", "delta")) {
      phy <- transf.branch.lengths(phy, model, parameters = p, check.pruningwise = F)$tree
    }
    sim <- vapply(seq_len(n),
                  FUN = function(i) simulateTipsCauchy(tree = phy, start = p$root.value, disp = p$disp),
                  FUN.VALUE = rep(0.0, length(phy$tip.label)))
  }
  rownames(sim) <- phy$tip.label
  if (n == 1) {
    sim <- as.vector(sim)
    names(sim) <- phy$tip.label
  }
  return(sim)
}


#' @title Simulate using the Cauchy Process
#'
#' @description
#' Simulate tip values with a Cauchy process
#' 
#' @param tree a phylogeny in \code{\link{ape}} \code{\link[ape]{phylo}} format.
#' @param start the initial root trait value.
#' @param disp the dispersion parameter of the Cauchy process.
#' 
#' @return a vector of simulated values.
#' 
#' @author Paul Bastide \email{paul.bastide@umontpellier.fr} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @keywords internal
#' 
#' 
simulateTipsCauchy <- function(tree, start, disp) {
	res <-.Call("SimulateTipsCauchy", tree, start, disp)
	res
}

#' @title Simulate using the OU Cauchy Process
#'
#' @description
#' Simulate tip values with a OU Cauchy process
#' 
#' @inheritParams simulateTipsCauchy
#' @param str the selection strength of the OU.
#' @param mean the optimal value of the OU.
#' 
#' @return a vector of simulated values.
#' 
#' @author Paul Bastide \email{paul.bastide@umontpellier.fr} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @keywords internal
#' 
#' 
simulateTipsCauchyOU <- function(tree, start, disp, str, mean) {
	res <-.Call("SimulateTipsCauchyOU", tree, start, disp, str, mean)
	res
}

#' @title Log Density of a Cauchy Process
#'
#' @description
#' Compute the log density of a Cauchy process on a phylogenetic tree.
#' 
#' @param tipTrait a names vector of tip trait values, with names matching the tree labels.
#' @param rootTip the tip used to re-root the tree, when the REML method is used. If \code{NULL}, the tip with the longest branch is used. Ignored in `method != "reml"`.
#' @inheritParams simulateTipsCauchy
#' @inheritParams fitCauchy
#' 
#' @return the density value.
#' 
#' @examples
#' phy <- ape::rphylo(5, 0.1, 0)
#' dat <- rTraitCauchy(n = 1, phy = phy, model = "cauchy", parameters = list(root.value = 0, disp = 1))
#' logDensityTipsCauchy(phy, dat, 0, 1, method = "fixed.root")
#' 
#' @author Paul Bastide \email{paul.bastide@umontpellier.fr} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @export
#' 
#' 
logDensityTipsCauchy <- function(tree, tipTrait, start = NULL, disp, method = c("reml", "random.root", "fixed.root"), rootTip = NULL) {
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
  if (method == "random.root" && is.null(start)) stop("Starting value must be specified for root node in the `random.root` method.")
  if (method == "random.root" && (is.null(tree$root.edge) || tree$root.edge == 0)) stop("In the random root model, the `root.edge` must be non NULL and non zero.")
  if ((method == "fixed.root") && (is.null(start))) stop ("Starting value must be specified for root node in the `fixed.root` method.")
  if (method == "reml" && !is.null(start)) stop("In the reml model, `start` cannot be specified.")
  tree$edge <- matrix(as.integer(tree$edge), ncol = 2)
  # rootTip
  if (!is.null(rootTip)) {
    rootTip <- rootTip - 1
  } else {
    rootTip <- which.min(colSums(cophenetic.phylo(tree))) - 1 #find_longest_tip_branch(tree) - 1
  }
  # likelihood
	res <-.Call("getLogDensityTipsCauchy", tree, tipTrait, names(tipTrait), start, disp, type, rootTip)
	return(res)
}


#' @title Find tip with longuest branch
#' 
#' @param tree a phylogenetic tree.
#' 
#' @return tip with the longuest branch.
#' 
#' @author Paul Bastide \email{paul.bastide@umontpellier.fr} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @keywords internal
#' 
#' 
find_longest_tip_branch <- function(tree) {
  n_tips <- length(tree$tip.label)
  tip_branches <- tree$edge[tree$edge[, 2] <= n_tips, ]
  tip_lengths <- tree$edge.length[tree$edge[, 2] <= n_tips]
  rootTip <- tip_branches[which.max(tip_lengths), 2]
  return(rootTip)
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
#' @examples
#' phy <- ape::rphylo(5, 0.1, 0)
#' dat <- rTraitCauchy(n = 1, phy = phy, model = "cauchy", parameters = list(root.value = 0, disp = 1))
#' posteriorDensityAncestral(7, 0.1, phy, dat, disp = 1)
#' 
#' @author Paul Bastide \email{paul.bastide@umontpellier.fr} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @export
#' 
#' 
posteriorDensityAncestral <- function(node, vals, tree, tipTrait, start = NULL, disp, method = c("reml", "random.root", "fixed.root")) {
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
  if ((method == "fixed.root") && (is.null(start))) stop ("Starting value must be specified for root node in the `fixed.root` method.")
  if (method == "random.root" && is.null(start)) stop("Starting value must be specified for root node in the `random.root` method.")
  if (method == "random.root" && (is.null(tree$root.edge) || tree$root.edge == 0)) stop("In the random root model, the `root.edge` must be non NULL and non zero.")
  if (method == "reml" && !is.null(start)) stop("In the reml model, `start` cannot be specified.")
  res <-.Call("getPosteriorLogDensityAncestralCauchy", node - 1, vals, tree, tipTrait, names(tipTrait), start, disp, type)
  return(exp(res))
}


##
#' @title Posterior density of a node
#'
#' @description
#' Compute the posterior density of a node value under a Cauchy process on a phylogenetic tree.
#' 
#' @param x an object of class \code{\link{fitCauchy}} or \code{\link{cauphylm}}.
#' @param node the vector of nodes for which to compute the posterior density. If not specified, the reconstruction is done on all the nodes.
#' @param values the vector of values where the density should be computed. If not specified, the reconstruction is done for a grid of \code{n_values} values between \code{1.5 * min(x$y)} and \code{1.5 * max(x$y)}.
#' @param n_values the number of point for the grid of values. Default to \code{100}. Ignored if \code{values} is provided.
#' @param n_cores number of cores for the parallelization. Default to 1.
# @param progress_bar if \code{TRUE}, show a progress bar for the computations.
#' @param ... other arguments to be passed to the method.
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
#' @author Paul Bastide \email{paul.bastide@umontpellier.fr} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @seealso \code{\link{fitCauchy}}, \code{\link{cauphylm}}, \code{\link{plot.ancestralCauchy}}, \code{\link{increment}}
#' @export
#'
##
ancestral <- function(x, ...) UseMethod("ancestral")

##
#' @describeIn ancestral \code{\link{cauphylm}} object
#' @export
##
ancestral.cauphylm <- function(x, node, values, n_values = 100, n_cores = 1, ...) {
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
                                                               start = x$coefficients, disp = x$disp,
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
                                                               start = x$x0, disp = x$disp, method = x$method, ...)
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
#' @author Paul Bastide \email{paul.bastide@umontpellier.fr} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @export
#' 
#' 
posteriorDensityIncrement <- function(node, vals, tree, tipTrait, start = NULL, disp, method = c("reml", "random.root", "fixed.root")) {
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
  if ((method == "fixed.root") && (is.null(start))) stop ("Starting value must be specified for root node in the `fixed.root` method.")
  if (method == "random.root" && is.null(start)) stop("Starting value must be specified for root node in the `random.root` method.")
  if (method == "random.root" && (is.null(tree$root.edge) || tree$root.edge == 0)) stop("In the random root model, the `root.edge` must be non NULL and non zero.")
  if (method == "reml" && !is.null(start)) stop("In the reml model, `start` cannot be specified.")
  
  if ((method == "fixed.root") && 
      (node <= length(tree$tip.label)) && # tip
      (node %in% tree$edge[tree$edge[, 1] == (length(tree$tip.label) + 1), 2])) { # descending from root
    warning(paste0("This branch ends at a tip, and the root is fixed: the posterior increment density is a Dirac in ", tipTrait[tree$tip.label[node]] - start, "."))
  }
  
  res <-.Call("getPosteriorLogDensityIncrementCauchy", node - 1, vals, tree, tipTrait, names(tipTrait), start, disp, type)
  return(exp(res))
}


##
#' @title Posterior density of an increment
#'
#' @description
#' Compute the posterior density of a branch increment under a Cauchy process on a phylogenetic tree.
#' 
#' @param x an object of class \code{\link{cauphylm}}.
#' @param node vector of nodes ending the branches for which to compute the posterior density of the increment. If not specified, the reconstruction is done on all the possible edges.
#' @param values the vector of values where the density should be computed. If not specified, the reconstruction is done for a grid of \code{n_values} values between \code{-1.5 * maxdiff} and \code{1.5 * maxdiff}, where \code{maxdiff} is the difference between the larger and smaller tip value.
#' @param n_values the number of point for the grid of values. Default to \code{100}. Ignored if \code{values} is provided.
#' @param n_cores number of cores for the parallelization. Default to 1.
# @param progress_bar if \code{TRUE}, show a progress bar for the computations.
#' @param ... other arguments to be passed to the method.
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
#' # Reconstruct the ancestral nodes
#' inc <- increment(fit)
#' plot_asr(fit, inc = inc)
#' plot(inc, node = c(2, 3, 4, 5, 8, 9), type = "l")
#' # Refine grid for edges ending at tips 3 and 8
#' inc2 <- increment(fit, node = c(3, 8), values = seq(-3, 3, 0.01))
#' plot(inc2, type = "l")
#' 
#' @author Paul Bastide \email{paul.bastide@umontpellier.fr} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @seealso \code{\link{fitCauchy}}, \code{\link{cauphylm}}, \code{\link{plot.ancestralCauchy}}, \code{\link{increment}}
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
                                                               start = x$coefficients, disp = x$disp,
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
                                                               start = x$x0, disp = x$disp,
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
#' \code{\link{ancestral}}, and plots the reconstructed states for given nodes.
#'
#' @param x an object of class \code{ancestralCauchy}, result of function
#' \code{\link{ancestral}}.
#' @param node the vector of nodes where to plot the ancestral reconstruction.
#' @param n.col the number of columns on which to display the plot. Can be left blank.
#' @param ... further arguments to be passed to \code{\link{plot}}.
#' 
#' @return
#' NULL
#' 
#' 
#' @seealso \code{\link{ancestral}}, \code{\link{fitCauchy}}
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
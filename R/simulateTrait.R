#' @importFrom stats reorder rnorm rpois runif

#' @title Simulate a Levy process
#'
#' @description
#' This function uses \code{\link[ape]{rTraitCont}} to simulate a Levy process
#' with Brownian variance \code{brownianVar}, 
#' jumps occurring as a Poisson process with parameter \code{lambda},
#' each jumps Gaussian with variance \code{jumpVar}.
#' 
#' @references 
#' Duchen P., Leuenberger C., Szilágyi S. M., Harmon L. J., Eastman J. M., Schweizer M., Wegmann D. 2017. Inference of Evolutionary Jumps in Large Phylogenies using Lévy Processes. Systematic Biology. 66:950–963.
#' 
#' @seealso
#' \code{\link[ape]{rTraitCont}}
#'
#' @inheritParams ape::rTraitCont
#' 
#' @param brownianVar the variance of the Brownian part of the process.
#' @param jumpVar the variance of the Gaussian jumps.
#' @param lambda  the parameter of the jump Poisson process.
#'
#' @return A named vector of simulated values.
#' 
#' @export
#'
simulateLevy <- function(phy, brownianVar, jumpVar, lambda, ancestor, root.value) {
  simLevy <- function(x, l, brownianVar, jumpVar, lambda) {
    ## BM
    x <- x + rnorm(1, 0, sd = sqrt(l * brownianVar))
    ## Number of jumps
    nj <- rpois(1, lambda * l)
    ## add gaussian jumps
    x <- x + sum(rnorm(nj, 0, sd = sqrt(jumpVar)))
    return(x)
  }
  return(ape::rTraitCont(phy, model = simLevy,
                         ancestor = ancestor, root.value = root.value,
                         brownianVar = brownianVar, jumpVar = jumpVar, lambda = lambda))
}

#' @title Simulate a Levy process
#'
#' @description
#' This function uses \code{\link[ape]{rTraitCont}} to simulate a Levy process
#' with Brownian variance \code{brownianVar}, 
#' conditionally on \code{nJumps} occurring on the phylogeny,
#' each jumps Gaussian with variance \code{jumpVar}.
#' 
#' @references 
#' Duchen P., Leuenberger C., Szilágyi S. M., Harmon L. J., Eastman J. M., Schweizer M., Wegmann D. 2017. Inference of Evolutionary Jumps in Large Phylogenies using Lévy Processes. Systematic Biology. 66:950–963.
#' 
#' @seealso
#' \code{\link[ape]{rTraitCont}}
#'
#' @inheritParams ape::rTraitCont
#' 
#' @param brownianVar the variance of the Brownian part of the process.
#' @param jumpVar the variance of the Gaussian jumps.
#' @param nJumps  total number of jumps.
#'
#' @return A named vector of simulated values.
#' 
#' @export
#'
simulateLevyConditional <- function(phy, brownianVar, jumpVar, nJumps, ancestor, root.value) {
  ## Distribution of jumps
  seg_lengths <- cumsum(phy$edge.length)
  jumpPos <- sort(runif(nJumps) * max(seg_lengths))
  jumpEdges <- sapply(jumpPos, function(x) sum(x > seg_lengths) + 1)
  
  ## Simulation function
  simLevy <- function(x, l, edge, brownianVar, jumpVar) {
    ## BM
    x <- x + rnorm(1, 0, sd = sqrt(l * brownianVar))
    ## Number of jumps
    nj <- sum(edge == jumpEdges)
    ## add gaussian jumps
    x <- x + sum(rnorm(nj, 0, sd = sqrt(jumpVar)))
    return(x)
  }
  
  ## Sim
  if (is.null(phy$edge.length)) stop("tree has no branch length")
  if (any(phy$edge.length < 0)) stop("at least one branch length negative")
  phy <- reorder(phy, "postorder")
  n <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- n + 1L
  x <- numeric(n + phy$Nnode)
  x[ROOT] <- root.value
  anc <- phy$edge[, 1]
  des <- phy$edge[, 2]
  el <- phy$edge.length
  for (i in N:1) x[des[i]] <- simLevy(x[anc[i]], el[i], i, brownianVar = brownianVar, jumpVar = jumpVar)
  
  ## Return
  if (ancestor) {
    if (is.null(phy$node.label)) phy <- makeNodeLabel(phy)
    names(x) <- c(phy$tip.label, phy$node.label)
  }
  else {
    x <- x[1:n]
    names(x) <- phy$tip.label
  }
  return(x)
}
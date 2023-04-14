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
#' @author Paul Bastide \email{paul.bastide@m4x.org} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @keywords internal programming
#' 
printRTreeTest <- function(tree) {
  res <-.Call("printRTree", tree)
  res
}

#' @title Invisible Plotting
#'
#' @description Invisible Plotting
#'
#' @param ... arguments to be passed to the plot function.
#' 
#' @details See https://stackoverflow.com/a/20363489
#'
#' @return the result of the plot, without displaying it.
#' 
#' @keywords internal
#'
plot.invisible <- function(...){
  ff <- tempfile()
  png(filename=ff)
  res <- plot(...)
  dev.off()
  unlink(ff)
  res
}

#' @title Re root tree at a tip
#'
#' @description Re root tree at a tip, taking care of the root length.
#' 
#' @param ... dots arguments to be captured
#'
#' @return a named list of the arguments in ...
#'
#' @keywords internal
#'
reroottip <- function(tree, tip) {
  ## Re root with ape
  tree2 <- root(tree, tip)
  
  ## Drop root tip
  tree3 <- drop.tip(tree2, tip)
  
  ## Case tip is an outlier
  if (isTRUE(all.equal(tree, tree2))) {
    anc <- getMRCA(tree, tree$tip.label[-tip])
    root_edge <- which(tree$edge[, 2] == anc)
    tree3$root.edge <- tree2$edge.length[root_edge]
  }
  
  ## Deal with root length
  if (is.null(tree3$root.edge)) tree3$root.edge <- 0.0
  tree3$root.edge <- tree2$edge.length[tree2$edge[, 2] == tip] + tree3$root.edge
  
  return(tree3)
}

#' @title Safely get element of a named vector
#'
#' @description
#' If the element does not exist, return NULL.
#'
#' @param x a named vector
#' @param name the name of the element to retrieve.
#' 
#' @keywords internal
#'
safe_get <- function(x, name) {
  if (!any(grepl(name, names(x)))) {
    return(NULL)
  } else {
    return(x[[name]])
  }
}

# See help(is.integer)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


slog1p <- Vectorize(function(x) {
  # if (x > 0) return(log(1 + x))
  # if (x < 0) return(-log(1 - x))
  # return(0)
  return(x)
})

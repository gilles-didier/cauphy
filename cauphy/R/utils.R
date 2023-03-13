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

#' @title Capture dot arguments
#'
#' @description Capture dot arguments
#' 
#' @details See: http://adv-r.had.co.nz/Computing-on-the-language.html#capturing-dots
#'
#' @param ... dots arguments to be captured
#'
#' @return a named list of the arguments in ...
#'
#' @keywords internal
#'
dots <- function(...) {
  eval(substitute(alist(...)))
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

slog <- Vectorize(function(x) {
  if (x >= 0) return(-log(x))
  if (x < 0) return(log(-x))
})

sexp1p <- Vectorize(function(x) {
  # if (x > 0) return(exp(x) - 1)
  # if (x < 0) return(-exp(-x) + 1)
  # return(0)
  return(x)
})

#' @importFrom graphics legend strwidth
NULL

#' @title Plot Ancestral States Reconstructions
#'
#' @description
#' Plot the ancestral states reconstructions from a fitted Cauchy model.
#'
#' @param x a \code{\link{cauphylm}} or \code{\link{fitCauchy}} object.
#' @param anc (optional) an object of class \code{ancestralCauchy}, obtained with \code{\link{ancestral}}.
#' @param inc (optional) an object of class \code{ancestralCauchy}, obtained with \code{\link{increment}}.
#' @param common_colorscale If both plotted, should the ancestral states and the increment be represented by the same color scale ? Default to \code{FALSE}.
#' @param x.legend,y.legend the x and y co-ordinates to be used to position the legend. They can be specified by keyword or in any way which is accepted by \code{\link[graphics]{legend}}.
# @param n_values if \code{anc} is not provided, the number of points in the grid discretizing the trait values. Default to 1000. Ignored if \code{anc} in not \code{NULL}.
# @param values if \code{anc} is not provided, an optional vector of grid values where the density will be estimated.
# If \code{NULL}, a grid with \code{n_values} going from \code{1.5 * min(x$y)} to \code{1.5 * max(x$y)}
# is chosen. If provided by the user, this overrides the \code{n_values} parameter.
# Ignored if \code{anc} in not \code{NULL}.
# @param mode_proba_method the method to compute the approximate relative importance of each mode.
# One of \code{"value"} (the default) or \code{"integral"}.
# Ignored it parameter \code{thermo} is provided. See details.
# @param plot_legend if \code{TRUE} (the default) a color legend for the trait values is plotted.
#' @inheritParams ape::nodelabels
#' @inheritParams ape::phydataplot
#' @inheritParams ape::plot.phylo
#' @param width.node,height.node,width.edge,height.edge parameters controlling the aspect of thermometers for the nodes and the edges; by default, their width and height are determined automatically.
#' @param ... other parameters to be passed on to \code{\link[ape]{plot.phylo}} or \code{\link[ape]{phydataplot}}.
#' @param x.intersp character interspacing factor for horizontal (x) spacing between symbol and legend text (see \code{\link[graphics]{legend}}).
#'
#' @details 
#' The main plot is done with \code{\link[ape]{plot.phylo}},
#' the node annotation use \code{\link[ape]{nodelabels}}, and the 
#' tip data plot use \code{\link[ape]{phydataplot}}.
#' Please refer to these functions for the details of the parameters.
#' 
#' Note that the thermo plot do not dive the true probabilities of each mode values.
#' Instead, they are simply proportional to the height of each mode, re-scaled so that the scores sum to one.
#' For an exact representation of a node posterior density, please plot it separately, using function \code{\link{ancestral}}.
#' 
# Parameter \code{mode_proba_method} controls what is plotted on each node.
# If set to \code{mode_proba_method="value"} (the default), then the importance of each node is taken as proportional
# to the height of the mode, re-scaled so that the scores sum to one.
# Not that this is not a real probability measure.
# If set to \code{mode_proba_method="integral"}, then the probability of each node is approximated by the integral of
# the density around the mode. The probabilities are still not exact and need to be re-scaled to sum to one,
# but might be more realistic. This method can be much slower on large datasets.
#' 
#' @return Plot of the tree with ancestral state reconstruction.
#' 
#' @seealso \code{\link{cauphylm}}, \code{\link{fitCauchy}}, \code{\link{ancestral}}, \code{\link{increment}},
#' \code{\link[ape]{plot.phylo}}, \code{\link[ape]{phydataplot}}, \code{\link[ape]{nodelabels}}
#' 
#' @export
#'
plot_asr <- function(x, 
                     anc = NULL, inc = NULL,
                     common_colorscale = FALSE,
                     x.legend = "topleft",
                     y.legend = NULL,
                     ## nodelabels
                     adj = c(0.5, 0.5),
                     piecol = NULL,
                     # horiz = FALSE,
                     width.node = NULL, height.node = NULL,
                     width.edge = NULL, height.edge = NULL,
                     ## phydataplot
                     style = "bars", offset = 1, scaling = 1,
                     x.lim = NULL, x.intersp = NULL, ...) UseMethod("plot_asr")

##
#' @export
##
plot_asr <- function(x, 
                     anc = NULL, inc = NULL,
                     common_colorscale = FALSE,
                     x.legend = "topleft",
                     y.legend = NULL,
                     ## nodelabels
                     adj = c(0.5, 0.5),
                     piecol = NULL,
                     width.node = NULL, height.node = NULL,
                     width.edge = NULL, height.edge = NULL,
                     ## phydataplot
                     style = "bars", offset = 1, scaling = 1,
                     x.lim = NULL, x.intersp = NULL, ...) {

  traitName <- "trait"
  if (length(x$formula) > 0) traitName <- all.vars(x$formula)[1]
  
  ## xlim
  lastPP <- plot.invisible(x$phy, ...)
  xlimtree <- lastPP$x.lim
  ## Space for data
  px <- pretty(c(0, x$y))
  range_data <- scaling * abs(max(px) - min(px)) + offset
  if (is.null(x.lim)) x.lim <- c(xlimtree[1], xlimtree[2] + range_data)
  ## Plot Tree
  plot(x$phy, x.lim = x.lim, ...)
  
  ## Method to compute proba of each node
  mode_proba_method <- "integral" # match.arg(mode_proba_method)
  trans_vals <- function(x) return(x)
  trans_inv_vals <- function(x) return(x)
  
  # ancestral nodes
  values_anc <- NULL
  if (!is.null(anc)) {
    if (!is(anc, "ancestralCauchy")) stop("'anc' must by of S3 class 'ancestralCauchy'.")
    node_anc <- as.numeric(rownames(anc))
    thermo_anc <- find_modes(anc, node_anc, mode_proba_method)
    values_anc <- as.numeric(colnames(anc))
  } else {
    values_anc <- seq(min(x$y), max(x$y), length.out = 1000)
  }
  
  # ancestral increments
  values_inc <- NULL
  increment.proba <- FALSE
  if (!is.null(inc)) {
    if (!is(inc, "ancestralCauchy")) stop("'inc' must by of S3 class 'ancestralCauchy'.")
    node_inc <- as.numeric(rownames(inc))
    thermo_inc <- find_modes(inc, node_inc, mode_proba_method)
    values_inc <- as.numeric(colnames(inc))
    edge_inc <- sapply(node_inc, function(nn) which(x$phy$edge[, 2] == nn))
    if (increment.proba) {
      trans_vals <- pracma::logit
      trans_inv_vals <- pracma::sigmoid
      thermo_inc <- get_p_value_thermo(x, thermo_inc, edge_inc)
      values_inc <-  as.numeric(colnames(thermo_inc))
      # codes_inc <- get_p_value_code(x, thermo_inc, edge_inc)
    }
  } else {
    common_colorscale <- TRUE
  }
  
  ## Color scales
  get_piecol <- function(vv, proba = FALSE) {
    if (!proba) {
      max_abs <- max(abs(vv))
      if ((min(vv) < 0) && (max(vv) > 0)) {
        scale_diverging <- TRUE
        vv <- seq(0, max_abs, length.out = 501)
        vv <- c(seq(-max_abs, -vv[2], length.out = 499), vv)
      } else {
        scale_diverging <- FALSE
        vv <- seq(floor(min(vv)), ceiling(max(vv)), length.out = 1000)
      }
    } else {
      scale_diverging <- TRUE
      vv <- seq(trans_vals(0.001), trans_vals(0.999), length.out = 1000)
    }
    
    if (is.null(piecol)) {
      if (scale_diverging) {
        piecol <- hcl.colors(length(vv), rev = TRUE, palette = "Tropic")
      } else {
        piecol <- hcl.colors(length(vv), rev = TRUE, palette = "inferno")
      }
    }
    return(list(values = vv,
                piecol = piecol))
  }
  
  if (common_colorscale) {
    all.values <- slog1p(c(values_inc, values_anc))
    res <- get_piecol(all.values)
    piecol.anc <- piecol.inc <- res$piecol
    all.values.anc <- all.values.inc <- res$values
  } else {
    res <- get_piecol(values_anc)
    piecol.anc <- res$piecol
    all.values.anc <- res$values
    res <- get_piecol(values_inc, proba = increment.proba)
    piecol.inc <- res$piecol
    all.values.inc <- res$values
  }
  
  ## Plot ancestral nodes
  if (!is.null(anc)) {
    nodelabels(node = node_anc, adj = adj,
               thermo = thermo_anc,
               piecol = piecol.anc[findInterval(values_anc, all.values.anc)],
               horiz = FALSE,
               width = width.node, height = height.node)
  }
  if (!is.null(inc)) {
    edgelabels(edge = edge_inc, adj = adj,
               thermo = thermo_inc,
               piecol = piecol.inc[findInterval(trans_vals(values_inc), all.values.inc)],
               horiz = TRUE,
               width = width.edge, height = height.edge)
    # if (increment.proba) {
    #   edgelabels(edge = edge_inc, adj = adj, pos = 3,
    #              text = codes_inc, frame = "none")
    # }
  }
  
  ## Plot data
  # get order colors
  col_data <- piecol.anc[findInterval(slog1p(x$y), all.values.anc)]
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  y1 <- lastPP$yy[seq_len(length(x$phy$tip.label))]
  o <- order(y1)
  col_data <- col_data[o]
  # offset
  offset_data <- scaling * abs(min(min(px), 0))
  # plot
  phydataplot(x$y, x$phy, col = col_data,
              style = style, offset = offset + offset_data, scaling = scaling,
              # continuous = continuous,
              # width = width,
              # legend = legend,
              # funcol = funcol,
              ...)
  
  ## legend
  get_legend <- function(vv, trans_inv = function(x) return(x)) {
    ticks <- unname(floor(quantile(seq_along(vv), prob = c(0, 0.25, 0.5, 0.75, 1))))
    lgd_ = rep(NA, length(vv))
    lgd_[ticks] <- formatC(trans_inv(vv)[ticks], digits = 1, format = "f") #signif(sexp1p(vv)[ticks], digits = 2)
    return(lgd_)
  }
  lgd_anc <- get_legend(all.values.anc)
  if (common_colorscale) {
    if (is.null(x.intersp)) {
      # lgd_width <- max(strwidth(lgd_anc))
      x.intersp <-  max(nchar(lgd_anc), na.rm = TRUE) - 2
    }
    legend(x.legend, y.legend,
           legend = lgd_anc,
           fill = piecol.anc,
           border = NA,
           y.intersp = 5 / length(all.values.anc),
           cex = 1, adj = 1, x.intersp = x.intersp, text.width = 0)
  } else {
    lgd_inc <- get_legend(all.values.inc, trans_inv_vals)
    if (is.null(x.intersp)) {
      x.intersp <-  max(nchar(c(lgd_anc, lgd_inc)), na.rm = TRUE) - 2 # max(strwidth(c(lgd_anc, lgd_inc), "user"))
    }
    legend(x.legend, y.legend,
           legend = c(lgd_anc, lgd_inc),
           fill = c(piecol.anc, piecol.inc),
           border = NA,
           y.intersp = 5 / (length(all.values.anc)),
           cex = 1, adj = 1, x.intersp = x.intersp, text.width = 0,
           ncol = 2)
  }

}


# #' @title Find the p values of increments
# #' 
# #' @description 
# #' This function fits a standard BM on the tree, and then use the infered variance
# #' to detect "abnormal" increments.
# #' For each branch i, the "pvalue" of a mode is given by the probability of observing
# #' the absolute value of this mode or higher assuming a Gaussian distribution
# #' with mean 0 and variance sigma_BM * t_i.
# #' The pvalue is signed to the sign of the mode.
# #'
# #' @param x a cauchy fitted object
# #' @param thermo_inc a thermo object
# #' @param edge_inc matching edges
# #'
# #' @return a new thermo object, with signed pvalues as values.
# #' 
# #' @keywords internal
# #'
# get_p_value_thermo <- function(x, inc, edge_inc) {
#   browser()
#   thermo_inc_pval <- matrix(0, ncol = 1001, nrow = nrow(inc))
#   pval_grid <- seq(0, 1, length.out = 1001)
#   colnames(thermo_inc_pval) <- pval_grid
#   for (i in seq_len(nrow(inc))) {
#     tt <- inc[i, ]
#     modes <- as.numeric(names(tt[tt > 0]))
#     scale_branch  <- x$disp * x$phy$edge.length[edge_inc[i]]
#     pvals <- stats::pcauchy(tt, location = 0, scale = scale_branch, lower.tail = TRUE)
#     pval_ind <- findInterval(pvals, pval_grid)
#     thermo_inc_pval[i, pval_ind] <- tt[tt > 0]
#   }
#   return(thermo_inc_pval)
#  }

#' @title Find the p values of increments
#' 
#' @description 
#' This function fits a standard BM on the tree, and then use the infered variance
#' to detect "abnormal" increments.
#' For each branch i, the "pvalue" of a mode is given by the probability of observing
#' the absolute value of this mode or higher assuming a Gaussian distribution
#' with mean 0 and variance sigma_BM * t_i.
#' The pvalue is signed to the sign of the mode.
#'
#' @param x a cauchy fitted object
#' @param thermo_inc a thermo object
#' @param edge_inc matching edges
#'
#' @return a new thermo object, with signed pvalues as values.
#' 
#' @keywords internal
#'
get_p_value_thermo <- function(x, thermo_inc, edge_inc) {
  thermo_inc_pval <- matrix(0, ncol = 1000, nrow = nrow(thermo_inc))
  pval_grid <- seq(0, 1, length.out = 1002)
  pval_grid <- pval_grid[-c(1, length(pval_grid))]
  colnames(thermo_inc_pval) <- pval_grid
  for (i in seq_len(nrow(thermo_inc))) {
    tt <- thermo_inc[i, ]
    modes <- as.numeric(names(tt[tt > 0]))
    disp_cau <- x$disp * x$phy$edge.length[edge_inc[i]]
    pvals <- stats::pcauchy(modes, location = 0, scale = disp_cau, lower.tail = TRUE)
    pval_ind <- findInterval(pvals, pval_grid)
    thermo_inc_pval[i, pval_ind] <- tt[tt > 0]
  }
  return(thermo_inc_pval)
}

# get_p_value_thermo <- function(x, thermo_inc, edge_inc) {
#   sigma2 <- phylolm::phylolm(x$y ~ 1, phy = x$phy)$sigma2
#   thermo_inc_pval <- matrix(0, ncol = 1000, nrow = nrow(thermo_inc))
#   pval_grid <- seq(-1, 1, length.out = 1001)
#   pval_grid <- pval_grid[-501]
#   colnames(thermo_inc_pval) <- pval_grid
#   for (i in seq_len(nrow(thermo_inc))) {
#     tt <- thermo_inc[i, ]
#     modes <- as.numeric(names(tt[tt > 0]))
#     var_bm <- sigma2 * x$phy$edge.length[edge_inc[i]]
#     pvals <- stats::pnorm(abs(modes), mean = 0, sd = sqrt(var_bm), lower.tail = FALSE) * sign(modes)
#     pval_ind <- findInterval(pvals, pval_grid)
#     thermo_inc_pval[i, pval_ind] <- tt[tt > 0]
#   }
# return(thermo_inc_pval)
# }

get_p_value_code <- function(x, thermo_inc, edge_inc) {
  codesignif <- c("***", "**", "*", ".", " ")
  codesignifvals <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
  sigma2 <- phylolm::phylolm(x$y ~ 1, phy = x$phy)$sigma2
  all_codes <- rep(NA, length(edge_inc))
  names(all_codes) <- edge_inc
  for (i in seq_len(nrow(thermo_inc))) {
    tt <- thermo_inc[i, ]
    modes <- as.numeric(names(tt[tt > 0]))
    var_bm <- sigma2 * x$phy$edge.length[edge_inc[i]]
    pvals <- stats::pnorm(abs(modes), mean = 0, sd = sqrt(var_bm), lower.tail = FALSE)
    all_codes[i] <- codesignif[findInterval(min(pvals), codesignifvals, left.open = TRUE)]
  }
  return(all_codes)
}

#' @title Find modes of a distribution
#'
#' @param anc an object of class \code{ancestralCauchy}
#' @param node a vector of node number in the tree
#'
#' @return a list, with \code{opt_vals} the abscissa values where the distribution is maximum,
#' and \code{prop_vals} the approximate density under each optimal value, normalized to sum to one.
#' 
#' @keywords internal
#'
find_modes <- function(anc, node, mode_proba_method = "integral") {
  mm <- match(node, rownames(anc))
  if (anyNA(mm)) message(paste0("Nodes ", paste(node[is.na(mm)], collapse = ", "), " have not been reconstructed, they will not be plotted."))
  values <- as.numeric(colnames(anc))
  fun <- function(nn) {
    peaks <- findmaxsgrid(anc[nn == rownames(anc), ], values)
    proba_values <- rep(0, ncol(anc))
    if (mode_proba_method == "value") {
      proba_values[peaks$opt_vals] <- peaks$dens_vals / sum(peaks$dens_vals)
    } else {
      proba_values[peaks$opt_vals] <- peaks$proba_vals / sum(peaks$proba_vals)
    }
    return(proba_values)
  }
  
  thermo <- t(vapply(node, FUN = fun, FUN.VALUE = rep(0.0, ncol(anc))))
  colnames(thermo) <- colnames(anc)
  rownames(thermo) <- node
  return(thermo)
}

findmaxsgrid <- function(dens, values) {
  pp <- pracma::findpeaks(dens)
  proba_vals <- apply(pp, 1, function(ppp) pracma::trapz(values[ppp[3]:ppp[4]], dens[ppp[3]:ppp[4]]))
  return(list(opt_vals = pp[, 2],
              dens_vals = pp[, 1],
              proba_vals = proba_vals))
}

#' @title Find closest value in a vector
#'
#' @param a a value to be matched
#' @param v a vector of grid values
#'
#' @return the index of the closest values of a in v.
#' 
#' @keywords internal
#'
find_closest <- function(a, v){
  which.min(abs(v-a))
}

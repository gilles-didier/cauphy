#' @title Posterior density of a node
#'
#' @description
#' Compute the posterior density of a set of node values under a Cauchy process on a phylogenetic tree.
#' 
#' @param tree the phylogenetic tree.
#' @param fit_lat a \code{\link{fitCauchy}} object for the latitude trait.
#' @param fit_long a \code{\link{fitCauchy}} object for the longitude trait.
#' @param anc_lat an object of class \code{ancestralCauchy}, obtained with \code{\link{ancestral}} on \code{fit_lat}.
#' @param anc_long an object of class \code{ancestralCauchy}, obtained with \code{\link{ancestral}} on \code{fit_long}.
#' @param level the level for the HPD intervals
#' 
#' @return An object of class \code{\link[tidytree]{treedata-class}}, that can
#' be exported with function \code{\link[treeio]{write.beast}}.
#' 
#' @seealso \code{\link{ancestral}}, \code{\link{fitCauchy}}
#' 
#' @author Paul Bastide \email{paul.bastide@m4x.org} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @keywords internal
#' 
#'
ancestral_to_treedata <- function(tree, 
                                  fit_lat, fit_long, 
                                  anc_lat, anc_long,
                                  level = 0.95) {
  
  if (!requireNamespace("matrixStats", quietly = TRUE) || 
      !requireNamespace("treeio", quietly = TRUE) || 
      !requireNamespace("tidytree", quietly = TRUE)) {
    stop(
      "Packages \"treeio\", \"tidytree\" and \"HDInterval\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  ## TODO: check entries
  
  n_nodes <- tree$Nnode
  n_tips <- length(tree$tip.label)
  
  ## Construct HPD intervals
  dens_lat <- anc_to_density(anc_lat)
  dens_long <- anc_to_density(anc_long)
  
  dens_lat_long <- vector(mode = "list", length = length(dens_long))
  for (i in 1:length(dens_long)) {
    dens_lat_long[[i]] <- list(x = dens_lat[[i]]$x,
                               y = dens_long[[i]]$x,
                               z = dens_lat[[i]]$y %*% t(dens_long[[i]]$y))
  }
  names(dens_lat_long) <- names(dens_long)
  
  contours_hpd <- lapply(dens_lat_long, function(dd) get_contour(dd, level))
  
  ## Create data list object
  rec_table <- list(node = seq_len(n_tips + n_nodes))
  
  max_comp <- max(sapply(contours_hpd, length))
  trait_name_lat <- paste0("location1_", level * 100, "%HPD")
  trait_name_long <- paste0("location2_", level * 100, "%HPD")
  
  for (i in seq_len(max_comp)) {
    rec_table[[paste0(trait_name_lat, "_", i)]] <- c(rep(NA, n_tips), lapply(contours_hpd, function(x) get_entry_safe(x, i, "x")))
    rec_table[[paste0(trait_name_long, "_", i)]] <- c(rep(NA, n_tips), lapply(contours_hpd,function(x) get_entry_safe(x, i, "y")))
  }
  
  rec_table[["location1"]] <- c(fit_lat$y, sapply(dens_lat, function(dd) stats::weighted.mean(dd$x, dd$y)))
  rec_table[["location2"]] <- c(fit_long$y, sapply(dens_long, function(dd) stats::weighted.mean(dd$x, dd$y)))
  
  rec_table[["location1_median"]] <- c(rep(NA, n_tips), sapply(dens_lat, function(dd) matrixStats::weightedMedian(dd$x, dd$y)))
  rec_table[["location2_median"]] <- c(rep(NA, n_tips), sapply(dens_long, function(dd) matrixStats::weightedMedian(dd$x, dd$y)))
  
  rec_table <- treeio::as_tibble(rec_table)
  
  ## Create treedata object
  tree_tibble <- treeio::as_tibble(tree)
  # tree_tibble$label <- rec_table$label
  
  tree_data <- treeio::full_join(tree_tibble, rec_table, by = 'node')
  
  return(treeio::as.treedata(tree_data))
  
}

#' @title Get entry from a list
#'
#' @description
#' Get the entry if the index is acceptable, otherwise return NA.
#' 
#' @param x a list
#' @param i the item number
#' @param nn the name to retrieve
#' 
#' @return entry nn from item i of list x, or NA
#' 
#' @keywords internal
#' 
get_entry_safe <- function(x, i, nn) {
  if (i > length(x)) return(NA)
  return(x[[i]][[nn]])
}

#' @title Get Contour
#'
#' @description
#' See \code{\link[emdbook]{HPDregionplot}}
#' 
#' @param den a two dimensional density object
#' @param prob the probaility level
#' 
#' @return contour plot coordinates
#' 
#' @keywords internal
#' 
get_contour <- function(den, prob) {
  dx <- diff(den$x[1:2])
  dy <- diff(den$y[1:2])
  sz <- sort(den$z)
  c1 <- cumsum(sz) * dx * dy
  levels <- sapply(prob, function(x) {
    stats::approx(c1, sz, xout = 1 - x)$y
  })
  cl <- grDevices::contourLines(den$x, den$y, den$z, levels = levels)
  return(cl)
}

#' @title Rectangular posterior density of a node
#'
#' @description
#' Compute the posterior density of a set of node values under a Cauchy process on a phylogenetic tree.
#' 
#' @param tree the phylogenetic tree.
#' @param fit_lat a \code{\link{fitCauchy}} object for the latitude trait.
#' @param fit_long a \code{\link{fitCauchy}} object for the longitude trait.
#' @param anc_lat an object of class \code{ancestralCauchy}, obtained with \code{\link{ancestral}} on \code{fit_lat}.
#' @param anc_long an object of class \code{ancestralCauchy}, obtained with \code{\link{ancestral}} on \code{fit_long}.
#' @param level the level for the HPD intervals
#' 
#' @return An object of class \code{\link[tidytree]{treedata-class}}, that can
#' be exported with function \code{\link[treeio]{write.beast}}.
#' 
#' @seealso \code{\link{ancestral}}, \code{\link{fitCauchy}}
#' 
#' @author Paul Bastide \email{paul.bastide@m4x.org} and Gilles Didier \email{gilles.didier@free.fr}
#' 
#' @keywords internal
#' 
#'
ancestral_to_treedata_rec <- function(tree, 
                                      fit_lat, fit_long, 
                                      anc_lat, anc_long,
                                      level = 0.95) {
  
  if (!requireNamespace("HDInterval", quietly = TRUE) || 
      !requireNamespace("treeio", quietly = TRUE) || 
      !requireNamespace("tidytree", quietly = TRUE)) {
    stop(
      "Packages \"treeio\", \"tidytree\" and \"HDInterval\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  ## TODO: check entries
  
  n_nodes <- tree$Nnode
  n_tips <- length(tree$tip.label)
  
  ## Construct HPD intervals
  dens_lat <- anc_to_density(anc_lat)
  dens_long <- anc_to_density(anc_long)
  
  cred_lat <- lapply(dens_lat, function(x) HDInterval::hdi(x, credMass = sqrt(level), allowSplit = TRUE))
  cred_long <- lapply(dens_long, function(x) HDInterval::hdi(x, credMass = sqrt(level), allowSplit = TRUE))
  
  trait_lat <- vector(mode = "list", length = n_nodes)
  trait_long <- vector(mode = "list", length = n_nodes)
  

  ## Product of intervals in two dimensions
  for (i in seq_len(n_nodes)) {
    cc_lat <- cred_lat[[i]]
    cc_long <- cred_long[[i]]
    
    ncomp_lat <- nrow(cc_lat)
    ncomp_long <- nrow(cc_long)
    ncomp_tot <- ncomp_lat * ncomp_long
    
    tt_lat <- matrix(NA, ncol = 4, nrow = ncomp_tot)
    tt_long <- matrix(NA, ncol = 4, nrow = ncomp_tot)
    
    j <- 1
    for (la in 1:ncomp_lat) {
      for (lo in 1:ncomp_long) {
        tt_lat[j, ] <- c(cc_lat[la, ], rev(cc_lat[la, ]))
        tt_long[j, ] <- c(rep(cc_long[lo, ], each = 2)) 
        j <- j + 1
      }
    }
    
    trait_lat[[i]] <- tt_lat
    trait_long[[i]] <- tt_long
    
  }
  
  ## Create data list object
  rec_table <- list(node = seq_len(n_tips + n_nodes))

  max_comp <- max(sapply(trait_lat, nrow))
  trait_name_lat <- paste0("location1_", level * 100, "%HPD")
  trait_name_long <- paste0("location2_", level * 100, "%HPD")

  for (i in seq_len(max_comp)) {
    rec_table[[paste0(trait_name_lat, "_", i)]] <- c(rep(NA, n_tips), lapply(trait_lat, function(x) get_line_safe(x, i)))
    rec_table[[paste0(trait_name_long, "_", i)]] <- c(rep(NA, n_tips), lapply(trait_long,function(x) get_line_safe(x, i)))
  }

  rec_table[["location1"]] <- c(fit_lat$y, sapply(dens_lat, function(dd) stats::weighted.mean(dd$x, dd$y)))
  rec_table[["location2"]] <- c(fit_long$y, sapply(dens_long, function(dd) stats::weighted.mean(dd$x, dd$y)))
  
  rec_table[["location1_median"]] <- c(rep(NA, n_tips), sapply(dens_lat, function(dd) matrixStats::weightedMedian(dd$x, dd$y)))
  rec_table[["location2_median"]] <- c(rep(NA, n_tips), sapply(dens_long, function(dd) matrixStats::weightedMedian(dd$x, dd$y)))
  
  rec_table <- treeio::as_tibble(rec_table)
  
  ## Create treedata object
  tree_tibble <- treeio::as_tibble(tree)
  # tree_tibble$label <- rec_table$label
  
  tree_data <- treeio::full_join(tree_tibble, rec_table, by = 'node')
  
  return(treeio::as.treedata(tree_data))
  
}

#' @title Get line of a matrix
#'
#' @description
#' Get the line if the index is acceptable, otherwise return NA.
#' 
#' @param x a matrix
#' @param i the line number
#' 
#' @return the line number i of matrix x, or NA
#' 
#' @keywords internal
#' 
get_line_safe <- function(x, i) {
  if (i <= nrow(x)) return(x[i, ])
  return(NA)
}

#' @title Create densities from ancestral
#'
#' @description
#' Transform an object of class \code{ancestralCauchy}
#' into a list with \code{density} objects.
#' 
#' @param anc an object of class \code{ancestralCauchy}, obtained with \code{\link{ancestral}}.
#' 
#' @return list of matrching density objects
#' 
#' @keywords internal
#' 
anc_to_density <- function(anc) {
  x = as.numeric(colnames(anc))
  tmp_func <- function(aa) {
    dens <- list(x = x, y = aa)
    class(dens) <- "density"
    return(dens)
  }
  anc_dens <- apply(anc, 1, tmp_func)
  return(anc_dens)
}

#' @title Write treedata in evolaps format
#' 
#' @param treedata an object of class \code{\link[tidytree]{treedata-class}}
#' @param file the path of the file
#' 
#' @keywords internal
#' 
write.evolaps <- function(treedata, file) {
  treeio::write.beast(treedata, file = file, tree.name = "TREE1")
  tree_nex <- readLines(file)
  tree_nex <- tree_nex[!grepl("[R-package treeio", tree_nex, fixed = TRUE)]
  tree_nex <- sub("BEGIN TAXA", "Begin taxa", tree_nex)
  tree_nex <- sub("DIMENSIONS NTAX", "Dimensions ntax", tree_nex)
  tree_nex <- sub("TAXLABELS", "Taxlabels", tree_nex)
  tree_nex <- sub("END", "End", tree_nex)
  tree_nex <- sub("BEGIN TREES;", "Begin trees;", tree_nex)
  tree_nex <- sub("TRANSLATE", "Translate", tree_nex)
  tree_nex <- sub("\tTREE *", "tree ", tree_nex, fixed = TRUE)
  # tree_nex <- sub("[&R] ", "[&R] \n", tree_nex, fixed = TRUE)
  write(tree_nex, file = file)
}

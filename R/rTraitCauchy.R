#' @importFrom stats reorder

#' @title Cauchy Trait Simulation
#'
#' @description
#' Simulate a continuous trait using the Cauchy Process
#' 
#' @param n number of independent replicates
#' @param phy a phylogeny in \code{\link{ape}} \code{\link[ape]{phylo}} format.
#' @param model a phylogenetic model. Default is "cauchy", for the Cauchy process. Alternative are "lambda", "kappa", and "delta".
#' @param parameters list of parameters for the model (see Details).
#' 
#' @return
#' If n=1, a numeric vector with names from the tip labels in the tree.
#' For more than 1 replicate, a matrix with the tip labels as row names, and one column per replicate.
#' 
#' @details
#' The default choice of parameters is as follow:
#' \describe{
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
#' }
#' 
#' @seealso \code{\link[phylolm]{rTrait}}, \code{\link[ape]{rTraitCont}}
#' 
#' @examples
#' set.seed(1289)
#' phy <- ape::rphylo(40, 0.01, 0)
#' # One trait
#' y <- rTraitCauchy(n = 1, phy = phy, model = "cauchy",
#'                   parameters = list(root.value = 0, disp = 0.1))
#' y
#' plot(phy, x.lim = c(0, 750))
#' phydataplot(y, phy, offset = 150)
#' # Many trait
#' y <- rTraitCauchy(n = 10, phy = phy, model = "cauchy",
#'                   parameters = list(root.value = 0, disp = 0.1))
#' head(y)
#' 
#' 
#' @export
#' 
rTraitCauchy <- function(n = 1, phy,
                         model = c("cauchy", "lambda", "kappa", "delta"),
                         parameters = NULL) {
  ## Check parameters
  if (is.null(n) || length(n) > 1) stop("n needs to be an integer (number of replicates)")
  n = as.numeric(n)
  check_tree(phy)
  phy <- reorder(phy, "pruningwise")
  ## Model and parameters
  model = match.arg(model)
  parameters.default = c(0, 1, 1, 1, 1)
  names(parameters.default) <- c("root.value", "disp", 
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
                    is.null(parameters$lambda), is.null(parameters$kappa), is.null(parameters$delta))
    parameters.user <- c(parameters$root.value, parameters$disp,
                         parameters$lambda, parameters$kappa, parameters$delta)
    parameters <- parameters.default
    parameters[specified] <- parameters.user
  }
  p <- as.list(parameters)
  if (model == "lambda" & p$lambda == 1) model <- "cauchy"
  if (model == "kappa" & p$kappa == 1) model <- "cauchy"
  if (model == "delta" & p$delta == 1) model <- "cauchy"
  # Simulations
  if (model %in% c("lambda", "kappa", "delta")) {
    phy <- transf.branch.lengths(phy, model, parameters = p, check.pruningwise = F)$tree
  }
  sim <- vapply(seq_len(n),
                FUN = function(i) simulateTipsCauchy(tree = phy, root.value = p$root.value, disp = p$disp),
                FUN.VALUE = rep(0.0, length(phy$tip.label)))
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
#' @param root.value the initial root trait value.
#' @param disp the dispersion parameter of the Cauchy process.
#' 
#' @return a vector of simulated values.
#' 
#' 
#' @keywords internal
#' 
#' 
simulateTipsCauchy <- function(tree, root.value, disp) {
	res <-.Call("SimulateTipsCauchy", tree, root.value, disp)
	res
}
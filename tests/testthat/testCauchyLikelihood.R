test_that("testLikelihoodTwoTipsTree", {
  
  ## cauchy dist
  cauchy <- function(x, m, d) {
    1/pi * d / ((x - m)^2 + d^2)
  }

  ## Parameters
  mu <- 1
  disp <- 0.1
  
  ## tree with two tips
  tree <- read.tree(text = "(A:1,B:0.5);")
  
  ## data
  trait <- c(0.5, 3)
  names(trait) <- c("A", "B")
  
  ## likelihood manual
  lmanA <- function(x0) 1/pi * (disp * tree$edge.length[1]) / ((trait[1] - x0)^2 + (disp * tree$edge.length[1])^2)
  lmanB <- function(x0) 1/pi * (disp * tree$edge.length[2]) / ((trait[2] - x0)^2 + (disp * tree$edge.length[2])^2)
  lman <- unname(lmanA(mu) * lmanB(mu))
  
  ## likelihood
  lalgolse <- logDensityTipsCauchy(tree, trait, mu, disp, method = "fixed.root")

  ## equal ?
  expect_equal(log(lman), lalgolse)
  
  ## Random root
  root.edge <- 10
  tree$root.edge <- root.edge
  lmanroot <- function(x) cauchy(x, 0, disp * root.edge)
  lmanfun <- Vectorize(function(x) {lmanA(x) * lmanB(x) * lmanroot(x)})
  lmanRandomRoot <- unname(integrate(lmanfun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  
  lalgolse <- logDensityTipsCauchy(tree, trait, 0, disp, method = "random.root")
  expect_equal(log(lmanRandomRoot), lalgolse)
})

test_that("testLikelihoodThreeTipsTree", {
  
  ## Parameters
  mu <- 0.8
  disp <- 0.6
  
  ## tree with three tips
  tree <- read.tree(text = "((A:1,B:0.5):0.3, C:1);")
  
  ## data
  trait <- c(0.4, 1.1, -1.2)
  names(trait) <- c("A", "B", "C")
  
  ## cauchy dist
  cauchy <- function(x, m, d) {
    1/pi * d / ((x - m)^2 + d^2)
  }
  
  ## likelihood manual (numerical integral)
  lmanC <- function(x0) cauchy(trait[3], x0, disp * tree$edge.length[4])
  
  lmanA <- function(x) cauchy(x, trait[1], disp * tree$edge.length[2])
  lmanB <- function(x) cauchy(x, trait[2], disp * tree$edge.length[3])
  lmanLat <- function(x, x0) cauchy(x, x0, disp * tree$edge.length[1])
  lmanAB <- function(x0) {
    lint <- function(x) lmanA(x) * lmanB(x) * lmanLat(x, x0)
    res <- integrate(lint, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)
    return(unname(res$value))
  }
  
  lman <- unname(lmanAB(mu) * lmanC(mu))
  
  ## likelihood
  lalgolse <- logDensityTipsCauchy(tree, trait, mu, disp, method = "fixed.root")
  
  ## equal ?
  expect_equal(log(lman), lalgolse)
  
  ## Random root
  root.edge <- 10
  tree$root.edge <- root.edge
  lmanroot <- function(x) cauchy(x, 0, disp * root.edge)
  lmanfun <- Vectorize(function(x) {lmanAB(x) * lmanC(x) * lmanroot(x)})
  lmanRandomRoot <- unname(integrate(lmanfun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  
  lalgolse <- logDensityTipsCauchy(tree, trait, 0, disp, method = "random.root")

  expect_equal(log(lmanRandomRoot), lalgolse)

})

test_that("testLikelihoodThreeTipsTreeUltra", {
  
  ## Parameters
  mu <- 0.8
  disp <- 0.6
  
  ## tree with three tips
  tree <- read.tree(text = "((A:1.5,B:1.5):0.3, C:1.8);")
  
  ## data
  trait <- c(0.4, 1.1, -1.2)
  names(trait) <- c("A", "B", "C")
  
  ## cauchy dist
  cauchy <- function(x, m, d) {
    1/pi * d / ((x - m)^2 + d^2)
  }
  
  ## likelihood manual (numerical integral)
  lmanC <- function(x0) cauchy(trait[3], x0, disp * tree$edge.length[4])
  
  lmanA <- function(x) cauchy(x, trait[1], disp * tree$edge.length[2])
  lmanB <- function(x) cauchy(x, trait[2], disp * tree$edge.length[3])
  lmanLat <- function(x, x0) cauchy(x, x0, disp * tree$edge.length[1])
  lmanAB <- function(x0) {
    lint <- function(x) lmanA(x) * lmanB(x) * lmanLat(x, x0)
    res <- integrate(lint, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)
    return(unname(res$value))
  }
  
  lman <- unname(lmanAB(mu) * lmanC(mu))
  
  ## likelihood
  lalgolse <- logDensityTipsCauchy(tree, trait, mu, disp, method = "fixed.root")
  
  ## equal ?
  expect_equal(log(lman), lalgolse)
  
  ## Random root
  root.edge <- 10
  tree$root.edge <- root.edge
  lmanroot <- function(x) cauchy(x, 0, disp * root.edge)
  lmanfun <- Vectorize(function(x) {lmanAB(x) * lmanC(x) * lmanroot(x)})
  lmanRandomRoot <- unname(integrate(lmanfun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  
  lalgolse <- logDensityTipsCauchy(tree, trait, 0, disp, method = "random.root")

  expect_equal(log(lmanRandomRoot), lalgolse)
  
})

test_that("testLikelihoodLSENormalisation", {

  ## Parameters
  mu <- 0.8
  disp <- 0.1

  ## tree with three tips
  set.seed(1289)
  n <- 100
  tree <- rphylo(n, 0.1, 0)
  tree_height <- max(diag(vcv(tree)))

  ## data
  trait <- rTraitCont(tree, model = "BM", sigma = 1)

  ## likelihood
  lalgolse <- logDensityTipsCauchy(tree, trait, mu, disp, method = "fixed.root")
  lalgolsenorm <- -length(tree$tip.label) * log(disp) + logDensityTipsCauchy(tree, trait / disp, mu / disp, 1, method = "fixed.root")
  lalgolsenorm2 <- -length(tree$tip.label) * log(disp/tree_height) + logDensityTipsCauchy(tree, trait / (disp/tree_height), mu / (disp/tree_height), disp /(disp/tree_height), method = "fixed.root")

  ## equal ?
  expect_equal(lalgolsenorm, lalgolse, tolerance = 1e-7)
  expect_equal(lalgolsenorm2, lalgolse, tolerance = 1e-7)

  # ## Plot
  # s <- seq(0.01, 0.5, 0.01)
  # densLSE <- sapply(s, function(x) logDensityTipsCauchy(tree, trait, mu, x, ultrametric = FALSE, lse = TRUE))
  # dens <- sapply(s, function(x) logDensityTipsCauchy(tree, trait, mu, x, ultrametric = FALSE, lse = FALSE))
  # densNorm <- sapply(s, function(x) -length(tree$tip.label) * log(x) + logDensityTipsCauchy(tree, trait / x, mu / x, 1, ultrametric = FALSE, lse = TRUE))
  # densNorm2 <- sapply(s, function(x) -length(tree$tip.label) * log(x/tree_height) + logDensityTipsCauchy(tree, trait / (x/tree_height), mu / (x/tree_height), x/(x/tree_height), ultrametric = FALSE, lse = TRUE))
  # plot(s, dens)
  # points(s, densLSE, pch = "+")
  # points(s, densNorm, pch = "-")
  # points(s, densNorm2, pch = "*")

})

test_that("testLikelihoodLSEBig", {

  ## Parameters
  mu <- 0.8
  disp <- 0.05

  ## tree with three tips
  set.seed(1289)
  n <- 1000
  tree <- rphylo(n, 0.1, 0)
  tree_height <- max(diag(vcv(tree)))

  ## data
  trait <- simulateTipsCauchy(tree, mu, disp)

  ## likelihood
  lalgolse <- logDensityTipsCauchy(tree, trait, mu, disp, method = "fixed.root")
  lalgolsenorm <- -length(tree$tip.label) * log(disp) + logDensityTipsCauchy(tree, trait / disp, mu / disp, 1,  method = "fixed.root")
  lalgolsenorm2 <- -length(tree$tip.label) * log(disp/tree_height) + logDensityTipsCauchy(tree, trait / (disp/tree_height), mu / (disp/tree_height), disp /(disp/tree_height),  method = "fixed.root")

  ## equal ?
  expect_equal(lalgolsenorm, lalgolse, tolerance = 1e-2)
  expect_equal(lalgolsenorm2, lalgolse, tolerance = 1e-2)

  # ## Plot
  # s <- seq(0.01, 0.5, 0.01)
  # densLSE <- sapply(s, function(x) logDensityTipsCauchy(tree, trait, mu, x, method = "fixed.root"))
  # densNorm <- sapply(s, function(x) -length(tree$tip.label) * log(x) + logDensityTipsCauchy(tree, trait / x, mu / x, 1,  method = "fixed.root"))
  # densNorm2 <- sapply(s, function(x) -length(tree$tip.label) * log(x/tree_height) + logDensityTipsCauchy(tree, trait / (x/tree_height), mu / (x/tree_height), x/(x/tree_height),  method = "fixed.root"))
  # # plot(s, dens)
  # plot(s, densLSE, pch = "+")
  # points(s, densNorm, pch = "-")
  # points(s, densNorm2, pch = "*")

})

test_that("testLikelihoodFunction", {
  
  ## Parameters
  mu <- 0.8
  disp <- 0.05
  
  ## tree with three tips
  set.seed(1289)
  n <- 50
  tree <- rphylo(n, 0.1, 0)
  
  root.edge <- 1
  treerr <- tree
  treerr$root.edge <- root.edge
  
  ## data
  trait <- simulateTipsCauchy(tree, mu, disp)
  
  ## Errors
  expect_error(logDensityTipsCauchy(tree, trait,root.value = NULL, disp = disp, method = "fixed.root"),
               "Starting value must be specified for root node in the `fixed.root` method.")
  expect_error(logDensityTipsCauchy(tree, trait,root.value = NULL, disp = disp, method = "random.root"),
               "Starting value must be specified for root node in the `random.root` method.")
  expect_error(logDensityTipsCauchy(tree, trait, root.value = 0.0, disp = disp, method = "random.root"),
               "In the random root model, the `root.edge` must be non NULL and non zero.")
  expect_error(logDensityTipsCauchy(tree, trait, root.value = 0.0, disp = disp, method = "random.root"),
               "In the random root model, the `root.edge` must be non NULL and non zero.")
  expect_error(logDensityTipsCauchy(tree, trait, root.value = mu, disp = disp, method = "reml"),
               "In the reml model, `root.value` cannot be specified.")
  
  ## Equalities
  expect_equal(logDensityTipsCauchy(treerr, trait, root.value = mu, disp = disp, method = "random.root"),
               logDensityTipsCauchy(treerr, trait - mu, root.value = 0.0, disp = disp, method = "random.root"))
  expect_equal(logDensityTipsCauchy(tree, trait, root.value = mu, disp = disp, method = "fixed.root"),
               logDensityTipsCauchy(tree, trait - mu, root.value = 0.0, disp = disp, method = "fixed.root"))
  
  ## Non binary tree
  tree$edge.length[sample(seq_len(nrow(tree$edge)), 10)] <- 0
  ff <- fitCauchy(tree, trait)
  expect_equal(ff$disp, 0.06785236, tolerance = 1e-6)
  tree_bis <- di2multi(tree)
  expect_error(check_binary_tree(tree_bis), "The tree must be binary.")
  
  ## Tree with no branch length
  tree <- stree(10)
  expect_error(check_binary_tree(tree), "the tree has no branch lengths.")
  expect_error(check_tree(tree), "the tree has no branch lengths.")
  
  ## Tree with no labels
  tree <- rtree(10)
  tree$tip.label <- NULL
  expect_error(check_binary_tree(tree), "the tree has no tip labels.")
  expect_error(check_tree(tree), "the tree has no tip labels.")
  
  ## Tree with no labels
  tree <- rtree(10)
  class(tree) <- NULL
  expect_error(check_binary_tree(tree), "tree object is not of class \"phylo\"")
  expect_error(check_tree(tree), "tree object is not of class \"phylo\"")
  
})


# test_that("cauphylm error", {
#   ## Phylogenetic Tree
#   data(lizards)
#   tree <- lizards$phy
#   # tree$edge.length <- tree$edge.length / max(vcv(tree))
# 
#   ###########################
#   ## Data Simulation
#   sigma <- 0.1
#   mu <- 0
# 
#   set.seed(1289)
#   datBM <- rTraitCont(tree, sigma = sigma, root.value = mu)
# 
#   ###########################
#   ## Likelihood
#   # Instability
#   expect_equal(logDensityTipsCauchy(tree, datBM, 0, 0.01, method = "fixed.root"),
#                logDensityTipsCauchy(tree, datBM, 0, exp(log(0.01)), method = "fixed.root"), tolerance = 1e-6)
# 
#   ## Plot
#   s <- seq(0.001, 0.05, 0.001)
#   densLSE <- sapply(s, function(x) logDensityTipsCauchy(tree, datBM, 0, x, method = "fixed.root"))
#   densLSENorm <- sapply(s, function(x) -length(tree$tip.label) * log(x) + logDensityTipsCauchy(tree, datBM / x, 0 / x, 1, method = "fixed.root"))
#   plot(s, densLSE, pch = "+")
#   points(s, densLSENorm, pch = "-")
# 
#   ## normalization
#   tree_height <- max(vcv(tree))
#   expect_equal(
#     -length(tree$tip.label) * log(exp(log(0.1)) * tree_height) + log(densityTipsCauchy(tree, datBM / (exp(log(0.1)) * tree_height), 0, 1 / tree_height)),
#     -length(tree$tip.label) * log(0.1 * tree_height) + log(densityTipsCauchy(tree, datBM / (0.1 * tree_height), 0, 1 / tree_height))
#     )
# 
#   expect_equal(
#     densityTipsCauchy(tree, datBM / (exp(log(0.1)) * tree_height), 0, 1 / tree_height),
#     densityTipsCauchy(tree, datBM / (0.1 * tree_height), 0, 1 / tree_height)
#     )
# 
#   (0.1 / tree_height)^(-length(tree$tip.label)) * densityTipsCauchy(tree, datBM * tree_height, 0, 0.1 * tree_height)
#   (exp(log(0.1)) / tree_height)^(-length(tree$tip.label)) * densityTipsCauchy(tree, datBM * tree_height, 0, exp(log(0.1)) * tree_height)
# 
#   (0.1 / tree_height)^(-length(tree$tip.label)) * densityTipsCauchy(tree, datBM / (0.1 / tree_height), 0, tree_height)
#   (exp(log(0.1)) / tree_height)^(-length(tree$tip.label)) * densityTipsCauchy(tree, datBM / (exp(log(0.1)) / tree_height), 0, tree_height)
# 
#   densityTipsCauchy(tree, datBM / 0.1, 0, 1)
#   densityTipsCauchy(tree, datBM / exp(log(0.1)), 0, 1)
# 
#   logDensityTipsCauchy(tree, datBM, 0, 0.1)
#   logDensityTipsCauchy(tree, datBM, 0, exp(log(0.1)))
# 
# })
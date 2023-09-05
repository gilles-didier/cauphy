test_that("testTwoTipsTree", {
  
  ## Parameters
  disp <- 0.1
  
  ## tree with two tips
  tree <- read.tree(text = "(A:1,B:0.5);")
  
  ## data
  trait <- c(0.5, 3)
  names(trait) <- c("A", "B")
  
  ## AIR - random root
  root.edge <- 1
  treerr <- tree
  treerr$root.edge <- root.edge
  
  expect_equal(sum(posteriorDensityIncrement(3, 20:25, treerr, trait, root.value = 0, disp = disp, method = "random.root")), 0.0)
  
  fun <- function(x) posteriorDensityIncrement(3, x, treerr, trait, root.value = 0, disp = disp, method = "random.root")
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  treerr$root.edge <- 0.01
  fun <- function(x) -posteriorDensityIncrement(3, x, treerr, trait, root.value = 0, disp = disp, method = "random.root")
  MAP <- optim(0.0, fun, method = "Brent", lower = -5, upper = 5)
  expect_equal(MAP$par,  0.5, tolerance = 1e-1)
  
})

test_that("testLikelihoodThreeTipsTree", {
  
  ## Parameters
  mu <- -2
  disp <- 0.6
  
  ## tree with three tips
  tree <- read.tree(text = "((A:1,B:0.5):0.3, C:1);")
  
  ## data
  trait <- c(0.4, 1.1, -1.2)
  names(trait) <- c("A", "B", "C")
  
  ## Random root
  root.edge <- 1
  treerr <- tree
  treerr$root.edge <- root.edge
  
  ## Errors
  expect_error(posteriorDensityIncrement(7, 1, tree, trait, root.value = mu, disp = disp, method = "fixed.root"),
               "This node does not exist in the tree.")
  expect_error(posteriorDensityIncrement(4, 1, tree, trait,root.value = mu, disp = disp, method = "fixed.root"),
               "Ancestral increment reconstruction is not allowed for the root branch with the fixed root model.")
  expect_error(posteriorDensityIncrement(4, 1, tree, trait,root.value = mu, disp = disp, method = "reml"),
               "Ancestral increment reconstruction is not allowed for the root branch with the reml model.")
  # expect_error(posteriorDensityIncrement(5, 1, tree, trait,root.value = mu, disp = disp, method = "reml"),
  #              "Ancestral increment reconstruction is not allowed for branches connected to the root with the reml model.")
  expect_error(posteriorDensityIncrement(5, 1, tree, trait,root.value = NULL, disp = disp, method = "fixed.root"),
               "Starting value must be specified for root node in the `fixed.root` method.")
  expect_error(posteriorDensityIncrement(5, 1, tree, trait,root.value = NULL, disp = disp, method = "random.root"),
               "Starting value must be specified for root node in the `random.root` method.")
  expect_error(posteriorDensityIncrement(4, 1, tree, trait, root.value = 0.0, disp = disp, method = "random.root"),
               "In the random root model, the `root.edge` must be non NULL and non zero.")
  expect_error(posteriorDensityIncrement(4, 1, tree, trait, root.value = 0.0, disp = disp, method = "random.root"),
               "In the random root model, the `root.edge` must be non NULL and non zero.")
  expect_error(posteriorDensityIncrement(2, 1, tree, trait, root.value = mu, disp = disp, method = "reml"),
               "In the reml model, `root.value` cannot be specified.")
  
  ## Fixed root
  expect_equal(posteriorDensityIncrement(5, 20, tree, trait, root.value = mu, disp = disp, method = "fixed.root"), 0.0)
  expect_equal(posteriorDensityIncrement(1, 1, tree, trait, root.value = mu, disp = disp, method = "fixed.root"), 0.06531257, tolerance = 1e-3)
  expect_equal(posteriorDensityIncrement(2, 1, tree, trait, root.value = mu, disp = disp, method = "fixed.root"), 0.2420315, tolerance = 1e-3)
  # dirac
  expect_warning(expect_equal(
    posteriorDensityIncrement(3, c(-10, 0.75, 1, trait[3] - mu), tree, trait, root.value = mu, disp = disp, method = "fixed.root"),
    c(0, 0, 0, 1)),
    "This branch ends at a tip, and the root is fixed: the posterior increment density is a Dirac in 0.8.")

  # Sum to 1
  fun <- Vectorize(function(x) {posteriorDensityIncrement(5, x, tree, trait, root.value = mu, disp = disp, method = "fixed.root")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- Vectorize(function(x) {posteriorDensityIncrement(2, x, tree, trait, root.value = mu, disp = disp, method = "fixed.root")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- Vectorize(function(x) {posteriorDensityIncrement(1, x, tree, trait, root.value = mu, disp = disp, method = "fixed.root")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)

  # MAP
  fun <- Vectorize(function(x) {-posteriorDensityIncrement(5, x, tree, trait, root.value = mu, disp = disp, method = "fixed.root")})
  MAP <- optim(1.0, fun, method = "Brent", lower = -10, upper = 10)
  expect_equal(MAP$par,  3.0, tolerance = 1e-2)

  ## Random root
  expect_equal(posteriorDensityIncrement(5, 25, treerr, trait, root.value = mu, disp = disp, method = "random.root"), 0.0)

  fun <- Vectorize(function(x) {posteriorDensityIncrement(5, x, treerr, trait, root.value = mu, disp = disp, method = "random.root")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- Vectorize(function(x) {posteriorDensityIncrement(4, x, treerr, trait, root.value = mu, disp = disp, method = "random.root")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- Vectorize(function(x) {posteriorDensityIncrement(3, x, treerr, trait, root.value = mu, disp = disp, method = "random.root")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- Vectorize(function(x) {posteriorDensityIncrement(2, x, treerr, trait, root.value = mu, disp = disp, method = "random.root")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- Vectorize(function(x) {posteriorDensityIncrement(1, x, treerr, trait, root.value = mu, disp = disp, method = "random.root")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)

  treerr$root.edge <- 0.001
  fun <- Vectorize(function(x) {-posteriorDensityIncrement(4, x, treerr, trait, root.value = mu, disp = disp, method = "random.root")})
  MAP <- optim(1.0, fun, method = "Brent", lower = -10, upper = 10)
  expect_equal(MAP$par,  0.0, tolerance = 1e-4)

  ## reml
  expect_equal(posteriorDensityIncrement(2, 20, tree, trait, root.value = NULL, disp = disp, method = "reml"), 0.0, tolerance = 1e-5)
  
  expect_equal(posteriorDensityIncrement(2, 0.4, tree, trait, root.value = NULL, disp = disp, method = "reml"),
               posteriorDensityIncrement(2, 0.4, reroottip(tree, 3), trait[-3], root.value = trait[3], disp = disp, method = "random.root"),
               tolerance = 1e-5)
  
  expect_equal(posteriorDensityIncrement(2, 0.4, tree, trait, root.value = NULL, disp = disp, method = "reml"),
               posteriorDensityIncrement(1, 0.4, reroottip(tree, 1), trait[-1], root.value = trait[1], disp = disp, method = "random.root"),
               tolerance = 1e-5)
  
  expect_equal(posteriorDensityIncrement(1, 0.4, tree, trait, root.value = NULL, disp = disp, method = "reml"),
               posteriorDensityIncrement(1, 0.4, reroottip(tree, 3), trait[-3], root.value = trait[3], disp = disp, method = "random.root"),
               tolerance = 1e-5)
  
  expect_equal(posteriorDensityIncrement(5, c(0.1, 0.4), tree, trait, root.value = NULL, disp = disp, method = "reml"),
               c(posteriorDensityIncrement(5, 0.1, tree, trait, root.value = NULL, disp = disp, method = "reml"),
                 posteriorDensityIncrement(5, 0.4, tree, trait, root.value = NULL, disp = disp, method = "reml")),
               tolerance = 1e-5)

  fun <- Vectorize(function(x) {posteriorDensityIncrement(2, x, tree, trait, root.value = NULL, disp = disp, method = "reml")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- Vectorize(function(x) {posteriorDensityIncrement(1, x, tree, trait, root.value = NULL, disp = disp, method = "reml")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- Vectorize(function(x) {posteriorDensityIncrement(3, x, tree, trait, root.value = NULL, disp = disp, method = "reml")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- Vectorize(function(x) {posteriorDensityIncrement(5, x, tree, trait, root.value = NULL, disp = disp, method = "reml")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
})

test_that("testAncestralIncrementFit", {

  ## Parameters
  mu <- 0.0
  disp <- 0.6

  ## tree with three tips
  set.seed(1289)
  n <- 10
  tree <- rphylo(n, 0.1, 0)

  ## data
  trait <- rTraitCont(tree, model = "BM", sigma = 1)

  ## Random root
  root.edge <- 1
  treerr <- tree
  treerr$root.edge <- root.edge

  ## ASR - fixed root
  fitfr <- fitCauchy(tree, trait, method = "fixed.root")
  expect_equal(unname(as.vector(t(increment(fitfr, c(12, 13), c(-0.3, 0.1))))),
               c(posteriorDensityIncrement(12, c(-0.3, 0.1), tree, trait, root.value = fitfr$x0, disp = fitfr$disp, method = "fixed.root"),
                 posteriorDensityIncrement(13, c(-0.3, 0.1), tree, trait, root.value = fitfr$x0, disp = fitfr$disp, method = "fixed.root")))

  anc_all <- increment(fitfr)
  expect_equal(dim(anc_all), c(18, 100))
  expect_equal(anc_all, increment(cauphylm(trait ~ 1, phy = tree)))
  
  fun <- function(x) increment(fitfr, 12, x)
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- function(x) increment(fitfr, 1, x)
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- function(x) increment(fitfr, 2, x)
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)

  fun <- function(x) -increment(fitfr, 12, x)
  MAP <- optim(1.0, fun, method = "Brent", lower = -10, upper = 10)
  expect_equal(MAP$par,  0.2811777, tolerance = 1e-4)

  fitlm <- cauphylm(trait ~ 1, phy = tree)
  expect_equal(increment(fitlm, 12, 2), increment(fitfr, 12, 2))

  ## ASR - random root
  fitrr <- fitCauchy(tree, trait, method = "random.root", root.edge = root.edge)
  expect_equal(unname(as.vector(t(increment(fitrr, c(12, 13), c(-0.3, 0.1))))),
               c(posteriorDensityIncrement(12, c(-0.3, 0.1), treerr, trait, root.value = 0.0, disp = fitrr$disp, method = "random.root"),
                 posteriorDensityIncrement(13, c(-0.3, 0.1), treerr, trait, root.value = 0.0, disp = fitrr$disp, method = "random.root")))
  
  anc_all <- increment(fitrr)
  expect_equal(dim(anc_all), c(19, 100))
  
  treerr$root.edge <- 0.00001
  fun <- function(x) -increment(fitrr, 11, x)
  MAP <- optim(1.0, fun, method = "Brent", lower = -10, upper = 10)
  expect_equal(MAP$par,  mu, tolerance = 1e-2)

  ## ASR - reml
  fitreml <- fitCauchy(tree, trait, method = "reml")
  expect_equal(unname(as.vector(t(increment(fitreml, c(14, 1, 10), c(-0.3, 0.1, 2))))),
               c(posteriorDensityIncrement(14, c(-0.3, 0.1, 2), tree, trait, root.value = fitreml$x0, disp = fitreml$disp, method = "reml"),
                 posteriorDensityIncrement(1, c(-0.3, 0.1, 2), tree, trait, root.value = fitreml$x0, disp = fitreml$disp, method = "reml"),
                 posteriorDensityIncrement(10, c(-0.3, 0.1, 2), tree, trait, root.value = fitreml$x0, disp = fitreml$disp, method = "reml")))

  anc_all <- increment(fitreml)
  expect_equal(dim(anc_all), c(18, 100))
  
  fun <- Vectorize(function(x) {posteriorDensityIncrement(10, x, tree, trait, root.value = NULL, disp = disp, method = "reml")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- Vectorize(function(x) {posteriorDensityIncrement(13, x, tree, trait, root.value = NULL, disp = disp, method = "reml")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  ## Error
  expect_error(ancestral(fitfr, c(0.2, 0.4)), "whole numbers")
  expect_error(ancestral(cauphylm(trait ~ 1, phy = tree), c(0.2, 0.4)), "whole numbers")
})

test_that("testHDI", {
  
  library(HDInterval)
  
  ## Parameters
  mu <- 0.0
  disp <- 0.6
  
  ## tree with three tips
  set.seed(1289)
  n <- 10
  tree <- rphylo(n, 0.1, 0)
  
  ## data
  trait <- rTraitCont(tree, model = "BM", sigma = 1)
  
  ## Random root
  root.edge <- 1
  treerr <- tree
  treerr$root.edge <- root.edge
  
  ## ASR - fixed root
  fitfr <- fitCauchy(tree, trait, method = "fixed.root")
  anc_all <- increment(fitfr)
  hdi_all <- hdi(anc_all)
  expect_equal(dim(hdi_all[[1]]), c(2, 2))
  expect_equal(dim(hdi_all[[2]]), c(1, 2))
  expect_equal(hdi_all[[7]], hdi(anc_all, node = 7)[[1]])
  
  anc_all <- increment(fitfr, n_values = 2000)
  cm <- 0.4
  hdi_all <- hdi(anc_all, credMass = cm)
  fun <- function(x) increment(fitfr, 12, x)
  total_mass <- unname(integrate(fun, hdi_all[["12"]][1, 1], hdi_all[["12"]][1, 2],
                                 rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, cm, tolerance = 1e-2)
  
  cm <- 0.8
  hdi_all <- hdi(anc_all, credMass = cm)
  fun <- function(x) increment(fitfr, 1, x)
  total_mass <- unname(integrate(fun, hdi_all[["1"]][1, 1], hdi_all[["1"]][1, 2],
                                 rel.tol = .Machine$double.eps^0.5)$value)
  total_mass <- total_mass + unname(integrate(fun, hdi_all[["1"]][2, 1], hdi_all[["1"]][2, 2],
                                              rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, cm, tolerance = 1e-2)
  
  fun <- function(x) increment(fitfr, 2, x)
  total_mass <- unname(integrate(fun, hdi_all[["2"]][1, 1], hdi_all[["2"]][1, 2],
                                 rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, cm, tolerance = 1e-2)
  
  ## Errors
  expect_error(hdi(anc_all, node = 0.5), "must be whole numbers.")
  expect_message(hdi(anc_all, node = c(2, 37)), "Nodes 37 are not in the ancestralCauchy reconstruction object. They will be ignored.")
  expect_message(expect_error(hdi(anc_all, node = c(37)), "There are no node left."),
                              "Nodes 37 are not in the ancestralCauchy reconstruction object. They will be ignored.")
  
  ## Plot
  expect_message(plot(anc_all, intervals = hdi_all, node = c(2, 37)),
                 "Nodes 37 are not in the ancestralCauchy reconstruction object. They will be ignored.")
})
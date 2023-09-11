test_that("testLikelihoodTwoTipsTree", {
  
  ## Parameters
  disp <- 0.1
  
  ## tree with two tips
  tree <- read.tree(text = "(A:1,B:0.5);")
  
  ## data
  trait <- c(0.5, 3)
  names(trait) <- c("A", "B")
  
  ## ASR - random root
  root.edge <- 1
  treerr <- tree
  treerr$root.edge <- root.edge
  
  expect_equal(sum(posteriorDensityAncestral(3, 20:25, treerr, trait, root.value = 0, disp = disp, method = "random.root")), 0.0)
  
  fun <- function(x) posteriorDensityAncestral(3, x, treerr, trait, root.value = 0, disp = disp, method = "random.root")
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  treerr$root.edge <- 0.01
  fun <- function(x) -posteriorDensityAncestral(3, x, treerr, trait, root.value = 0, disp = disp, method = "random.root")
  MAP <- optim(0.0, fun, method = "Brent", lower = -0.1, upper = 0.1)
  expect_equal(MAP$par,  0.0, tolerance = 1e-4)
  
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
  
  ## Random root
  root.edge <- 1
  treerr <- tree
  treerr$root.edge <- root.edge

  ## ASR - Errors
  expect_error(posteriorDensityAncestral(1, 1, tree, trait, root.value = mu, disp = disp, method = "fixed.root"),
               "Ancestral reconstruction is only allowed for ancestral nodes.")
  expect_error(posteriorDensityAncestral(7, 1, tree, trait, root.value = mu, disp = disp, method = "fixed.root"),
               "This node does not exist in the tree.")
  expect_error(posteriorDensityAncestral(4, 1, tree, trait,root.value = mu, disp = disp, method = "fixed.root"),
               "Ancestral state reconstruction is not allowed for the root with the fixed root model.")
  expect_error(posteriorDensityAncestral(5, 1, tree, trait,root.value = NULL, disp = disp, method = "fixed.root"),
               "Starting value must be specified for root node in the `fixed.root` method.")
  expect_error(posteriorDensityAncestral(5, 1, tree, trait,root.value = NULL, disp = disp, method = "random.root"),
               "Starting value must be specified for root node in the `random.root` method.")
  expect_error(posteriorDensityAncestral(4, 1, tree, trait, root.value = 0.0, disp = disp, method = "random.root"),
               "In the random root model, the `root.edge` must be non NULL and non zero.")
  expect_error(posteriorDensityAncestral(4, 1, tree, trait, root.value = 0.0, disp = disp, method = "random.root"),
               "In the random root model, the `root.edge` must be non NULL and non zero.")
  expect_error(posteriorDensityAncestral(4, 1, tree, trait, root.value = mu, disp = disp, method = "reml"),
               "In the reml model, `root.value` cannot be specified.")
  
  ## ASR - fixed root
  expect_equal(posteriorDensityAncestral(5, 20, tree, trait, root.value = mu, disp = disp, method = "fixed.root"), 0.0)
  
  fun <- Vectorize(function(x) {posteriorDensityAncestral(5, x, tree, trait, root.value = mu, disp = disp, method = "fixed.root")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- Vectorize(function(x) {-posteriorDensityAncestral(5, x, tree, trait, root.value = mu, disp = disp, method = "fixed.root")})
  MAP <- optim(1.0, fun, method = "Brent", lower = -10, upper = 10)
  expect_equal(MAP$par,  0.8289036, tolerance = 1e-4)
  
  ## ASR - random root
  expect_equal(posteriorDensityAncestral(5, 20, treerr, trait, root.value = mu, disp = disp, method = "random.root"), 0.0)
  
  fun <- Vectorize(function(x) {posteriorDensityAncestral(5, x, treerr, trait, root.value = mu, disp = disp, method = "random.root")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  treerr$root.edge <- 0.001
  fun <- Vectorize(function(x) {-posteriorDensityAncestral(4, x, treerr, trait, root.value = mu, disp = disp, method = "random.root")})
  MAP <- optim(1.0, fun, method = "Brent", lower = -10, upper = 10)
  expect_equal(MAP$par,  mu, tolerance = 1e-4)
  
  treerr$root.edge <- 100000
  fun <- Vectorize(function(x) {-posteriorDensityAncestral(4, x, treerr, trait, root.value = 0, disp = disp, method = "random.root")})
  MAP_rr <- optim(1.0, fun, method = "Brent", lower = -10, upper = 10)
  
  ## ASR - reml
  expect_equal(posteriorDensityAncestral(5, 20, tree, trait, root.value = NULL, disp = disp, method = "reml"), 0.0, tolerance = 1e-5)
  
  fun <- Vectorize(function(x) {posteriorDensityAncestral(4, x, tree, trait, root.value = NULL, disp = disp, method = "reml")})
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- Vectorize(function(x) {-posteriorDensityAncestral(4, x, tree, trait, root.value = NULL, disp = disp, method = "reml")})
  MAP_reml <- optim(1.0, fun, method = "Brent", lower = -10, upper = 10)
  
  expect_equal(MAP_rr$par, MAP_reml$par, tolerance = 1e-3)
  
})

test_that("testAncestralFit", {
  
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
  expect_equal(unname(as.vector(t(ancestral(fitfr, c(12, 13, 15), c(-0.3, 0.1))))),
               c(posteriorDensityAncestral(12, c(-0.3, 0.1), tree, trait, root.value = fitfr$x0, disp = fitfr$disp, method = "fixed.root"),
                 posteriorDensityAncestral(13, c(-0.3, 0.1), tree, trait, root.value = fitfr$x0, disp = fitfr$disp, method = "fixed.root"),
                 posteriorDensityAncestral(15, c(-0.3, 0.1), tree, trait, root.value = fitfr$x0, disp = fitfr$disp, method = "fixed.root")))
  
  anc_all <- ancestral(fitfr)
  expect_equal(dim(anc_all), c(8, 100))
  expect_equal(anc_all, ancestral(cauphylm(trait ~ 1, phy = tree)))
  
  fun <- function(xx) ancestral(fitfr, 12, xx)
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- function(xx) ancestral(fitfr, 15, xx)
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- function(x) -ancestral(fitfr, 12, x)
  MAP <- optim(1.0, fun, method = "Brent", lower = -10, upper = 10)
  expect_equal(MAP$par,  -0.5086788, tolerance = 1e-4)
  
  fitlm <- cauphylm(trait ~ 1, phy = tree)
  expect_equal(ancestral(fitlm, 12, 2), ancestral(fitfr, 12, 2))
  
  ## ASR - random root
  fitrr <- fitCauchy(tree, trait, method = "random.root", root.edge = root.edge)
  expect_equal(unname(as.vector(t(ancestral(fitrr, c(12, 13), c(-0.3, 0.1))))),
               c(posteriorDensityAncestral(12, c(-0.3, 0.1), treerr, trait, root.value = 0.0, disp = fitrr$disp, method = "random.root"),
                 posteriorDensityAncestral(13, c(-0.3, 0.1), treerr, trait, root.value = 0.0, disp = fitrr$disp, method = "random.root")))
  
  fun <- function(xx) ancestral(fitrr, 12, xx)
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- function(xx) ancestral(fitrr, 15, xx)
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  anc_all <- ancestral(fitrr)
  expect_equal(dim(anc_all), c(9, 100))
  
  treerr$root.edge <- 0.00001
  fun <- function(x) -ancestral(fitrr, 11, x)
  MAP <- optim(1.0, fun, method = "Brent", lower = -10, upper = 10)
  expect_equal(MAP$par,  mu, tolerance = 1e-2)
  
  treerr$root.edge <- 100000
  fun <- function(x) -ancestral(fitrr, 11, x)
  MAP_rr <- optim(1.0, fun, method = "Brent", lower = -10, upper = 10)
  
  ## ASR - reml
  fitreml <- fitCauchy(tree, trait, method = "reml")
  expect_equal(unname(as.vector(t(ancestral(fitreml, c(12, 13), c(-0.3, 0.1))))),
               c(posteriorDensityAncestral(12, c(-0.3, 0.1), tree, trait, root.value = fitreml$x0, disp = fitreml$disp, method = "reml"),
                 posteriorDensityAncestral(13, c(-0.3, 0.1), tree, trait, root.value = fitreml$x0, disp = fitreml$disp, method = "reml")))
  
  anc_all <- ancestral(fitreml)
  expect_equal(dim(anc_all), c(9, 100))
  
  fun <- function(xx) ancestral(fitreml, 15, xx)
  total_mass <- unname(integrate(fun, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value)
  expect_equal(total_mass, 1.0)
  
  fun <- function(x) -ancestral(fitrr, 11, x)
  MAP_reml <- optim(1.0, fun, method = "Brent", lower = -10, upper = 10)
  
  expect_equal(MAP_rr$par, MAP_reml$par, tolerance = 1e-3)
  
  ## Error
  expect_error(ancestral(fitfr, c(0.2, 0.4)), "whole numbers")
  expect_error(ancestral(cauphylm(trait ~ 1, phy = tree), c(0.2, 0.4)), "whole numbers")
  
  ## lambda
  fit <- fitCauchy(tree, trait, method = "reml", model = "lambda")
  expect_error(ancestral(fit), "Ancestral reconstruction is only available for the Cauchy process.")
  expect_error(increment(fit), "Ancestral reconstruction is only available for the Cauchy process.")
  fit <- cauphylm(trait ~ 1, phy = tree, model = "lambda")
  expect_error(ancestral(fit), "Ancestral reconstruction is only available for the Cauchy process.")
  expect_error(increment(fit), "Ancestral reconstruction is only available for the Cauchy process.")
  
  ## lm
  reg <- ape::rTraitCont(tree, model = "BM", sigma = 0.1, root.value = 0)
  fit <- cauphylm(trait ~ reg, phy = tree)
  expect_error(ancestral(fit), "Ancestral reconstruction is only available for the Cauchy regression agains the intercept.")
  expect_error(increment(fit), "Ancestral reconstruction is only available for the Cauchy regression agains the intercept.")
  
})

test_that("testAncestralFitPlot", {
  
  ## Parameters
  mu <- 0.0
  sigma <- 1
  
  ## tree with three tips
  set.seed(1289)
  n <- 10
  tree <- rphylo(n, 0.1, 0)
  
  ## data
  trait <- rTraitCont(tree, model = "BM", sigma = sigma, root.value = 0)
  
  ## ASR - fixed root
  fitfr <- fitCauchy(tree, trait, method = "fixed.root")
  anc_fr <- ancestral(fitfr)
  inc_fr <- increment(fitfr)
  plot_asr(fitfr, anc = anc_fr)
  plot_asr(fitfr, inc = inc_fr)
  plot_asr(fitfr, inc = inc_fr, common_colorscale = TRUE)
  plot_asr(fitfr, anc = anc_fr, inc = inc_fr, width.edge = 1, scaling = 0.01)
  plot_asr(fitfr, anc = anc_fr, inc = inc_fr, width.edge = 1, scaling = 0.01, x.intersp = 4)
  plot_asr(fitfr, anc = anc_fr, inc = inc_fr, width.edge = 1, scaling = 10)
  plot_asr(fitfr, anc = anc_fr, inc = inc_fr, width.edge = 1, common_colorscale = TRUE)
  
  ## ASR - reml
  fitreml <- fitCauchy(tree, trait, method = "reml")
  anc_reml <- ancestral(fitreml)
  inc_reml <- increment(fitreml)
  plot_asr(fitreml, anc = anc_reml)
  plot_asr(fitreml, inc = inc_reml)
  plot_asr(fitreml, anc = anc_reml, inc = inc_reml)

  expect_message(plot(anc_fr, node = c(9, 11, 12)), 
                 "Nodes 9, 11 are not in the ancestralCauchy reconstruction object.")
  expect_message(expect_error(plot(anc_fr, node = c(9, 11)), "no node left"))
  expect_error(plot(anc_fr, node = c(9.2)), "whole number")
  expect_error(plot_asr(fitreml, anc = 3), "'anc' must by of S3 class 'ancestralCauchy'")
  expect_error(plot_asr(fitreml, inc = 3), "'inc' must by of S3 class 'ancestralCauchy'")
  
})
test_that("testREMLThreeTipsTree", {
  
  ## Parameters
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
  lint <- Vectorize(function(mu){
    unname(lmanAB(mu) * lmanC(mu))
  })
  # Compute the restricted likelihood
  lman <- integrate(lint, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)$value
  
  ## REML
  lalgolse <- logDensityTipsCauchy(tree, trait, NULL, disp, method = "reml")
  
  ## equal ?
  expect_equal(log(lman), lalgolse)
  
  ## Other tip as root
  retree <- reroottip(tree, 3)
  lalgolse <- logDensityTipsCauchy(retree, trait[-3], trait[3], disp, method = "random.root")
  expect_equal(log(lman), lalgolse)
  
})


test_that("testREMLFourTipsTree", {
  
  ## Parameters
  disp <- 1
  
  ## tree with three tips
  tree <- read.tree(text = "(((A:1.3,B:1.3):0.3, C:1.6):0.8, D:2.4);")
  expect_true(is.ultrametric(tree))
  
  ## data
  trait <- c(3.9, 5.9, 1.9, 2.9)
  names(trait) <- c("A", "B", "C", "D")
  
  ## cauchy dist
  cauchy <- function(x, m, d) {
    1/pi * d / ((x - m)^2 + d^2)
  }
  
  ## likelihood manual (numerical integral)
  lmanD <- function(x0) cauchy(trait[4], x0, disp * tree$edge.length[6])
  lmanC <- function(x0) cauchy(trait[3], x0, disp * tree$edge.length[5])
  lmanA <- function(x) cauchy(x, trait[1], disp * tree$edge.length[3])
  lmanB <- function(x) cauchy(x, trait[2], disp * tree$edge.length[4])
  lmanLat <- function(x, x0) cauchy(x, x0, disp * tree$edge.length[2])
  lmanAB <- Vectorize(function(x0) {
    lint <- function(x) lmanA(x) * lmanB(x) * lmanLat(x, x0)
    res <- integrate(lint, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)
    return(unname(res$value))
  }
  )
  lmanBranch <- function(x, x0) cauchy(x, x0, disp * tree$edge.length[1])
  lmanABC <- function(x0) {
    lint <- function(x) lmanAB(x) * lmanC(x) * lmanBranch(x, x0)
    res <- integrate(lint, -Inf, +Inf, rel.tol = .Machine$double.eps^0.5)
    return(unname(res$value))
  }
  lfin <- Vectorize(function(mu){
    unname(lmanABC(mu) * lmanD(mu))
  })
  # Compute the restricted likelihood
  lman <- integrate(lfin, -Inf, +Inf, rel.tol = 0.1)$value
  
  ## REML
  lalgolse <- logDensityTipsCauchy(tree, trait, NULL, disp, method = "reml")
  
  ## equal ?
  expect_equal(log(lman), lalgolse, tolerance = 1e-4)
  
  ## Other tips as root
  lalgolse2 <- logDensityTipsCauchy(reroottip(tree, 3), trait[-3], trait[3], disp, method = "random.root")
  expect_equal(lalgolse, lalgolse2)
  lalgolse2 <- logDensityTipsCauchy(reroottip(tree, 4), trait[-4], trait[4], disp, method = "random.root")
  expect_equal(lalgolse, lalgolse2)
  
})

test_that("testRootingStrategies", {
  
  ## Parameters
  disp <- 0.1
  rootTip <- 1
  
  ## tree with three tips
  set.seed(1289)
  n <- 100
  tree <- rphylo(n, 0.1, 0)
  # tree <- stree(n, type = "balanced")
  # tree$edge.length <- rep(0.1, nrow(tree$edge))

  ## data
  trait <- rTraitCont(tree, model = "BM", sigma = 1)
  
  ## likelihood
  ll1 <- logDensityTipsCauchy(tree, trait, NULL, disp, method = "reml")
  
  # Reroot
  tip <- rootTip
  tipName <- tree$tip.label[rootTip]
  retree <- reroottip(tree, tip)
  # root value
  root.value <- trait[tipName]
  tipTraitBis <- trait[names(trait) != tipName]
  # Force integers
  stopifnot(all.equal(matrix(as.integer(retree$edge), ncol = 2), retree$edge))
  retree$edge <- matrix(as.integer(retree$edge), ncol = 2)
  # likelihood
  ll2 <- logDensityTipsCauchy(tree = retree, tipTrait = tipTraitBis, root.value = root.value, disp = disp, method = "random.root")
  # likelihood
  ll3 <- logDensityTipsCauchy(tree = retree, tipTrait = tipTraitBis - root.value, root.value = 0.0, disp = disp, method = "random.root")
  
  ## equal ?
  expect_equal(ll1, ll2)
  expect_equal(ll1, ll3)
  
  ## disp 0
  ll1 <- logDensityTipsCauchy(tree, trait, NULL, 0, method = "reml", rootTip = rootTip)
  ll2 <- logDensityTipsCauchy(tree = retree, tipTrait = tipTraitBis, root.value = root.value, disp = 0, method = "random.root")
  expect_equal(ll1, ll2, tolerance = 1e-5)
  
  ## disp 1.3
  ll1 <- logDensityTipsCauchy(tree, trait, NULL, 1.3, method = "reml", rootTip = rootTip)
  ll2 <- logDensityTipsCauchy(tree = retree, tipTrait = tipTraitBis, root.value = root.value, disp = 1.3, method = "random.root")
  expect_equal(ll1, ll2, tolerance = 1e-5)
  
  # dd <- seq(0.1, 2, 0.01)
  # ll1 <- sapply(dd, function(d) logDensityTipsCauchy(tree, trait, NULL, d, method = "reml", rootTip = rootTip))
  # ll2 <- sapply(dd, function(d) logDensityTipsCauchy(tree = retree, tipTrait = tipTraitBis, root.value = root.value, disp = d, method = "random.root"))
  # plot(dd, ll1)
  # points(dd, ll2, col = "red")
})

test_that("testREMLError", {
  tree <- read.tree(text = "(((((t7:0.4421316155,t6:0.4421316155):0.3741413588,((t4:0.3812102035,t8:0.3812102035):0.1468429299,t3:0.5280531334):0.2882198409):0.04669787571,t5:0.8629708501):0.02858242796,(t9:0.2390721964,t2:0.2390721964):0.6524810816):0.108446722,(t10:0.1825693151,t1:0.1825693151):0.8174306849);")
  trait <- 1:10 / 10
  names(trait) <- tree$tip.label
  
  expect_equal(logDensityTipsCauchy(tree, trait, disp = 0.1, method = "reml"),
               logDensityTipsCauchy(reroottip(tree, 5), trait[-5], trait[5], disp = 0.1, method = "random.root"))
  
  expect_equal(logDensityTipsCauchy(tree, trait, disp = 0.1),
               logDensityTipsCauchy(reroottip(tree, 10), trait[-10], trait[10], disp = 0.1, method = "random.root"))

})
test_that("testLikelihoodMulti", {
  
  ## Parameters
  n <- 20
  p <- 4
  
  ## tree
  set.seed(1289)
  tree <- rphylo(n, 0.1, 0)
  tree_height <- max(diag(vcv(tree)))
  
  # data
  mu <- 1.3
  disp <- 1.0
  dat <- rTraitCauchy(n = p, phy = tree, model = "cauchy", parameters = list(root.value = mu, disp = disp))
  
  # transformation
  A <- matrix(rnorm(p^2), p, p)
  dispVec <- rlnorm(p)
  Id <- diag(dispVec)
  M <- rep(mu, p)
  
  # test transformations
  expect_equal(sweep(dat, 2, M) %*% t(solve(A)), t(solve(A, t(dat) - M)))
  expect_equal(t(t(dat) - M) %*% t(solve(A)), t(solve(A, t(dat) - M)))
  
  # microbenchmark::microbenchmark(sweep(dat, 2, M) %*% t(solve(A)),
  #                                t(t(dat) - M) %*% t(solve(A)),
  #                                t(solve(A, t(dat) - M)),
  #                                solve(A) %*% (t(dat) - M),
  #                                solve(A, t(dat) - M))
  
  # microbenchmark::microbenchmark(dat %*% t(solve(A)),
  #                                solve(A, t(dat)))
  
  ## Fixed Root
  # likelihood
  ll1 <- logDensityTipsCauchyMulti(tree, dat, root.value = M, transMat = A, method = "fixed.root")
  ll2 <- logDensityTipsCauchyMulti(tree, dat, root.value = M, transMat = Id, method = "fixed.root")
  expect_true(ll1 <= ll2)
  
  # independent
  llind <- 0
  for (k in 1:p) {
    llind <- llind + logDensityTipsCauchy(tree, dat[, k], root.value = M[k], disp = dispVec[k], method = "fixed.root")
  }
  llind2 <- 0
  for (k in 1:p) {
    llind2 <- llind2 + logDensityTipsCauchy(tree, dat[, k] - M[k], root.value = 0.0, disp = dispVec[k], method = "fixed.root")
  }
  expect_equal(ll2, llind)
  expect_equal(llind, llind2)
  
  ## REML
  # likelihood
  ll1 <- logDensityTipsCauchyMulti(tree, dat, transMat = A, method = "reml")
  ll2 <- logDensityTipsCauchyMulti(tree, dat, transMat = Id, method = "reml")
  expect_true(ll1 <= ll2)
  
  # independent
  llind <- 0
  for (k in 1:p) {
    llind <- llind + logDensityTipsCauchy(tree, dat[, k], disp = dispVec[k], method = "reml")
  }
  expect_equal(ll2, llind)
  
  ## random root
  # likelihood
  tree$root.edge <- 10
  ll1 <- logDensityTipsCauchyMulti(tree, dat, root.value = M[k], transMat = A, method = "random.root")
  ll2 <- logDensityTipsCauchyMulti(tree, dat, root.value = M[k], transMat = Id, method = "random.root")
  expect_true(ll1 <= ll2)
  
  # independent
  llind <- 0
  for (k in 1:p) {
    llind <- llind + logDensityTipsCauchy(tree, dat[, k], root.value = M[k], disp = dispVec[k], method = "random.root")
  }
  expect_equal(ll2, llind)

})

test_that("testLikelihoodBi", {
  
  ## Parameters
  n <- 50
  p <- 2
  
  ## tree
  set.seed(1289)
  tree <- rphylo(n, 0.1, 0)
  tree_height <- max(diag(vcv(tree)))
  
  # data
  mu <- 1.3
  disp <- 1.0
  dat <- rTraitCauchy(n = p, phy = tree, model = "cauchy", parameters = list(root.value = mu, disp = disp))
  
  # transformation
  theta <- c(pi / 3, pi / 3 + pi / 4)
  dispVec <- c(0.1, 2.3)
  A <- sapply(theta, function(x) c(cos(x), sin(x)))
  A <- A %*% diag(dispVec)
  M <- rep(mu, p)
  
  ## Fixed Root
  ll1 <- logDensityTipsCauchyMulti(tree, dat, root.value = M, transMat = A, method = "fixed.root")
  ll2 <- logDensityTipsCauchyBi(tree, dat, root.value = M, angle = theta, disp = dispVec, method = "fixed.root")
  expect_equal(ll1, ll2)
  
  ## REML
  ll1 <- logDensityTipsCauchyMulti(tree, dat, transMat = A, method = "reml")
  ll2 <- logDensityTipsCauchyBi(tree, dat, angle = theta, disp = dispVec, method = "reml")
  expect_equal(ll1, ll2)
  
  ## random root
  tree$root.edge <- 10
  ll1 <- logDensityTipsCauchyMulti(tree, dat, root.value = M, transMat = A, method = "random.root")
  ll2 <- logDensityTipsCauchyBi(tree, dat, root.value = M, angle = theta, disp = dispVec, method = "random.root")
  expect_equal(ll1, ll2)
  
})
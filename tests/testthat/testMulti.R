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

# test_that("testLikelihoodSurfaceBi", {
# 
#   ## Parameters
#   n <- 50
#   p <- 2
# 
#   ## tree
#   set.seed(1289)
#   tree <- rphylo(n, 0.1, 0)
#   tree_height <- max(diag(vcv(tree)))
# 
#   # True parameters
#   theta <- c(pi / 3, pi / 3 + pi / 4)
#   dispVec <- c(0.1, 2.3)
#   A <- sapply(theta, function(x) c(cos(x), sin(x)))
#   A <- A %*% diag(dispVec)
#   mu <- 0.0
#   M <- rep(mu, p)
# 
#   # data
#   dat <- rTraitCauchy(n = p, phy = tree, model = "cauchy", parameters = list(root.value = 0.0, disp = 1.0))
#   dat <- dat %*% t(A)
# 
#   ## REML, fixed disp
#   # grid of angle values
#   lim_angle <- c(0, pi)
#   ang <- seq(lim_angle[1], lim_angle[2], pi / 50)
#   ang_vals <- expand.grid(ang, ang)
#   ang_vals <- ang_vals[((ang_vals[, 2] - ang_vals[, 1]) > pi / 100) & (abs(ang_vals[, 1] - ang_vals[, 2]) < pi), ]
# 
#   method <- "reml"
#   densang <- apply(ang_vals, 1, function(aa) logDensityTipsCauchyBi(tree, dat, angle = aa, disp = dispVec, method = method))
#   data_plot <- data.frame(theta1 = ang_vals[, 1], theta2 = ang_vals[, 2], dens = densang)
#   ang_hat <- ang_vals[which.max(densang), ]
# 
#   library(ggplot2)
#   ggplot(data_plot, aes(x = theta1, y = theta2, z = dens)) +
#     geom_contour_filled(binwidth = 10) +
#     theme_bw() + coord_fixed(xlim = lim_angle, ylim = lim_angle) +
#     geom_point(aes(x = theta[1], y = theta[2]))
#   # zoom
#   ggplot(subset(data_plot, theta2 > 1.75 & theta2 < 1.9 & theta1 < 1.2), aes(x = theta1, y = theta2, z = dens)) +
#     geom_contour_filled(binwidth = 5) +
#     theme_bw() + coord_cartesian(xlim = c(0, 1.2), ylim = c(1.75, 1.9)) +
#     geom_point(aes(x = theta[1], y = theta[2])) +
#     geom_point(aes(x = ang_hat[, 1], y = ang_hat[, 2]), color = "blue")
# 
#   ## REML, fixed angles
#   # grid of disp values
#   dd_vals <- expand.grid(seq(0.05, 0.5, 0.01), seq(1, 3, 0.1))
#   dd_vals <- dd_vals[dd_vals[, 2] > dd_vals[, 1], ]
#   method <- "reml"
#   densgrid <- apply(dd_vals, 1, function(dd) logDensityTipsCauchyBi(tree, dat, angle = theta, disp = dd, method = method))
#   data_plot <- data.frame(disp1 = dd_vals[, 1], disp2 = dd_vals[, 2], dens = densgrid)
#   disp_hat <- dd_vals[which.max(densgrid), ]
# 
#   library(ggplot2)
#   ggplot(data_plot, aes(x = disp1, y = disp2, z = dens)) +
#     geom_contour_filled(binwidth = 1) +
#     theme_bw() +
#     geom_point(aes(x = dispVec[1], y = dispVec[2])) +
#     geom_point(aes(x = disp_hat[, 1], y = disp_hat[, 2]), color = "blue")
# 
#   ## Init
#   # SVD
#   svd_init <- svd(dat)
#   theta_init_1 <- atan2(svd_init$v[2, 1], svd_init$v[1, 1])
#   theta_init_2 <- (theta_init_1 + pi / 2) %% pi
#   theta_init <- sort(c(theta_init_1, theta_init_2))
#   theta_init
#   theta
#   ang_hat
# 
#   # disp
#   rot_mat <- sapply(theta_init, function(x) c(cos(x), sin(x)))
#   dat_trans <- dat %*% rot_mat
#   t(dat_trans[, 2]) %*% dat_trans[, 1]
#   initDispersionParameter(tree, dat_trans[, 1], 0.0, method.init.disp = "Qn", method)
#   initDispersionParameter(tree, dat_trans[, 2], 0.0, method.init.disp = "Qn", method)
#   dispVec
#   disp_hat
# 
# })

test_that("testTransform", {
  
  ## log sorted transform
  ori <- c(2.3, 0.1)
  trans <- log_sorted_transform(ori)
  back <- log_sorted_back_transform(trans)
  expect_equal(ori, back)
  # derivative
  dd <- gradient_sorted_transform(ori)
  ddnum <- nloptr::nl.jacobian(ori, sorted_transform)
  expect_equal(dd, ddnum)
  # derivative
  dd <- gradient_log_sorted_transform(ori)
  ddnum <- nloptr::nl.jacobian(ori, log_sorted_transform)
  expect_equal(dd, ddnum)
  
  ## logit transform
  ori <- c(pi / 3, pi / 3 + pi / 4)
  trans <- logit_transform(ori, 0, pi)
  back <- logit_back_transform(trans, 0, pi)
  expect_equal(ori, back)
  # derivative
  dd <- gradient_logit_transform(ori, lower_bound = 0, upper_bound = pi)
  ddnum <- nloptr::nl.jacobian(ori, logit_transform, lower_bound = 0, upper_bound = pi)
  expect_equal(dd, diag(ddnum))
  
  ## All transforms
  ori <- c(-3.4, 2.5, pi / 3, pi / 3 + pi / 4, 2.3, 0.1)
  names(ori) <- c("coef1", "coef2", "angle1", "angle2", "disp1", "disp2")
  trans <- transform_values(ori)
  back <- back_transform_values(trans)
  expect_equal(ori, back)
  # derivative
  dd <- gradient_transform_values(ori)
  ddnum <- nloptr::nl.jacobian(ori, transform_values)
  expect_equal(unname(dd), ddnum)
  
})

test_that("testFitBi", {

  ## Parameters
  n <- 20
  p <- 2
  ## tree
  set.seed(1289)
  tree <- rphylo(n, 0.1, 0)
  tree_height <- max(diag(vcv(tree)))
  # True parameters
  theta <- c(pi / 3, pi / 3 + pi / 4)
  dispVec <- c(0.1, 2.3)
  A <- sapply(theta, function(x) c(cos(x), sin(x)))
  A <- A %*% diag(dispVec)
  mu <- 0.0
  M <- rep(mu, p)
  # data
  dat <- rTraitCauchy(n = p, phy = tree, model = "cauchy", parameters = list(root.value = 0.0, disp = 1.0))
  dat <- dat %*% t(A)

  ## Fixed Root
  fitfr <- fitCauchyBi(tree, dat, method = "fixed.root")
  fitfr
  # Compare with true value
  expect_true(logDensityTipsCauchyBi(tree, dat, M, dispVec, theta, method = "fixed.root") <= fitfr$logLik)
  # Compare with independent fits
  fitfr1 <- fitCauchy(tree, dat[, 1], method = "fixed.root")
  fitfr2 <- fitCauchy(tree, dat[, 2], method = "fixed.root")
  expect_true(fitfr$logLik >= fitfr1$logLik + fitfr2$logLik)
  expect_equal(1 - pchisq(fitfr$logLik - (fitfr1$logLik + fitfr2$logLik), df = 2), 0.0)
  # vcov
  fitfr <- compute_vcov(fitfr)
  expect_equal(dim(fitfr$vcov), c(6, 6))
  expect_equal(colnames(fitfr$vcov), c("x01", "x02", "angle1", "angle2", "disp1", "disp2" ))
  expect_equal(fitfr$vcov[6,6], 0.007512049, tolerance = 1e-4)
  expect_equal(unname(diag(fitfr$vcov)),
               c(2.353479e+01, 3.348212e+02, 8.433166e-08, 9.369576e-04, 3.633829e-01, 7.512049e-03),
               tolerance = 1e-5)
  # confint
  expect_message(ii <- confint(fitfr))
  expect_equal(rownames(ii), rownames(fitfr$vcov))
  expect_equal(ii[5,1], 0.7054907, tolerance = 1e-6)
  for (i in 1:6) {
    expect_true(ii[i, 1] <= fitfr$all_params[i] && ii[i, 2] >= fitfr$all_params[i])
  }
  # profile
  pr <- profile(fitfr)
  # plot(pr)
  expect_equal(max(pr$x01$profLogLik), fitfr$logLik)
  expect_equal(max(pr$x02$profLogLik), fitfr$logLik)
  expect_equal(max(pr$angle1$profLogLik), fitfr$logLik, tolerance = 1e-5)
  expect_equal(max(pr$angle2$profLogLik), fitfr$logLik, tolerance = 1e-2)
  expect_equal(max(pr$disp1$profLogLik), fitfr$logLik, tolerance = 1e-4)
  expect_equal(max(pr$disp2$profLogLik), fitfr$logLik, tolerance = 1e-5)
  
  
  ## REML
  fitreml <- fitCauchyBi(tree, dat, method = "reml")
  # Compare with true value
  expect_true(logDensityTipsCauchyBi(tree, dat, M, dispVec, theta, method = "reml") <= fitreml$logLik)
  # Compare with independent fits
  fitreml1 <- fitCauchy(tree, dat[, 1], method = "reml")
  fitreml2 <- fitCauchy(tree, dat[, 2], method = "reml")
  expect_true(fitreml$logLik >= fitreml1$logLik + fitreml1$logLik)
  expect_equal(1 - pchisq(fitreml$logLik - (fitreml1$logLik + fitreml2$logLik), df = 2), 0.0)
  # vcov
  fitreml <- compute_vcov(fitreml)
  expect_equal(dim(fitreml$vcov), c(4, 4))
  expect_equal(colnames(fitreml$vcov), c("angle1", "angle2", "disp1", "disp2" ))
  expect_equal(unname(diag(fitreml$vcov)),
               c(9.282913e-08, 8.612623e-04, 4.293138e-01, 8.710660e-03),
               tolerance = 1e-5)
  # confint
  expect_message(ii <- confint(fitreml))
  expect_equal(rownames(ii), rownames(fitreml$vcov))
  expect_equal(ii[4,1], 0.06478452, tolerance = 1e-6)
  for (i in 1:4) {
    expect_true(ii[i, 1] <= fitreml$all_params[i] && ii[i, 2] >= fitreml$all_params[i])
  }
  # profile
  pr <- profile(fitreml)
  # plot(pr)
  expect_equal(max(pr$angle1$profLogLik), fitreml$logLik, tolerance = 1e-5)
  expect_equal(max(pr$angle2$profLogLik), fitreml$logLik, tolerance = 1e-2)
  expect_equal(max(pr$disp1$profLogLik), fitreml$logLik, tolerance = 1e-4)
  expect_equal(max(pr$disp2$profLogLik), fitreml$logLik, tolerance = 1e-5)
  
  ## random root
  root.edge <- 100
  fitrr <- fitCauchyBi(tree, dat, method = "random.root", root.edge = root.edge)
  # Compare with true value
  treerr <- res$phy
  treerr$root.edge <- root.edge
  expect_true(logDensityTipsCauchyBi(treerr, dat, M, dispVec, theta, method = "random.root") <= fitrr$logLik)
  # Compare with independent fits
  fitrr1 <- fitCauchy(tree, dat[, 1], method = "random.root", root.edge = root.edge)
  fitrr2 <- fitCauchy(tree, dat[, 2], method = "random.root", root.edge = root.edge)
  expect_true(fitrr$logLik >= fitrr1$logLik + fitrr1$logLik)
  expect_equal(1 - pchisq(fitrr$logLik - (fitrr1$logLik + fitrr2$logLik), df = 2), 0.0)
  # vcov
  fitrr <- compute_vcov(fitrr)
  expect_equal(dim(fitrr$vcov), c(4, 4))
  expect_equal(colnames(fitrr$vcov), c("angle1", "angle2", "disp1", "disp2" ))
  expect_equal(unname(diag(fitrr$vcov)),
               c(4.819142e-08, 1.555542e-03, 5.880968e-01, 2.360420e-03),
               tolerance = 1e-5)
  # confint
  expect_message(ii <- confint(fitrr))
  expect_equal(rownames(ii), rownames(fitrr$vcov))
  expect_equal(ii[4,1], 0.03607927, tolerance = 1e-6)
  for (i in 1:4) {
    expect_true(ii[i, 1] <= fitrr$all_params[i] && ii[i, 2] >= fitrr$all_params[i])
  }
  # profile
  pr <- profile(fitrr)
  # plot(pr)
  expect_equal(max(pr$angle1$profLogLik), fitrr$logLik, tolerance = 1e-3)
  expect_equal(max(pr$angle2$profLogLik), fitrr$logLik, tolerance = 1e-2)
  expect_equal(max(pr$disp1$profLogLik), fitrr$logLik, tolerance = 1e-2)
  expect_equal(max(pr$disp2$profLogLik), fitrr$logLik, tolerance = 1e-2)
  
})
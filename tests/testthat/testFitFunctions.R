test_that("optim", {
  set.seed(1289)
  
  ## parameters
  ntips <- 20
  rateTest <- 0.1
  
  ## Tree
  tree <- rphylo(ntips, 0.1, 0)
  tree$edge.length <-  tree$edge.length / max(vcv(tree))
  
  ## root.value
  mu <- 1
  disp <- 0.1
  
  ## Sim
  trait <- simulateTipsCauchy(tree, mu, disp)
  
  ## Opt
  res <- fitCauchy(tree, trait, method = "fixed.root")
  res
  
  ## Change init
  expect_equal(res$logLik,
               fitCauchy(tree, trait, method = "fixed.root", method.init.disp = "Sn")$logLik)
  expect_equal(res$logLik,
               fitCauchy(tree, trait, method = "fixed.root", method.init.disp = "MAD")$logLik)
  expect_equal(res$logLik,
               fitCauchy(tree, trait, method = "fixed.root", method.init.disp = "IQR")$logLik)
  
  ## Compare with true value
  expect_true(logDensityTipsCauchy(tree, trait, mu, disp, method = "fixed.root") <= res$logLik)
  
  ## Global
  resGlob <- fitCauchy(tree, trait, method = "fixed.root", optim = "global")
  expect_equal(res$logLik, resGlob$logLik)
  
  ## vcov
  expect_equal(colnames(vcov(res)), c("x0", "disp"))
  expect_equal(unname(sqrt(diag(vcov(res)))), c(0.03243854, 0.01808950), tolerance = 1e-2)
  
  ## lambda
  reslam <- fitCauchy(tree, trait, method = 'fixed.root', root.edge = root.edge, model = "lambda")
  expect_true(reslam$logLik >= res$logLik)
  reslamglob <- fitCauchy(tree, trait, method = 'fixed.root', root.edge = root.edge, model = "lambda", optim = "global")
  expect_equal(reslam$logLik, reslamglob$logLik)
  
  ## vcov
  expect_equal(colnames(vcov(reslam)), c("x0", "disp", "lambda"))
  expect_equal(unname(diag(vcov(reslam))), c(1.457096e-03, 6.289735e-04, 1.100595e-05), tolerance = 1e-2)
})

test_that("optim, random root", {
  set.seed(1289)
  
  ## parameters
  ntips <- 20
  rateTest <- 0.1
  
  ## Tree
  tree <- rphylo(ntips, 0.1, 0)
  tree$edge.length <-  tree$edge.length / max(vcv(tree))
  
  ## params
  disp <- 0.1
  root.edge <- 100
  
  ## Sim
  mu <- rcauchy(1, 0, disp * root.edge)
  trait <- simulateTipsCauchy(tree, mu, disp)
  
  ## Opt
  res <- fitCauchy(tree, trait, method = 'random.root', root.edge = root.edge)
  res
  
  ## Change init
  expect_equal(res$logLik,
               fitCauchy(tree, trait, method = "random.root", method.init.disp = "Sn")$logLik)
  expect_equal(res$logLik,
               fitCauchy(tree, trait, method = "random.root", method.init.disp = "IQR")$logLik)
  
  ## Compare with true value
  treerr <- res$phy
  treerr$root.edge <- root.edge
  expect_true(logDensityTipsCauchy(treerr, trait, 0, disp, method = 'random.root') <= res$logLik)
  
  ## Global
  resGlob <- fitCauchy(tree, trait, method = 'random.root', root.edge = root.edge, optim = "global")
  expect_equal(res$logLik, resGlob$logLik)
  
  ## vcov
  expect_equal(colnames(vcov(res)), c("disp"))
  expect_equal(unname(sqrt(diag(vcov(res)))), c(0.021), tolerance = 1e-1)
  
  ## lambda
  reslam <- fitCauchy(tree, trait, method = 'random.root', root.edge = root.edge, model = "lambda")
  expect_true(reslam$logLik > res$logLik)
  
  ## vcov
  expect_equal(colnames(vcov(reslam)), c("disp", "lambda"))
  expect_equal(unname(sqrt(diag(vcov(reslam)))), c(0.0285, 0.0415), tolerance = 1e-1)
  
  ## lambda fixed
  reslam2 <- fitCauchy(tree, trait, method = 'random.root', root.edge = root.edge, model = "lambda",
                       starting.value = list(lambda = 1), lower.bound = list(lambda = 1), upper.bound = list(lambda = 1))
  expect_equal(reslam2$lambda, 1.0, ignore_attr = TRUE)
  expect_equal(reslam2$logLik, res$logLik)
  
})

test_that("optim, reml", {
  set.seed(1289)
  
  ## parameters
  ntips <- 20
  rateTest <- 0.1
  
  ## Tree
  tree <- rphylo(ntips, 0.1, 0)
  tree$edge.length <-  tree$edge.length / max(vcv(tree))
  
  ## params
  disp <- 0.1
  root.edge <- 1
  
  ## Sim
  mu <- rcauchy(1, 0, disp * root.edge)
  trait <- simulateTipsCauchy(tree, mu, disp)
  
  ## Opt
  res <- fitCauchy(tree, trait, method = 'reml')
  res
  
  ## Change init
  expect_equal(res$logLik,
               fitCauchy(tree, trait, method = "reml", method.init.disp = "Sn")$logLik)
  expect_equal(res$logLik,
               fitCauchy(tree, trait, method = "reml", method.init.disp = "MAD")$logLik)
  expect_equal(res$logLik,
               fitCauchy(tree, trait, method = "reml", method.init.disp = "IQR")$logLik)
  
  ## Compare with true value
  expect_true(logDensityTipsCauchy(res$phy, trait, NULL, disp, method = "reml") <= res$logLik)
  
  # logDensityTipsCauchy(res$phy, trait, NULL, res$disp, method = "reml")
  # dd <- seq(0.05, 0.2, 0.005)
  # plot(dd,
  #      sapply(dd, function(d) logDensityTipsCauchy(res$phy, trait, NULL, d, method = "reml")))
  
  ## Global
  resGlob <- fitCauchy(tree, trait, method = 'reml', optim = "global")
  expect_equal(res$logLik, resGlob$logLik)
  
  ## vcov
  expect_equal(colnames(vcov(res)), c("disp"))
  expect_equal(unname(sqrt(diag(vcov(res)))), c(0.02016588), tolerance = 1e-2)
  
  ## lambda
  reslam <- fitCauchy(tree, trait, method = 'reml', model = "lambda")
  expect_true(reslam$logLik > res$logLik)
  
  ## vcov
  expect_equal(colnames(vcov(reslam)), c("disp", "lambda"))
  expect_equal(unname(sqrt(diag(vcov(reslam)))), c(0.0283, 0.0574), tolerance = 1e-2)
  
  # ## comparison with random root
  # root.edge <- 1000000000000
  # resran <- fitCauchy(tree, trait, method = 'random.root', root.edge = root.edge)
  # expect_equal(resran$disp, res$disp)
  
})

test_that("cauphylm", {
  set.seed(1289)
  
  ## parameters
  ntips <- 20
  rateTest <- 0.1
  
  ## Tree
  tree <- rphylo(ntips, 0.1, 0)
  tree$edge.length <-  tree$edge.length / max(vcv(tree))
  
  ## root.value
  mu <- 1
  disp <- 0.1
  
  ## Sim
  trait <- simulateTipsCauchy(tree, mu, disp)
  
  ## Reg
  reg <- ape::rTraitCont(tree, model = "BM", sigma = 0.1, root.value = 0)
  dat <- data.frame(tt = trait, rr = reg)
  
  ## Opt
  reslm <- cauphylm(trait ~ 1, phy = tree)
  reslmdat <- cauphylm(tt ~ 1, data = dat, phy = tree)
  res <- fitCauchy(tree, trait, method = "fixed.root")
  
  ## Compare fits
  expect_equal(res$logLik, reslm$logLik)
  expect_equal(res$aic, reslm$aic)
  expect_equal(res$x0, reslm$coefficients, ignore_attr = TRUE, tolerance = 1e-4)
  expect_equal(res$disp, reslm$disp, ignore_attr = TRUE, tolerance = 1e-3)
  
  expect_equal(reslmdat[!(names(reslmdat) %in% c("formula", "call"))],
               reslm[!(names(reslm) %in% c("formula", "call"))])
  
  ## Compare vcov
  expect_equal(unname(vcov(reslm)), unname(vcov(res)), tolerance = 1e-3)
  
  ## Regression
  trait <- -3 + 2 * reg + simulateTipsCauchy(tree, 0, disp)
  dat <- data.frame(trait = trait, reg = reg)
  
  reslm <- cauphylm(trait ~ reg, data = dat, phy = tree)
  reslmnull <- cauphylm(trait ~ 1, data = dat, phy = tree)
  
  expect_true(reslmnull$logLik <= reslm$logLik)
  
})

test_that("cauphylm lambda", {
  set.seed(1289)
  
  ## parameters
  ntips <- 20
  rateTest <- 0.1
  
  ## Tree
  tree <- rphylo(ntips, 0.1, 0)
  tree$edge.length <-  tree$edge.length / max(vcv(tree))
  tree_lambda <- phylolm::transf.branch.lengths(tree, model = "lambda", parameters = list(lambda = 0.6))$tree
  
  ## root.value
  mu <- 1
  disp <- 0.1
  
  ## Sim
  trait <- simulateTipsCauchy(tree_lambda, mu, disp)
  
  ## Reg
  reg <- ape::rTraitCont(tree, model = "BM", sigma = 0.1, root.value = 0)
  dat <- data.frame(tt = trait, rr = reg)
  
  ## Opt
  reslm <- cauphylm(trait ~ 1, phy = tree, model = "lambda")
  reslmdat <- cauphylm(tt ~ 1, data = dat, phy = tree, model = "lambda")
  res <- fitCauchy(tree, trait, method = "fixed.root", model = "lambda")
  
  ## Compare fits
  expect_equal(reslmdat$logLik, reslm$logLik)
  expect_equal(reslmdat$coefficient, reslm$coefficients, ignore_attr = TRUE)
  expect_equal(reslmdat$disp, reslm$disp, ignore_attr = TRUE)
  expect_equal(reslmdat$lambda, reslm$lambda, ignore_attr = TRUE)
  
  expect_equal(res$logLik, reslm$logLik)
  expect_equal(res$aic, reslm$aic)
  expect_equal(res$x0, reslm$coefficients, ignore_attr = TRUE, tolerance = 1e-4)
  expect_equal(res$disp, reslm$disp, ignore_attr = TRUE, tolerance = 1e-3)
  expect_equal(res$lambda, reslm$lambda, ignore_attr = TRUE, tolerance = 1e-3)
  
  expect_equal(reslmdat[!(names(reslmdat) %in% c("formula", "call"))],
               reslm[!(names(reslm) %in% c("formula", "call"))])
  
  ## Compare without lambda
  reslmcau <- cauphylm(trait ~ 1, phy = tree, model = "cauchy")
  expect_true(reslmcau$logLik <= reslm$logLik)
  
  ## Compare vcov
  expect_equal(unname(vcov(reslm)), unname(vcov(res)), tolerance = 1e-2)
  
  
  ## Regression
  trait <- -3 + 2 * reg + simulateTipsCauchy(tree, 0, disp)
  dat <- data.frame(trait = trait, reg = reg)
  
  reslmstart <- cauphylm(trait ~ reg, data = dat, phy = tree, model = "lambda", starting.value = list(lambda = 1))
  reslm <- cauphylm(trait ~ reg, data = dat, phy = tree, model = "lambda")
  reslmnull <- cauphylm(trait ~ 1, data = dat, phy = tree, model = "lambda")
  reslmcau <- cauphylm(trait ~ reg, data = dat, phy = tree, model = "cauchy")

  expect_true(reslmnull$logLik <= reslm$logLik)
  expect_true(reslmcau$logLik <= reslmstart$logLik)
  # expect_true(reslmstart$logLik <= reslm$logLik)
  
})

test_that("cauphylm helper functions", {
  set.seed(1289)
  
  ## parameters
  ntips <- 20
  rateTest <- 0.1
  
  ## Tree
  tree <- rphylo(ntips, 0.1, 0)
  tree$edge.length <-  tree$edge.length / max(vcv(tree))
  tree_lambda <- phylolm::transf.branch.lengths(tree, model = "lambda", parameters = list(lambda = 0.6))$tree
  
  ## root.value
  mu <- 1
  disp <- 0.1
  
  ## Sim
  trait <- simulateTipsCauchy(tree_lambda, mu, disp)
  
  ## Reg
  reg <- ape::rTraitCont(tree, model = "BM", sigma = 0.1, root.value = 0)
  dat <- data.frame(tt = trait, rr = reg)
  
  ## Opt
  reslmdat <- cauphylm(tt ~ 1, data = dat, phy = tree, model = "lambda", hessian = TRUE)
  
  ## Functions
  reslmdat
  logLik(reslmdat)
  expect_equal(logLik(reslmdat)$logLik, 13.296526)
  expect_equal(AIC(reslmdat), -20.5930517)
  # expect_equal(extractAIC(reslmdat)[2], -20.5930517)
  # expect_equal(nobs(reslmdat), ntips)
  expect_equal(unique(predict(reslmdat)), unname(reslmdat$coefficients))
  expect_equal(unname(sqrt(diag(vcov(reslmdat)))), c(0.04988679, 0.06517500, 0.18057543), tolerance = 1e-2)
  expect_message(cc <- confint(reslmdat, level = 0.9), "Approximated asymptotic confidence interval using the Hessian")
  expect_true(cc[1, 1] <= mu && cc[1, 2] >= mu)
  expect_true(cc[2, 1] <= disp && cc[2, 2] >= disp)
  expect_true(cc[3, 1] <= 1.0 && cc[3, 2] >= 1.0)
  
  ## No Hessian
  reslmdatbis <- cauphylm(tt ~ 1, data = dat, phy = tree, model = "lambda", hessian = FALSE)
  expect_equal(vcov(compute_vcov(reslmdatbis)), vcov(reslmdat))
  
  ## Compare vcov with fitCauchy
  resfitcau <- fitCauchy(tree, trait, model = "lambda", method = "fixed.root", hessian = TRUE)
  expect_equal(vcov(resfitcau), vcov(reslmdat), ignore_attr = TRUE, tolerance = 1e-2)
  expect_equal(vcov(compute_vcov(fitCauchy(tree, trait, model = "lambda", method = "fixed.root", hessian = FALSE))),
               vcov(resfitcau))
  
  ## Functions fitCauchy
  resfitcau
  expect_equal(logLik(resfitcau)$logLik, logLik(reslmdat)$logLik)
  expect_equal(AIC(resfitcau), AIC(reslmdat))
  expect_equal(unname(sqrt(diag(vcov(resfitcau)))), unname(sqrt(diag(vcov(reslmdat)))), tolerance = 1e-2)
  expect_message(cc <- confint(resfitcau, level = 0.9), "Approximated asymptotic confidence interval using the Hessian")
  expect_true(cc[1, 1] <= mu && cc[1, 2] >= mu)
  expect_true(cc[2, 1] <= disp && cc[2, 2] >= disp)
  expect_true(cc[3, 1] <= 1.0 && cc[3, 2] >= 1.0)
  
})

# test_that("testFitBig", {
#   
#   ## Parameters
#   mu <- 0.8
#   disp <- 0.1
#   
#   ## tree with three tips
#   set.seed(1289)
#   n <- 500
#   tree <- rphylo(n, 0.1, 0)
#   tree_height <- max(diag(vcv(tree)))
#   
#   ## data
#   trait <- simulateTipsCauchy(tree, mu, disp)
#   
#   ## Opt
#   reslm <- cauphylm(trait ~ 1, phy = tree)
#   res <- fitCauchy(tree, trait, method = "fixed.root")
#   expect_equal(res$logLik, reslm$logLik)
#   
# })
test_that("testSimulationCauchyTree", {
  
  simCauchy <- function(phy, mu, disp) {
    cauchy_fun <- function(x, l, disp) {
      return(rcauchy(1, location = x, scale = disp * l))
    }
    return(rTraitCont(phy = phy,
                      model = cauchy_fun,
                      root.value = mu,
                      disp = disp))
  }
  
  set.seed(1289)
  
  phy <- rphylo(10, 0.1, 0)
  phy$edge.length <- phy$edge.length / max(vcv(phy))
  mu <- 1
  disp <- 0.1
  
  Nrep <- 100000
  trait1 <- sapply(1:Nrep, function(x) simulateTipsCauchy(phy, root.value = mu, disp = disp))
  expect_equal(unname(trait1[1, 1]), 0.93183503)
  expect_equal(unname(trait1[2, 1]), 1.29706603)
  expect_equal(mean(apply(trait1, 1, median)), mu, tolerance = 1e-4)
  expect_equal(mean(apply(trait1, 1, function(x) 0.5 * IQR(x))), disp, tolerance = 1e-4)
  expect_equal(ks.test(trait1[1, ], pcauchy, location = mu, scale = disp)$p.value, 0.307100108)
  
  Nrep <- 50000
  trait2 <- sapply(1:Nrep, function(x) simCauchy(phy, mu = mu, disp = disp))
  expect_equal(unname(trait2[1, 1]), 1.0967437)
  expect_equal(unname(trait2[2, 1]), 1.02307893)
  expect_equal(ks.test(trait1[5, 1:Nrep], trait2[5, ])$p.value, 0.43600472)
  
})

test_that("testrTraitCauchy", {
  
  phy <- rphylo(10, 0.1, 0)
  phy$edge.length <- phy$edge.length / max(vcv(phy))
  mu <- 1
  disp <- 0.1
  
  Nrep <- 1
  set.seed(1289)
  trait1 <- rTraitCauchy(n = Nrep, phy, model = "cauchy", parameters = list(root.value = mu, disp = disp))
  set.seed(1289)
  trait2 <- simulateTipsCauchy(phy, root.value = mu, disp = disp)
  expect_equal(trait1, trait2)
  
  Nrep <- 10
  set.seed(1289)
  trait1 <- rTraitCauchy(n = Nrep, phy, model = "cauchy", parameters = list(root.value = mu, disp = disp))
  set.seed(1289)
  trait2 <- sapply(1:Nrep, function(x) simulateTipsCauchy(phy, root.value = mu, disp = disp))
  expect_equal(trait1, trait2)
  
  Nrep <- 1
  set.seed(1289)
  traitlam <- rTraitCauchy(n = Nrep, phy, model = "lambda", parameters = list(root.value = mu, disp = disp, lambda = 0.5))
  traitkappa <- rTraitCauchy(n = Nrep, phy, model = "kappa", parameters = list(root.value = mu, disp = disp, kappa = 0.5))
  traitdelta <- rTraitCauchy(n = Nrep, phy, model = "delta", parameters = list(root.value = mu, disp = disp, delta = 0.5))
  
  set.seed(1289)
  trait1 <- rTraitCauchy(n = 1, phy, model = "cauchy", parameters = list(root.value = 0, disp = 1))
  set.seed(1289)
  trait2 <- rTraitCauchy(n = 1, phy)
  expect_equal(trait1, trait2)
  
  set.seed(1289)
  trait1 <- rTraitCauchy(n = Nrep, phy, model = "lambda", parameters = list(root.value = mu, disp = disp, lambda = 1))
  set.seed(1289)
  trait2 <- rTraitCauchy(n = 1, phy, parameters = list(root.value = mu, disp = disp, lambda = 1))
  expect_equal(trait1, trait2)
  
  set.seed(1289)
  trait1 <- rTraitCauchy(n = Nrep, phy, model = "kappa", parameters = list(root.value = mu, disp = disp, kappa = 1))
  expect_equal(trait1, trait2)
  set.seed(1289)
  trait1 <- rTraitCauchy(n = Nrep, phy, model = "delta", parameters = list(root.value = mu, disp = disp, delta = 1))
  expect_equal(trait1, trait2)
  
  
  
})

test_that("test rTraitCauchy errors", {
  
  phy <- rphylo(10, 0.1, 0)
  phy$edge.length <- phy$edge.length / max(vcv(phy))
  mu <- 1
  disp <- 0.1
  Nrep <- 10
  
  expect_error(rTraitCauchy(n = Nrep, phy, model = "cauchy", parameters = list(root.state = mu, disp = disp)),
               "Parameters root.state are unknown.")
  expect_error(rTraitCauchy(n = Nrep, phy, model = "cauchy", parameters = list(root.state = mu, disper = disp)),
               "Parameters root.state, disper are unknown.")
  
})
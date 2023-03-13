test_that("testMoments", {
  
  set.seed(128910)
  
  ## Parameters
  mu <- 0
  disp <- 1
  lambda <- 1
  
  ## tree with three tips
  tree <- read.tree(text = "((A:1,B:0.5):0.3, C:1);")
  
  ## Trait simulation
  trait <- rTraitCauchy(100000, tree, "cauchy", parameters = list(root.value = mu, disp = disp, lambda = lambda))
  
  ## tips
  expect_equal(mad(trait[1, ], constant = 1), vcv(tree)[1, 1], tolerance = 1e-2)
  expect_equal(mad(trait[2, ], constant = 1), vcv(tree)[2, 2], tolerance = 1e-2)
  expect_equal(mad(trait[3, ], constant = 1), vcv(tree)[3, 3], tolerance = 1e-2)
  
  ## Linear combinations
  expect_equal(mad(0.5 * trait[1, ] + 0.5 * trait[2, ], constant = 1),
               0.5 * vcv(tree)[1, 1] + 0.5 * vcv(tree)[2, 2],
               tolerance = 1e-2)
  expect_equal(mad(0.5 * trait[1, ] + 0.5 * trait[3, ], constant = 1),
               0.5 * vcv(tree)[1, 1] + 0.5 * vcv(tree)[3, 3],
               tolerance = 1e-2)
  expect_equal(mad(0.5 * trait[2, ] + 0.5 * trait[3, ], constant = 1),
               0.5 * vcv(tree)[2, 2] + 0.5 * vcv(tree)[3, 3],
               tolerance = 1e-2)
  
  a <- 0.4
  b <- -a
  expect_equal(mad(a * trait[1, ] + b * trait[2, ], constant = 1),
               (abs(a + b) - abs(a) - abs(b)) * vcv(tree)[1, 2] + abs(a) * vcv(tree)[1, 1] + abs(b) * vcv(tree)[2, 2],
               tolerance = 1e-2)
  
  ## With lambda
  lambda <- 0.7
  trait <- rTraitCauchy(100000, tree, "lambda", parameters = list(root.value = mu, disp = disp, lambda = lambda))
  expect_equal(median((abs(trait[1, ] - trait[2, ]) - vcv(tree)[1, 1] - vcv(tree)[2, 2]) / (- 2 * vcv(tree)[1, 2])),
               lambda, tolerance = 1e-1)
})

test_that("testEstimLambda", {
  
  set.seed(128910)
  
  ## Parameters
  mu <- 0
  disp <- 0.5
  lambda <- 0.7
  
  ## tree with three tips
  tree <- rphylo(50, 0.1, 0.02, fossils = TRUE)
  
  ## Trait simulation
  trait <- rTraitCauchy(1000, tree, "lambda", parameters = list(root.value = mu, disp = disp, lambda = lambda))
  
  ## disp
  estim_disp <- function(trait) {
    norm_trait <- (trait - median(trait)) / diag(ape::vcv(tree))
    return(mad(norm_trait, constant = 1))
  }
  
  all_disp <- apply(trait, 2, function(y) initDispersionParameter(tree, y, center = median(y), method.init.disp = "MAD"))
  expect_equal(mean(all_disp), disp, tolerance = 1e-1)
  
  ## lambda
  all_lambda_true_disp <- sapply(1:ncol(trait), function(i) initLambdaParameter(tree, trait[, i], disp))
  expect_equal(mean(all_lambda_true_disp), lambda, tolerance = 1e-1)
  
})
  
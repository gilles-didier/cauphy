test_that("reRootTree", {
  
  tree <- read.tree(text = "((A:1,B:0.5):0.3, C:1);")
  tip <- 1
  
  retree <- reroottip(tree, tip)
  
  ## Visual checks
  plot(tree)
  edgelabels(tree$edge.length)
  plot(retree, root.edge = T)
  edgelabels(retree$edge.length)
  
  ## check lengths
  revcv <- vcv(retree)
  expect_equal(0.5, revcv["B", "B"])
  expect_equal(1.3, revcv["C", "C"])
  expect_equal(retree$root.edge, 1.0)
  
  trimtree <- tree
  trimtree$edge.length[trimtree$edge[, 2] %in% which(trimtree$tip.label %in% c("A", "B"))] <- 0
  cc <- cophenetic(trimtree)
  expect_equal(cc["A", "C"], revcv["C", "C"])
  
  
  # Outlier tip
  tip <- 3
  retree <- reroottip(tree, tip)
  revcv <- vcv(retree)
  expect_equal(0.5, revcv["B", "B"])
  expect_equal(1.0, revcv["A", "A"])
  expect_equal(retree$root.edge, 1.3)
  
})

test_that("reRootTree bigger", {
  
  set.seed(1289)
  ntips <- 10
  tree <- rphylo(ntips, 0.1, 0)
  tree$edge.length <- round(tree$edge.length, 0)
  
  
  tip <- 3
  retree <- reroottip(tree, tip)
  
  ## Visual checks
  plot(tree)
  edgelabels(tree$edge.length)
  plot(retree, root.edge = T)
  edgelabels(retree$edge.length)
  
  ## check lengths
  revcv <- vcv(retree)
  cc <- cophenetic(tree) - retree$root.edge
  expect_equal(cc["t8", "t3"], revcv["t8", "t8"])
  expect_equal(cc["t5", "t3"], revcv["t5", "t5"])
  expect_equal(cc["t2", "t3"], revcv["t2", "t2"])
  expect_equal(retree$root.edge, tree$edge.length[tree$edge[, 2] %in% which(tree$tip.label %in% c("t3"))])
  
  ## Outlier tip
  tip <- 2
  treebis <- root(tree, outgroup = tip)
  treebis <- drop.tip(treebis, c("t1", "t7", "t10", "t6"))
  plot(treebis)
  edgelabels(treebis$edge.length)
  
  tip <- which(treebis$tip.label == "t2")
  retree <- reroottip(treebis, tip)
  expect_equal(retree$root.edge, 30)
  
})

test_that("read tree", {
  
  set.seed(1289)
  ntips <- 10
  tree <- rphylo(ntips, 0.1, 0)
  
  ## no root edge
  
  prtree <- printRTreeTest(tree)
  prtree <- substring(prtree, 1, nchar(prtree)-1)
  tree2 <- ape::read.tree(text = prtree)
  tree2$tip.label <- sub("'", "", tree2$tip.label)
  tree2$tip.label <- sub("'", "", tree2$tip.label)
  
  expect_true(all.equal(tree, tree2, tolerance = 1e-4))
  
  ## with root edge
  tree$root.edge <- 10
  
  prtree <- printRTreeTest(tree)
  prtree <- substring(prtree, 1, nchar(prtree)-1)
  tree2 <- ape::read.tree(text = prtree)
  tree2$tip.label <- sub("'", "", tree2$tip.label)
  tree2$tip.label <- sub("'", "", tree2$tip.label)
  
  expect_true(all.equal(tree, tree2, tolerance = 1e-4))
  
})


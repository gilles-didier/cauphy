# test_that("testExport", {
# 
#   ## Parameters
#   mu <- 0.0
#   disp <- 0.6
# 
#   ## tree with three tips
#   set.seed(1289)
#   n <- 20
#   tree <- rphylo(n, 0.1, 0.02, fossils = TRUE)
#   ntaxa <- length(tree$tip.label)
# 
#   ## data
#   trait_lat <- 34 + rTraitCont(tree, model = "BM", sigma = 0.1)
#   trait_lat[c(2, 18, 21)] <- trait_lat[c(2, 18, 21)] + 2
#   trait_lat[c(1, 13, 14, 20)] <- trait_lat[c(1, 13, 14, 20)] - 4
#   trait_long <- -80 + rTraitCont(tree, model = "BM", sigma = 0.1)
# 
#   ## ASR - reml
#   fit_lat <- fitCauchy(tree, trait_lat, method = "reml")
#   fit_long <- fitCauchy(tree, trait_long, method = "reml")
# 
#   anc_lat <- ancestral(fit_lat, n_values = 1000)
#   anc_long <- ancestral(fit_long, n_values = 1000)
# 
#   plot(tree); nodelabels(); tiplabels()
#   plot(anc_lat, node = 38)
#   plot(anc_lat, node = 36)
#   plot(anc_lat, node = 30)
# 
#   tree_with_anc <- ancestral_to_treedata(tree, fit_lat, fit_long, anc_lat, anc_long, level = 0.80)
#   tree_with_anc_rec <- ancestral_to_treedata_rec(tree, fit_lat, fit_long, anc_lat, anc_long, level = 0.80)
# 
#   # treeio::write.beast(tree_with_anc, file = "inst/test_tree.tree")
# 
#   # test <- read.nexus(file = "inst/test_tree.tree")
#   # plot(test)
# 
# })
# 
# test_that("testExport", {
# 
#   ## data
#   n <- 20
#   tree <- read.tree(text = "((AF404755_Bu_2000:0.3618982163,(DQ164202_Hs_2002:1.76071759,DQ080071_Ec_2002:1.516730968):0.7023872207):0.3711544571,(DQ431701WG214_Hs_2004.53:3.858441141,(((WG142_Hs_2006.64:5.652539455,(DQ431705WG233_Hs_2004.52:1.827524824,(WG024_Hs_2003.53:0.3930328678,DQ431694WG022_Hs_2003.62:0.4830328678):0.4444919562):1.705014631):0.1357036686,(AY712946_Cc_2003:2.217929827,(DQ080059_Pn_2003:0.9713525334,WG080_Hs_2004.56:1.652915653):1.578836395):0.4764910763):0.01470759248,(AY712948_Cq_2003:2.349029222,((DQ666450_Hs_2005:4.370836638,DQ164198_Hs_2002:1.424190227):0.1074749768,(DQ431699WG124_Hs_2003.73:2.725470244,(DQ164203_Ph_2003:2.254305112,(DQ666452_Hs_2005:3.947260488,(DQ080065_Qq_2003:1.599417224,(DQ666448_Hs_2004:1.706006917,DQ431704WG219_Hs_2004.57:1.702127249):1.006918897):0.3282027029):0.2679864214):0.2602349738):0.00590988911):0.01571614976):0.1458544329):0.1654904246):0.9074364635);")
# 
#   lat <- c(43.46000, 29.94453, 29.85942, 38.58000, 30.40262, 28.70000, 31.19032, 40.19000, 39.24414, 34.65000, 27.97657, 39.14801, 39.13490, 44.07000, 33.98030, 31.35336, 44.22000, 33.81000, 34.22000, 41.26000)
#   long <- c(-76.24000, -95.45249, -95.14223, -121.49000, -115.07045, -82.00000, -99.86344, -82.67000, -105.75979, -102.75000, -82.58716, -108.31023, -108.65065, -103.21000, -112.05157, -99.91597, -100.25000, -117.83000, -118.44000, -96.12000)
# 
#   names(lat) <- names(long) <- tree$tip.label
# 
#   ## ASR - reml
#   fit_lat <- fitCauchy(tree, lat, method = "reml")
#   fit_long <- fitCauchy(tree, long, method = "reml")
# 
#   anc_lat <- ancestral(fit_lat)
#   anc_long <- ancestral(fit_long)
# 
#   plot(anc_lat, node = 22)
# 
#   tree_with_anc <- ancestral_to_treedata(tree, fit_lat, fit_long, anc_lat, anc_long, level = 0.8)
#   tree_with_anc_rec <- ancestral_to_treedata_rec(tree, fit_lat, fit_long, anc_lat, anc_long, level = 0.80)
# 
#   # write.evolaps(tree_with_anc, file = "inst/test_tree.tree")
# })
# 

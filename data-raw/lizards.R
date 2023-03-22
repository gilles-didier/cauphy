## Data:
## Mahler, D. Luke; Ingram, Travis; Revell, Liam J.; Losos, Jonathan B. (2013),
## Data from: Exceptional convergence on the macroevolutionary landscape in island lizard radiations,
## Dryad, Dataset, https://doi.org/10.5061/dryad.9g182

## Phylogenetic Tree
phy <- read.tree(file = file.path("data-raw", "Mahler_et_al_2013_Data", "GA_Anolis_MCC.tre"))
## Data
dat <- read.csv(file = file.path("data-raw", "Mahler_et_al_2013_Data", "GA_Anolis_traits.csv"))
## ecomorphs
ecomorphs <- read.csv(file.path("data-raw", "Mahler_et_al_2013_Data", "GA_Anolis_trad_ecomorph_class.csv"),
                      header = FALSE)
ecomorphs <- ecomorphs[match(phy$tip.label, ecomorphs[, 1]), ]
ecomorphs <- factor(ecomorphs[, 2])
levels(ecomorphs) <- c("Trunk-Ground",
                       "Trunk-Crown",
                       "Crown-Giant",
                       "Twig",
                       "Grass-Bush",
                       "Trunk",
                       "Unique")
ecomorphs <- factor(ecomorphs, level = sort(levels(ecomorphs)))

## Keep only svl
svl <- dat[, "AVG.SVL"]
names(svl) <- dat$species
svl <- suppressWarnings(cauphy:::checkTraitTree(svl, phy))

## Dataset
lizards <- list(phy = phy,
                svl = svl,
                ecomorph = ecomorphs)

usethis::use_data(lizards, overwrite = TRUE)

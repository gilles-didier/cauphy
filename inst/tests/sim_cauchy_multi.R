##
#' @title Transformation
#'
#' @description
#' Transform iid bivariate vectors to transformed vectors in the space
#' defined by axes theta1 and theta2, with scales sigma1 and sigma2.
#'
#' @param Z matrix of size d times nrep, with iid draws.
#' @param theta vector of angles
#' @param sigma vector of associated scales
#' 
#' @return
#' The transformed vectors
#' 
#'
transform_multi <- function(Z, theta, sigma = rep(1, length(theta))) {
  M <- sapply(theta, function(x) c(cos(x), sin(x)))
  Y <- M %*% diag(sigma) %*% Z
  return(Y)
}

library(ggplot2)
nrep <- 1e6
d <- 2

## Gaussians
set.seed(1289)
Z <- matrix(rnorm(d * nrep), nrow = d)

## Transform
YG <- transform_multi(Z, c(pi / 3, pi / 3 + pi / 2), c(0.5, 1))

# ggplot(data.frame(x = YG[1, ], y = YG[2, ]), aes(x = x, y = y)) + geom_hex(bins = 100) +
#   scale_fill_continuous(type = "viridis") + theme_bw() +
#   coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10))

# ggplot(data.frame(x = Y[1, ], y = Y[2, ]), aes(x = x, y = y)) +
  # geom_density_2d_filled(contour_var = "ndensity")

## Cauchy
set.seed(1289)
th <- 10
Z <- matrix(rcauchy(d * nrep), nrow = d)

## Transform orthogonal
angle <- c(pi / 3, pi / 3 + pi / 2)
disp <- c(0.05, 0.1)
YC1 <- transform_multi(Z, angle, disp)
get_cut <- function(Y, th) Y[1, ] <= th & Y[1, ] >= -th & Y[2, ] <= th & Y[2, ] >= -th

cutC1 <- get_cut(YC1, th)

ggplot(data.frame(x = YC1[1, cutC1], y = YC1[2, cutC1]), aes(x = x, y = y)) + geom_hex(bins = 200) +
  scale_fill_continuous(type = "viridis") + theme_bw() +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10))

x <- seq(-10, 10, 0.1)
vals <- expand.grid(x, x)
dens <- logDensityCauchyBi(vals, center = c(0, 0), angle = angle, disp = disp)

ggplot(data.frame(x = vals[, 1], y = vals[, 2], z = dens),
       aes(x = x, y = y, z = z)) + 
  geom_contour() +
  scale_fill_continuous(type = "viridis") + theme_bw() +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10))

## Transform non orthogonal
YC2 <- transform_multi(Z, c(pi / 3, pi / 3 + pi / 4), c(0.05, 0.1))

cutC2 <- get_cut(YC2, th)

# ggplot(data.frame(x = YC2[1, cutC2], y = YC2[2, cutC2]), aes(x = x, y = y)) + geom_hex(bins = 200) + 
#   scale_fill_continuous(type = "viridis") + theme_bw() +
#   coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10))

# ggplot(data.frame(x = Y[1, cut], y = Y[2, cut]), aes(x = x, y = y)) +
#   geom_density_2d_filled(contour_var = "density") +
#   coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10))

## Plot
data_plot <- rbind(
  data.frame(x = YG[1, ], y = YG[2, ], alpha = 2, theta1 = 60, theta2 = 150),
  data.frame(x = YC1[1, cutC1], y = YC1[2, cutC1], alpha = 1, theta1 = 60, theta2 = 150),
  data.frame(x = YC2[1, cutC2], y = YC2[2, cutC2], alpha = 1, theta1 = 60, theta2 = 105)
)

ggplot(data_plot, aes(x = x, y = y)) + geom_hex(bins = 200) + 
  scale_fill_continuous(type = "viridis") + theme_bw() +
  coord_fixed(xlim = c(-10, 10), ylim = c(-10, 10)) +
  facet_wrap(~alpha + theta1 + theta2, labeller = "label_both") +
  geom_abline(aes(slope = sin(theta1 * pi / 180) / cos(theta1 * pi / 180), intercept = 0)) +
  geom_abline(aes(slope = sin(theta2 * pi / 180) / cos(theta2 * pi / 180), intercept = 0))

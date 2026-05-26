library(burnout)
library(dplyr)
library(ggplot2); theme_set(theme_bw())

calc_bp <- function(R0_range = c(1.1, 5),
                    eps_range = c(0.02, 0.25)) {
  cc <- emdbook::curve3d(
    burnout_prob(R0, epsilon, N= 1e6),
    varnames = c("R0", "epsilon"),
    from = c(R0_range[1], eps_range[1]),
    to = c(R0_range[2], eps_range[2]),
    n = c(51,51))
  ## convert to long format
  dd <- tibble(
    R0 = rep(cc$x, length(cc$y)),
    epsilon = rep(cc$y, each = length(cc$x)),
    bp = c(cc$z))
  dd
}

dd1 <- calc_bp()
dd2 <- calc_bp(eps_range = c(0.001, 0.02))
ggplot(dd2, aes(epsilon, R0)) +
  geom_raster(aes(fill = 1-bp)) +
  scale_fill_viridis_c() ## trans = "log10"

plot_P1(epsilon=0.125, N=10^4)

xvec <- seq(1.1, 32, length = 51
curve(1- fizzle_prob(x), from = 1.1, to = 32)
curve(1-1/x, col = 2, add = TRUE)

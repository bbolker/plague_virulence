library(dplyr)
library(tidyr)
library(odin)
library(dde) ## odin insists on this
library(ggplot2); theme_set(theme_bw())
library(parallel)

source(here::here("metapop/odin", "discrete_run.R"))

nsim   <- 20
ncores <- max(1L, detectCores() - 1L)

first_zero <- function(x) which(x == 0)[1]

## parameter note
## plague generation time ~ 10 - 20 days
## rat $r$ ~ 3 / year?  ~ 0.125?
sumfun <- function(beta1, beta2, K = 1e6, I_init = c(10, 10), cl = NULL) {
  runs <- discrete_run(beta_vec = c(beta1, beta2), K = K, r = 0.125,
                       n_patch = 1, nt = 1000, alpha = 0, I_init = I_init,
                       nsim = nsim, cl = cl, platform = "odin")
  if (!is.list(runs)) runs <- list(runs)
  ext <- sapply(runs, function(traj) {
    c(I1 = first_zero(traj$value[traj$state == "I1"]),
      I2 = first_zero(traj$value[traj$state == "I2"]))
  })
  rowMeans(ext, na.rm = TRUE)
}

R0vec <- seq(1, 5, by = 0.025)
dd <- expand.grid(R01 = R0vec, R02 = R0vec)

set.seed(101)
cl <- makeCluster(ncores)
clusterSetRNGStream(cl)
system.time(
  res <- apply(dd, 1, \(x) sumfun(x[1], x[2], cl = cl))
)
stopCluster(cl)

res_df <- t(res) |> as.data.frame() |> bind_cols(dd)
dd2    <- res_df |> pivot_longer(cols = c("I1", "I2"))
dd2_i1 <- dd2 |> dplyr::filter(name == "I1")

ggplot(dd2_i1, aes(R01, R02, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10", na.value = adjustcolor("lightyellow", alpha.f = 0.5)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  labs(x = expression(R[0]*' of strain 1'),
       y = expression(R[0]*' of strain 2'),
       title = expression('extinction time when '*list(K==10^6, r==0.125))) +
  annotate(x = 2.5, y = 1.5, size = 10, label = "time > 1000", geom = "label")

ggsave(width = 6.5, height = 6, "twostrain_onepatch.png")

library(dplyr)
library(tidyr)
library(odin)
library(dde) ## odin insists on this
library(ggplot2); theme_set(theme_bw())

nsim <- 20

## furrr/future/parallel? don't play nicely with odin
## library(furrr)
## plan(multisession, workers = 5)

odin_fn <- here::here("metapop/odin", "odin_twostrain1.R")
twostrain_generator <- odin::odin(odin_fn)

first_zero <- function(x) {
  which(x==0)[1]
}

## parameter note
## plague generation time ~ 10 - 20 days
## rat $r$ ~ 3 / year?  ~ 0.125?
sumfun <- function(beta1 = 3, beta2 = 3, K = 1e6, I_init = c(10, 10), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  odin_fn <- twostrain_generator$new(beta = c(beta1, beta2),
                               r = 0.125,
                               K = K,
                               I_ini = I_init,
                               S_ini = K - sum(I_init))
  res <- odin_fn$run(0:1000)
  res2 <- res[, !colnames(res) %in% c("step", "S"), drop = FALSE]
  apply(res2, 2, first_zero)
}

R0vec <- seq(1, 5, by = 0.025)
dd <- expand.grid(R01 = R0vec, R02 = R0vec)


set.seed(101)
## odin doesn't work well with future
system.time(
  res <- purrr::map(seq.int(nsim),
                    \(y) apply(dd, 1, \(x) sumfun(x[1], x[2])) |> t() |> bind_cols(dd)) |>
    bind_rows(.id = "run")
)

res2 <- res |> select(-run) |> summarise(across(starts_with("I"), mean), .by = starts_with("R"))
                         
dd2 <- cbind(res2) |>
  pivot_longer(cols = !starts_with("R0"))

dd2_i1 <- (dd2 |> dplyr::filter(name == "I[1]"))

ggplot(dd2_i1, aes(R01, R02, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10", na.value = adjustcolor("lightyellow", alpha.f = 0.5))  +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_abline( intercept = 0, slope = 1, colour = "red") +
  labs(x = expression(R[0]*' of strain 1'),
       y = expression(R[0]*' of strain 2'),
       title = expression('extinction time when '*list(K==10^6,r==0.125))) +
  annotate(x = 2.5, y = 1.5, size = 10, label = "time > 1000", geom = "label")

ggsave(width = 6.5, height = 6, "twostrain_onepatch.png")
       
  
## ext <- which(rowSums(is.na(res))==0 & res[,"I[1]"]+res[,"I[2]"] == 0)[1]
## res <- res[1:(ext-1),]

## pdf("odin_twostrain_run.Rout.pdf")
## matplot(res[,1], res[,-1], type = "l", log = "y", lwd = 2.5, col = c(1,2,3,4), lty = 1)
## legend("right", legend = colnames(res)[-1], lty = 1, col = 1:4)
## dev.off()

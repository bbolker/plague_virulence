library(dplyr)
library(ggplot2)
library(patchwork)
library(optparse)
zmargin <- theme(panel.spacing=grid::unit(0,"lines"))

h <- function(x) here::here("metapop/odin", x)

fits <- c("logistic_continuous", "logistic_reedfrost",
           "linear_continuous",   "linear_reedfrost")

outfn <- h("euler_onepatch_onestrain_comb")
fn_vec <- sprintf("%s/euler_onepatch_onestrain_extinct_%s.rds",
                  h("sharcnet/outputs"), fits)

dat <- fn_vec |>
  setNames(fits) |>
  lapply(readRDS) |>
  dplyr::bind_rows(.id = "model") |>
  mutate(log10K = log10(K)) |>
  tidyr::separate(model, sep = "_",
                  into = c("demography", "time-step")) |>
  mutate(`time-step` = if_else(`time-step` == "reedfrost", "discrete", "continuous"))
                          
theme_set(theme_bw(base_size=16))

## FIXME: convert to long and facet rather than using patchwork
## (will free scales on the fill guide work? maybe not worth it?)
## Hmm, this doesn't seem easy: https://stackoverflow.com/q/45109293/190277

gg0 <- ggplot(dat, aes(log10K, R0, fill = ext_prob.I1)) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10", na.value = "grey80",
                       name  = "extinction probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = 3:6,
                     labels = scales::label_math(10^.x)) + 
  labs(y = expression(R[0]), x = "K") +
  facet_grid(`time-step` ~ demography, labeller = label_both) +
  zmargin

dd_ann <- tibble::tibble(demography = "linear", `time-step` = "continuous",
                         log10K = 5.5, R0 = 3.5,
                         ext_prob.I1 = NA,
                         label = expression("extprob" < 10^{-3}))

gg0 + geom_label(data = dd_ann, aes(label = label),
                 fill = "white")

ggsave(paste0(outfn, ".png"), width = 12, height = 8)
ggsave(paste0(outfn, ".pdf"), width = 12, height = 8)

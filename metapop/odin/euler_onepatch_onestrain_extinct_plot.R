library(dplyr)
library(ggplot2)
library(patchwork)
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-c", "--combo"), type = "character",
              default = "logistic_reedfrost",
              help = "combination to plot: logistic_continuous, logistic_reedfrost, linear_continuous, linear_reedfrost [default: %default]")
)))

valid <- c("logistic_continuous", "logistic_reedfrost",
           "linear_continuous",   "linear_reedfrost")
if (!opt$combo %in% valid)
  stop("--combo must be one of: ", paste(valid, collapse = ", "))

base_fn <- paste0("euler_onepatch_onestrain_extinct_", opt$combo)
fn <- here::here("metapop/odin", "sharcnet/outputs",
                 paste0(base_fn, ".rds"))
dat <- readRDS(fn) |>
  mutate(log10K = log10(K))

theme_set(theme_bw())

## FIXME: convert to long and facet rather than using patchwork
## (will free scales on the fill guide work? maybe not worth it?)
## Hmm, this doesn't seem easy: https://stackoverflow.com/q/45109293/190277

make_raster <- function(fill_var, fill_label) {
  ggplot(dat, aes(R0, log10K, fill = .data[[fill_var]])) +
    geom_raster() +
    scale_fill_viridis_c(trans = "log10", na.value = "grey80",
                         name = fill_label) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = 3:6,
                       labels = scales::label_math(10^.x)) +
    labs(x = expression(R[0]), y = "K")
}

p1 <- make_raster("ext_prob.I1",      "extinction\nprobability")
p2 <- make_raster("mean_ext_time.I1", "mean extinction\ntime (steps)")

print(p1 + p2)
ggsave(paste0(base_fn, ".png"), width = 10, height = 5)
ggsave(paste0(base_fn, ".pdf"), width = 10, height = 5)

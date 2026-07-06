library(plagueMetapop)
library(dplyr)
library(arrow)
library(ggplot2); theme_set(theme_minimal())
## download big file:
##  scp nibi.sharcnet.ca:~/project/bolker/plague_virulence/metapop/odin/sharcnet/outputs/euler_twostrain_examples_task_000003.rds .

fn <- "metapop/odin/outputs/euler_twostrain_examples_task_000003.rds"
system.time(x <- readRDS(here::here(fn)))

dt <- 0.1
strain2_delay <- 50

xa <- arrow_table(x)

## delay
xamin <- xa |> select(step:run) |>
    filter(step > strain2_delay, abs(step %% 1)<1e-6) |>
    mutate(across(run, as.integer))
    
mean_nonzero_arrow <- function(x) {
  # Compose Arrow compute ops: filter zeros then compute mean
  # Replace 0s with NA, then Arrow's mean() ignores NAs by default
  x_masked <- if_else(x == 0, NA_real_, x)
  mean(x_masked, na.rm = TRUE)
}

sd_nonzero_arrow <- function(x) {
  # Compose Arrow compute ops: filter zeros then compute mean
  # Replace 0s with NA, then Arrow's mean() ignores NAs by default
  x_masked <- if_else(x == 0, NA_real_, x)
  sd(x_masked, na.rm = TRUE)
}

xsum <- xamin |> summarise(patches = sum(value>0),
                           mean = mean_nonzero_arrow(value),
                           sd = sd_nonzero_arrow(value),
                           .by = c(state, step, run))

m <- collect(xsum)

mL <- m |> tidyr::pivot_longer(c(patches, mean, sd))

gg0 <- ggplot(mL, aes(step, value)) +
  geom_line(aes(colour = state), alpha = 0.25) +
  facet_wrap(~name, scale = "free")
ggsave(gg0, height = 8, width = 16, file = "test_twostrain_example.png")
  

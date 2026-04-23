dd <- readRDS("outputs_2strain/finalsize.rds")
library(gsl)
library(tidyverse)
library(ggrastr)
theme_set(theme_bw())
finalsize <- function(R0) {
  1+1/R0*lambert_W0(-R0*exp(-R0))
}

## I’m just using (R01+R02)/2 to calculate the final size in of the patch and partitioning the final death toll proportionally based on their initial number.

dd2 <- (dd
  |> as_tibble()
  |>   mutate(
         finalsize  = 1 -S,
         finalsize_yy = finalsize((R01+R02)/2),
         I1tot_yy = finalsize_yy*I10/(I10+I20),
         I2tot_yy = finalsize_yy*I20/(I10+I20))
)

simpars <- names(dd)[1:4]
respars <- grepv("^(finalsize|I[12]tot)(_yy)?$", names(dd2))

dd3 <- dd2 |>
  select(any_of(simpars), any_of(respars)) |>
  pivot_longer(any_of(respars)) |>
  mutate(approx = ifelse(grepl("_yy$", name), "approx", "sim"),
         across(name, ~ stringr::str_remove(name, "_yy$"))
         ) |>
  pivot_wider(names_from = approx, values_from = value)

gg0 <- ggplot(dd3, aes(sim, approx)) +
  facet_wrap(~name) +
  geom_abline(intercept=0, slope= 1, colour = "red") +
  labs(x = "Exact (ODE) value", y = "Approximation",
       title = "Comparison of approximate and exact final sizes")

gg1 <- gg0 + rasterise(geom_point(aes(colour = R01))) 

print(gg1)

## trying to figure out other ways to plot ...

## https://teunbrand.github.io/ggchromatic/
## remotes::install_github("teunbrand/ggchromatic")
## gg2 <- gg0 + geom_point(aes(colour = hcl_spec(R01, I10, l = 0.5)))

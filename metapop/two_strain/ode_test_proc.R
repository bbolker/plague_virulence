dd <- readRDS("outputs_2strain/finalsize.rds")
library(gsl)
library(tidyverse)
library(ggrastr)
library(biscale)
library(cowplot)
theme_set(theme_bw())
finalsize <- function(R0) {
  1+1/R0*lambert_W0(-R0*exp(-R0))
}

## Yuyang:
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

bi_data <- bi_class(dd3, x = R01, y = R02, style = "quantile", dim = 3)

gg2 <- gg0 + bi_data + rasterise(geom_point(aes(colour = bi_class), show.legend = FALSE)) +
  bi_scale_color(pal = "GrPink", dim = 3)

gg_legend <- bi_legend(pal = "GrPink",
                    dim = 3,
                    xlab = "R01",
                    ylab = "R02",
                    size = 8)

ggdraw() +
   draw_plot(gg2, 0, 0, 1, 1) +
   draw_plot(gg_legend, x = 0.8, y=0.1, 0.2, 0.2)

dd4 <- mutate(dd2,
              logit_ratio = qlogis(I1tot/(I1tot+I2tot)),
              logit_finalsize = qlogis(finalsize))

vars <- c("R01", "R02","log(I10)","log(I20)")
mk_form <- function(vars, rawpoly = FALSE, resp = "logit_ratio") {
  reformulate(sprintf("poly(%s, degree = 2, raw = %s)",
                      paste(vars, collapse = ","),
                      as.character(rawpoly)),
              response = resp)
}

fit1 <- lm(mk_form(vars), data = dd4)
fit1R <- lm(mk_form(vars, rawpoly = TRUE), data = dd4)

fit2 <- lm(mk_form(vars, resp = "logit_finalsize"), data = dd4)
fit2R <- lm(mk_form(vars, rawpoly = TRUE, resp = "logit_finalsize"), data = dd4)

print(summary(fit1)$adj.r.squared)
print(summary(fit2)$adj.r.squared)

stopifnot(all.equal(summary(fit1)$adj.r.squared, summary(fit1R)$adj.r.squared))
cc1R <- coef(fit1R)
names(cc1R) <- stringr::str_extract(names(cc1R), "(.Intercept.|([0-2]\\.){3}[0-2])")

## augment names with names of vars -- split, add vars, collapse ... ?

plot(sort(cc1R))
abline(h=c(-0.2, 0.2),lty = 2)
## pick out large values?
cc1R[abs(cc1R)>0.2]

## phenomenological but could work ... ??

## trying to figure out other ways to plot ...
## https://teunbrand.github.io/ggchromatic/
## remotes::install_github("teunbrand/ggchromatic")
## gg2 <- gg0 + geom_point(aes(colour = hcl_spec(R01, I10, l = 0.5)))

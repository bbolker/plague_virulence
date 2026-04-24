dd <- readRDS("outputs_2strain/finalsize.rds")
library(gsl)
library(tidyverse)
library(ggrastr)
library(biscale)
library(cowplot)
theme_set(theme_bw())
zmargin <- theme(panel.spacing=grid::unit(0,"lines"))

finalsize <- function(R0) {
  1+1/R0*lambert_W0(-R0*exp(-R0))
}

## Yuyang:
## I’m just using (R01+R02)/2 to calculate the final size in of the patch and partitioning the final death toll proportionally based on their initial number.

dd2 <- (dd
  |> as_tibble()
  ## don't bother with I2tot (or yy analogue), determined by
  ## finalsize + I1tot
  |> select(-I2tot)
  |>   mutate(
         finalsize  = 1 -S,
         finalsize_yy = finalsize((R01+R02)/2),
         I1tot_yy = finalsize_yy*I10/(I10+I20)
  )
)

## check approximations
dd4 <- mutate(dd2,
              logit_ratio = qlogis(I1tot/finalsize),
              logit_finalsize = qlogis(finalsize))

vars <- c("R01", "R02","log(I10)","log(I20)")

## construct 2d-order poly formula based on vars
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

## raw and orthogonal polynomials give same result
stopifnot(all.equal(summary(fit1)$adj.r.squared, summary(fit1R)$adj.r.squared))
cc1R <- coef(fit1R)
names(cc1R) <- stringr::str_extract(names(cc1R), "(.Intercept.|([0-2]\\.){3}[0-2])")

## augment names with names of vars -- split, add vars, collapse ... ?

plot(sort(cc1R))
abline(h=c(-0.2, 0.2),lty = 2)
## pick out large values?
cc1R[abs(cc1R)>0.2]


## add poly-regression predictions
dd3 <- dd2 |>
  mutate(finalsize_poly = plogis(predict(fit2, newdata = dd2)),
         I1tot_poly = finalsize_poly*plogis(predict(fit1, newdata = dd2))
         ) |>
  rename(I1tot_sim = "I1tot", finalsize_sim = "finalsize")


## visualize YY approximation
simpars <- names(dd)[1:4]
respars <- grepv("^(finalsize|I[12]tot)(_(yy|poly|sim))?$", names(dd3))

dd4 <- dd3 |>
  select(any_of(simpars), any_of(respars)) |>
  pivot_longer(any_of(respars)) |>
  mutate(approx = stringr::str_extract(name,"(?<=_)(yy|poly|sim)$"),
         across(name, ~ stringr::str_remove(name, "_.*$"))
         ) |>
  pivot_wider(names_from = approx, values_from = value)

gg0 <- ggplot(dd4, aes(sim, yy)) +
  facet_wrap(~name) +
  zmargin + 
  geom_abline(intercept=0, slope= 1, colour = "red") +
  labs(x = "Exact (ODE) value", y = "Approximation",
       title = "Comparison of approximate (YY) and exact final sizes")

gg1 <- gg0 + rasterise(geom_point(aes(colour = R01))) 

## no longer care about this, the one with a bivariate colour
## scheme is nicer
## print(gg1)

bi_data <- bi_class(dd4, x = R01, y = R02, style = "quantile", dim = 3)

gg2 <- gg0 + bi_data + rasterise(geom_point(aes(colour = bi_class), show.legend = FALSE)) +
  bi_scale_color(pal = "GrPink", dim = 3)

gg_legend <- bi_legend(pal = "GrPink",
                    dim = 3,
                    xlab = "R01",
                    ylab = "R02",
                    size = 8)

draw_fun <- function(gg) {
  ggdraw() +
          draw_plot(gg, 0, 0, 1, 1) +
          ## position hand-tweaked (top left of right-hand facet)
          draw_plot(gg_legend, x = 0.55, y=0.7, 0.2, 0.2)
}

gg3 <- draw_fun(gg2)
print(gg3)

gg2B <- gg2 +
  aes(y = poly) +
  labs(title = "Comparison of approximate (poly) and exact final sizes")

gg4 <- draw_fun(gg2B)

print(gg4)
## phenomenological but could work ... ??

#' generate predicted outcomes from existing polynomial fits
#' @examples
#' pred_outcomes_poly(1.2, 1.3, 0.1, 0.08)
pred_outcomes_poly <- function(R01, R02, I10, I20) {
  nd <- data.frame(R01, R02, I10, I20)
  pfun <- function(f) unname(plogis(predict(f, newdata = nd)))
  finalsize <- pfun(fit2)
  I1tot <- finalsize*pfun(fit1)
  I2tot <- finalsize-I1tot
  unlist(tibble::lst(finalsize, I1tot, I2tot))
}

saveRDS(pred_outcomes_poly, file = "polyfit.rds")


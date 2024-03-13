library(deSolve)
library(ggplot2); theme_set(theme_bw())
library(tidyverse)
library(patchwork)

gradfun <- function(t, y, parms) {
    comb <- as.list(c(y, parms))
    attach(comb)
    on.exit(detach(comb))
    beta1 <- R01*gamma1; beta2 <- R02*gamma2
    inc1 <- beta1*S*I1
    inc2 <- beta2*S*I2
    g <- c(
        S = -(inc1 + inc2),
        I1 = inc1 - gamma1*I1,
        I2 = inc2 - gamma2*I2,
        R = gamma1*I1 + gamma2*I2)
    list(g)
}

p0 <- c(R01=2, R02=2, gamma1 = 2, gamma2 = 1)
y0 <- c(S = 1, I1 = 0.01, I2 = 0.01, R= 0)
tvec <- seq(0, 20, by = 0.1)
res1 <- ode(y0, times = tvec, func = gradfun, parms = p0)
matplot(tvec, res1[,-1], type = "l", log = "y", ylim = c(1e-4, 1), type = "l")

res1T <- (res1
    |> as.data.frame()
    |> as_tibble()
    |> mutate(ratio = I1/(I1+I2))
    |> pivot_longer(-time, names_to = "var")
    |> mutate(across(var, \(x) factor(x,
                                      levels=c("S", "I1", "I2", "R", "ratio"))
                     ))
)

gg0 <- filter(res1T, var!="ratio") |> droplevels() |>
    ggplot(aes(time, value, colour = var)) +
    geom_line() +
    scale_y_log10(limits = c(1e-5, NA))

gg1 <- filter(res1T, var=="ratio") |> droplevels() |>
    ggplot(aes(time, value)) +
           geom_line() +
    labs(y="ratio I1/(I1+I2)") +
    geom_hline(yintercept =  0.5, lty = 2)

gg0 + gg1
ggsave("fastslow.png")

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
        R = gamma1*I1 + gamma2*I2,
        cumfoi1 = beta1*I1,
        cumfoi2 = beta2*I2)
    list(g)
}

p0 <- c(R01=2, R02=2, gamma1 = 2, gamma2 = 1.5)
varstr <- with(as.list(p0),
               sprintf("R01=%1.1f, R02=%1.1f, g1=%1.1f, g2=%1.1f",
                       R01, R02, gamma1, gamma2))
y0 <- c(S = 1, I1 = 0.01, I2 = 0.01, R= 0, cumfoi1 = 0, cumfoi2 = 0)
tvec <- seq(0, 20, by = 0.1)
res1 <- ode(y0, times = tvec, func = gradfun, parms = p0)
matplot(tvec, res1[,-1], log = "y", ylim = c(1e-4, 10), type = "l")

levs <- c("S", "I1", "I2", "R", "ratio", "cumfoi1", "cumfoi2")
res1T <- (res1
    |> as.data.frame()
    |> as_tibble()
    |> mutate(ratio = I1/(I1+I2))
    |> pivot_longer(-time, names_to = "var")
    |> mutate(across(var, \(x) factor(x,levels=levs)))
)

gg0 <- ggplot(res1T, aes(time, value)) + geom_line()

dat1 <- filter(res1T, var %in% levs[1:4]) |> droplevels()
gg1 <- (gg0 %+% dat1) +
    aes(colour = var) +
    scale_y_log10(limits = c(1e-5, NA)) +
    labs(title = "state variables")

dat2 <- filter(res1T, var=="ratio") |> droplevels()
gg2 <- (gg0 %+% dat2) +
    labs(y="ratio I1/(I1+I2)") +
    geom_hline(yintercept =  0.5, lty = 2) +
    labs(title = "prevalence ratio")

dat3 <- filter(res1T, startsWith(as.character(var), "cum")) |> droplevels()
gg3 <- (gg0 %+% dat3) +
    aes(colour = var) +
    scale_y_log10() +
    labs(title = "cumulative FOI")

ggcomb <- gg1 + gg2 + gg3 + plot_annotation(title = varstr)
ggsave(ggcomb, file = "fastslow.pdf", width = 16)

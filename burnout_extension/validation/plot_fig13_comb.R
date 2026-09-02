library(here)
source(here('graphics_utils.R'))
setwd(here('burnout_extension'))
source('validation/R/theory.R')
library(data.table)
library(ggplot2); theme_set(theme_bw())
ylims <- c(0.0001, 0.9999)

analytic <- fread('validation/data/stochastic_scan_analytic.csv')
analytic[,BI:=vapply(seq_len(.N),function(i)
  bi_quantities(R0[i],rho[i],theta[i],K[i])$P_conditional,0.)]
analytic[,large_K:=FALSE]

large_analytic <- unique(fread('validation/data/fig11_scan_analytic.csv')[,.(rho,theta,K,R0)])
large_analytic[,BI:=vapply(seq_len(.N),function(i)
  bi_quantities(R0[i],rho[i],theta[i],K[i])$P_conditional,0.)]
large_analytic[,large_K:=TRUE]

dd_analytical <- rbind(analytic[,.(rho,theta,K,R0,BI,large_K)],
  large_analytic[,.(rho,theta,K,R0,BI,large_K)])

sim <- fread('validation/data/stochastic_scan_results.csv')
sim[,large_K:=FALSE]
large_sim <- fread('validation/data/fig11_scan_results.csv')
large_sim[,large_K:=TRUE]

dd_sim <- rbind(sim,large_sim,fill=TRUE)

gg0 <- ggplot(dd_analytical, aes(R0, BI, colour = factor(rho))) +
  geom_line(aes(linetype = factor(theta)), lwd = 1, alpha = 0.7) +
  facet_wrap(~K, labeller = label_both) +
  geom_pointrange(data = dd_sim,
                  aes(y = P_conditional,
                      ymin = cond_low,
                      ymax = cond_high,
                      pch = factor(theta))) +
  scale_colour_okabeito(guide = guide_legend(reverse=TRUE)) +
  scale_y_continuous(transform = "logit", breaks = logit_breaks(),
                     labels = function(b)
                       format(b, trim = TRUE, drop0trailing = TRUE),
                     limits = ylims) +
  scale_x_log10_shifted() +
  labs(x = expression(R[0]), y = "conditional persistence probability")

print(gg0)
ggsave(gg0, file = "validation/figures/fig13_comb.pdf",
       height = 9, width = 12)

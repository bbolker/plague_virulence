source('validation/R/scan.R')
library(data.table)
s <- fread('validation/data/stochastic_results.csv')
for(i in seq_len(nrow(s))) {
  pars <- as.list(s[i,.(R0,rho,theta,K)])
  p <- do.call(run_point,pars)
  p_not_fizzle <- 1-1/s$R0[i]
  s$P1_uncond_K_ODE[i] <- p_not_fizzle*p$P1_ref
  s$P1_uncond_K_first[i] <- p_not_fizzle*p$P1_K_first
  s$P1_uncond_K_second[i] <- p_not_fizzle*p$P1_K_second
  s$P1_uncond_L_second[i] <- p_not_fizzle*p$P1_L_second
}
setcolorder(s,c('R0','rho','theta','K','attempts','established','persistent',
  'P1_hat','CI_low','CI_high','P1_uncond_K_ODE','P1_uncond_K_first',
  'P1_uncond_K_second','P1_uncond_L_second'))
fwrite(s,'validation/data/stochastic_results.csv')

library("tinytest")
using("tinysnapshot")

pars0 <- list(
  n_patches = 100,
  n_years = 500,    
  K = 1e6,
  r = 0.5,
  c0 = 0.2,              
  nu = 5,                
  rho = 3,               
  alpha = 5e-6,          
  R01 = 2.0,             
  R02 = 1.8,            
  invade_year = 100,     
  initial_inf_ratio_1 = 0.1,
  initial_inf_ratio_2 = 0.05
)

source("sim_fun.R")
sim1 <- simulate_metapopulation_2strain(n_patches = 2,
                                n_years  = 10,
                                R01 = 4, R02 = 4,
                                alpha = 3.2e-4,
                                rho = 8,
                                r = 0.5,
                                seed = 101)
expect_snapshot_print(sim1, label = "sim1")

## FIXME: more basic sims (e.g. single strain, very high R0, should persist
##  in a single patch? What goes wrong? No persistence -- do we hit a
## colonization bottleneck?

sim2 <- simulate_metapopulation_2strain(n_patches = 2,
                                n_years  = 1000,
                                R01 = 10, R02 = 4,
                                alpha = 3.2e-4,
                                rho = 8,
                                r = 0.5,
                                seed = 101)

load("polyfit.rda")
p1 <-pred_outcomes_poly(1.5, 2.0, 0.2, 0.1)
p2 <- pred_outcomes_poly(2.0, 1.5, 0.1, 0.2)
expect_equal(p1[1], p2[1])
expect_equal(p1[2], 1-p2[2])
expect_snapshot_print(p1, label = "pred1")

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

source("simulation_funs.R")
sim1 <- simulate_metapopulation_2strain(n_patches = 2,
                                n_years  = 10,
                                R01 = 4, R02 = 4,
                                alpha = 3.2e-4,
                                rho = 8,
                                r = 0.5,
                                seed = 101,
                                coinf_approx = "yy")
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

## Distributional equivalence: C++ and R versions should produce the same
## distribution of outcomes.  Uses two-sample t-tests on log-total-infections
## across 100 replications.  Each test has ~1% chance of spurious failure.
dist_params <- list(
    n_patches = 50, n_years = 150, K = 1e6,
    r = 0.5, c0 = 0.2, nu = 5, rho = 3, alpha = 5e-6, D = 1,
    R01 = 2.0, R02 = 1.8, invade_year = 30,
    initial_inf_ratio_1 = 0.1, initial_inf_ratio_2 = 0.05
)
n_dist_sim <- 100

for (approx in c("polyfit", "yy")) {
    p <- c(dist_params, list(coinf_approx = approx))

    set.seed(42)
    r_inf <- t(vapply(seq_len(n_dist_sim), function(i) {
        res <- do.call(simulate_metapopulation_2strain, p)
        log1p(c(inf1 = sum(res$total_inf1), inf2 = sum(res$total_inf2)))
    }, numeric(2)))

    set.seed(42)
    cpp_inf <- t(vapply(seq_len(n_dist_sim), function(i) {
        res <- do.call(simulate_metapopulation_2strain_cpp, p)
        log1p(c(inf1 = sum(res$total_inf1), inf2 = sum(res$total_inf2)))
    }, numeric(2)))

    for (stat in c("inf1", "inf2")) {
        tt <- t.test(r_inf[, stat], cpp_inf[, stat])
        expect_true(
            tt$p.value > 0.01,
            info = sprintf("coinf_approx=%s %s: distributions differ (p=%.4f, R_mean=%.2f, cpp_mean=%.2f)",
                           approx, stat, tt$p.value,
                           mean(r_inf[, stat]), mean(cpp_inf[, stat]))
        )
    }
}

library(plagueMetapop)

## shared parameters
R0      <- 1.5
K       <- 1e4
r       <- 0.125
gamma   <- 1.0
dt      <- 0.1
nt      <- round(200 / dt)   ## 2000 steps
I_init  <- 10L
nsim    <- 100L

gen1 <- compile_odin("euler_onestrain_odin_def.R")
gen2 <- compile_odin("euler_odin_def.R")

make_mod1 <- function() {
  gen1$new(beta = R0, gamma = gamma, dt = dt,
           I_ini = I_init, S_ini = K - I_init,
           alpha = 0, r = r, K = K, n_patch = 1L)
}

make_mod2 <- function() {
  gen2$new(beta = c(R0, 0), gamma = c(gamma, gamma), dt = dt,
           I_ini = matrix(c(I_init, 0L), nrow = 1L),
           S_ini = K - I_init,
           alpha = 0, strain2_delay = 0L, I2_ini = 0,
           r = r, K = K, n_patch = 1L)
}

## --- structure tests ---

set.seed(42)
out1 <- make_mod1()$run(seq(0L, nt))

expect_true(is.matrix(out1),
            info = "single-strain output is a matrix")
expect_equal(ncol(out1), 3L,
             info = "single-strain output: 3 columns (step, S[1], I[1])")
expect_equal(nrow(out1), nt + 1L,
             info = "single-strain output: nt+1 rows")
expect_equal(colnames(out1), c("step", "S[1]", "I[1]"),
             info = "single-strain column names are correct")

## --- sanity checks ---

expect_true(all(out1[, "S[1]"] >= 0),
            info = "S never negative")
expect_true(all(out1[, "I[1]"] >= 0),
            info = "I never negative")
expect_true(all(out1[, "S[1]"] + out1[, "I[1]"] <= K * 1.1),
            info = "S + I stays within 10% of K")
expect_equal(as.numeric(out1[1L, "step"]), 0,
             info = "first row is step 0 (initial conditions)")
expect_equal(as.numeric(out1[1L, "S[1]"]) + as.numeric(out1[1L, "I[1]"]), K,
             info = "S + I == K at t=0")

## --- benchmark: single-strain vs two-strain (strain 2 absent) ---

run_nsim <- function(make_fn, n, steps) {
  mod <- make_fn()
  for (i in seq_len(n)) mod$run(seq(0L, steps))
  invisible(NULL)
}

set.seed(42); t1 <- system.time(run_nsim(make_mod1, nsim, nt))
set.seed(42); t2 <- system.time(run_nsim(make_mod2, nsim, nt))

cat(sprintf(
  "\nBenchmark (%d sims, n_patch=1, nt=%d, dt=%.1f):\n  euler_onestrain_odin_def: %.2f s\n  euler_odin_def (strain2=0): %.2f s\n  speedup: %.2fx\n",
  nsim, nt, dt, t1["elapsed"], t2["elapsed"],
  t2["elapsed"] / t1["elapsed"]
))

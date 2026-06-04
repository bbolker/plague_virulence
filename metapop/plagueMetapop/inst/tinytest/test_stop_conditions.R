library(plagueMetapop)

## minimal odin-style row: t column + S[1] + I[1,1] + I[1,2]
make_row <- function(I1, I2) {
  m <- matrix(c(1, 100, I1, I2), nrow = 1)
  colnames(m) <- c("t", "S[1]", "I[1,1]", "I[1,2]")
  m
}

## --- stop_either_extinct() unit tests ---

f <- stop_either_extinct()
expect_true(is.function(f),
            info = "stop_either_extinct() returns a function (factory pattern)")

## both strains present: no stop
f <- stop_either_extinct()
expect_false(f(make_row(10, 10)),
             info = "both strains present: returns FALSE")

## strain 1 extinct with strain 2 present: stop
f <- stop_either_extinct()
expect_true(f(make_row(0, 10)),
             info = "strain 1 extinct, strain 2 present: returns TRUE")

## strain 2 zero but never seen (pre-seeding): must NOT stop
f <- stop_either_extinct()
expect_false(f(make_row(10, 0)),
             info = "I2=0 never observed (pre-seeding with strain2_delay): returns FALSE")

## strain 2 seen then extinct, strain 1 still present: stop
f <- stop_either_extinct()
f(make_row(10, 10))           ## strain 2 observed
expect_true(f(make_row(10, 0)),
            info = "I2 was present then went extinct: returns TRUE")

## strain 2 never seen, then strain 1 also goes extinct: stop via strain-1 rule
f <- stop_either_extinct()
f(make_row(10, 0))            ## pre-seeding: strain 2 not seen
expect_true(f(make_row(0, 0)),
            info = "strain 2 never seen, strain 1 goes extinct: returns TRUE")

## both extinct after strain 2 was seen: stop
f <- stop_either_extinct()
f(make_row(10, 10))
expect_true(f(make_row(0, 0)),
            info = "both extinct after strain 2 seen: returns TRUE")

## --- stop_both_extinct() unit tests (regression) ---

expect_false(stop_both_extinct(make_row(10, 10)),
             info = "stop_both_extinct: both present returns FALSE")
expect_false(stop_both_extinct(make_row(0, 10)),
             info = "stop_both_extinct: only strain 1 extinct returns FALSE")
expect_false(stop_both_extinct(make_row(10, 0)),
             info = "stop_both_extinct: only strain 2 extinct returns FALSE")
expect_true(stop_both_extinct(make_row(0, 0)),
            info = "stop_both_extinct: both extinct returns TRUE")

## --- integration test ---
## With strain2_delay=30 and chunk=5, the stop condition is evaluated at steps
## 5, 10, 15, 20, 25, 30, ... Before step 30, I2=0 everywhere.  Without the
## fix, the run would stop at step 5.  With the fix it must run at least to
## step 30 (when strain 2 is seeded).  I_init[1]=50 with K=1e4 makes strain-1
## extinction in 30 steps essentially impossible.  n_patch=1 so I_ini is
## deterministic (round(I_init)), not Poisson.
set.seed(42)
mod <- make_simulator_odin(
  beta_vec      = c(1.5, 2.5),
  K             = 1e4,
  n_patch       = 1L,
  nt            = 500L,
  I_ini_mat     = make_I_ini_mat(c(50, 20), n_patch = 1L, method = "fixed"),
  strain2_delay = 30
)
run <- run_simulator_odin(mod, chunk = 5, stop_cond = stop_either_extinct())
expect_true(max(run[, 1]) >= 30,
            info = "with strain2_delay=30, simulation runs at least to step 30")

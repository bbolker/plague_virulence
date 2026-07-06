library(plagueMetapop)
library(methods)

all.equal.nocheck <- function(x, y, ..., check.attributes = FALSE, check.class = FALSE) {
  if (is(x, "Matrix")) x <- matrix(x)
  if (is(y, "Matrix")) y <- matrix(y)
  all.equal(x, y, ..., check.attributes = check.attributes, check.class = check.class)
}
expect_equal_nocheck <- function(..., tolerance = 5e-5) {
  expect_true(isTRUE(all.equal.nocheck(..., tolerance = tolerance)))
}

## Shared parameters for both runs
common <- list(
  beta_vec      = c(2, 0),
  K             = 500,
  r             = 0.125,
  n_patch       = 3,
  nt            = 30,
  alpha         = 0,
  I_init        = c(5, 0),
  strain2_delay = 0,
  stop_cond     = NULL,
  seed          = 42
)

run_discrete <- do.call(discrete_run,
                        c(common, list(def_file = "discrete_odin_def.R")))

run_euler_rf <- do.call(discrete_run,
                        c(common, list(def_file    = "euler_odin_def.R",
                                       reedfrost   = 1,
                                       gamma       = c(1, 1),
                                       dt          = 1)))

expect_equal_nocheck(run_discrete, run_euler_rf,
                     info = "euler with reedfrost=1, gamma=1, dt=1 matches discrete model")

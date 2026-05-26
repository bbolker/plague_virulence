## HPC array version of discrete_onepatch_twostrain_extinct_run.R
## The full grid (stepR0=0.025) has 25921 rows — over the 1000-job limit.
## We batch: each of 1000 jobs processes ~26 rows.
## Combine outputs with discrete_onepatch_twostrain_extinct_combine.R

library(dplyr)
library(tidyr)
library(odin)
library(dde)
library(optparse)

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) stop("SLURM_ARRAY_TASK_ID not set")

source(here::here("metapop/odin", "discrete_run.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-s", "--stepR0"), type = "double",  default = 0.025,
              help = "R0 grid step size [default %default]"),
  make_option(c("-n", "--nsim"),   type = "integer", default = 20L,
              help = "simulations per grid point [default %default]"),
  make_option(c("-j", "--njobs"),  type = "integer", default = 1000L,
              help = "total number of array jobs [default %default]")
)))

nsim   <- opt$nsim
n_jobs <- opt$njobs

R0vec <- seq(1, 5, by = opt$stepR0)
dd    <- expand.grid(R01 = R0vec, R02 = R0vec)

batch_size <- ceiling(nrow(dd) / n_jobs)
row_start  <- (task_id - 1L) * batch_size + 1L
row_end    <- min(task_id * batch_size, nrow(dd))

if (row_start > nrow(dd)) {
  message("task ", task_id, ": no rows to process (grid smaller than n_jobs)")
  quit(save = "no", status = 0)
}

rows <- dd[seq(row_start, row_end), ]

first_zero <- function(x) which(x == 0)[1]

result_list <- lapply(seq_len(nrow(rows)), function(j) {
  beta1 <- rows$R01[j]; beta2 <- rows$R02[j]
  runs <- discrete_run(beta_vec = c(beta1, beta2),
                       K        = 1e6,
                       r        = 0.125,
                       n_patch  = 1,
                       nt       = 1000,
                       alpha    = 0,
                       I_init   = c(10, 10),
                       stop_cond = NULL,
                       nsim     = nsim,
                       platform = "odin")
  if (!is.list(runs)) runs <- list(runs)
  ext <- sapply(runs, function(traj)
    c(I1 = first_zero(traj$value[traj$state == "I1"]),
      I2 = first_zero(traj$value[traj$state == "I2"])))
  result <- rowMeans(ext, na.rm = TRUE)
  bind_cols(rows[j, ], as.data.frame(as.list(result)))
})

out <- bind_rows(result_list)

dir.create("outputs", showWarnings = FALSE)
saveRDS(out, sprintf("outputs/twostrain_extinct_task_%04d.rds", task_id))

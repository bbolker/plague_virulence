## HPC array version of euler_twostrain_run.R
## Submit via submit_euler_twostrain.sh (full) or submit_euler_twostrain_mini.sh (mini).
## Combine outputs with euler_twostrain_combine.R

library(plagueMetapop)
library(future)
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-m", "--mini"), action = "store_true", default = FALSE,
              help = "run mini test grid (coarser, fewer sims)"),
  make_option("--mini2", action = "store_true", default = FALSE,
              help = "run mini2 grid (R01 x R02 at fixed K=1e4, alpha=10^-5.5; stop_both_extinct)"),
  make_option("--singlepatchintro", action = "store_true", default = FALSE,
              help = "seed strain 2 into patch 1 only (Poisson mean 10); all other patches start at 0")
)))

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) stop("SLURM_ARRAY_TASK_ID not set")

base_fn <- "euler_twostrain"
if (opt$mini)            base_fn <- paste0(base_fn, "_mini")
if (opt$mini2)           base_fn <- paste0(base_fn, "_mini2")
if (opt$singlepatchintro) base_fn <- paste0(base_fn, "_singlepatchintro")

if (opt$mini) {
  dd <- expand.grid(R01   = seq(1.1, 3, by = 0.5),
                    R02   = seq(1.1, 3, by = 0.5),
                    K     = 10^seq(3, 5, by = 1),
                    alpha = 10^seq(-5, -4))
  n_patch   <- 50
  n_sim     <- 50L
  dt        <- 0.2
  stop_cond <- stop_either_extinct()
} else if (opt$mini2) {
  dd <- expand.grid(R01   = seq(1.1, 4, by = 0.1),
                    R02   = seq(1.1, 4, by = 0.1),
                    K     = 1e4,
                    alpha = 10^-5.5)
  n_patch   <- 200
  n_sim     <- 200L
  dt        <- 0.1
  stop_cond <- stop_both_extinct
} else {
  dd <- expand.grid(R01   = seq(1.1, 4, by = 0.1),
                    R02   = seq(1.1, 4, by = 0.1),
                    K     = 10^seq(3, 5),
                    alpha = 10^seq(-5, -3))
  n_patch   <- 200
  n_sim     <- 200L
  dt        <- 0.1
  stop_cond <- stop_both_extinct
}

row           <- dd[task_id, ]
strain2_delay <- round(100 / dt)

## For --singlepatchintro, seed strain 2 only into patch 1 (Poisson mean 10);
## all other patches have lambda=0 so rpois always gives 0.
if (opt$singlepatchintro) {
  n_sim  <- 500L
  I_init <- matrix(c(rep(10, n_patch), c(10, rep(0, n_patch - 1L))),
                   nrow = n_patch, ncol = 2L)
} else {
  I_init <- c(10, 10)
}

plan(multicore(workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1L))))
runs <- discrete_run(beta_vec      = c(row$R01, row$R02),
                     K             = row$K,
                     r             = 0.125,
                     n_patch       = n_patch,
                     nt            = round(500 / dt),
                     alpha         = row$alpha,
                     I_init        = I_init,
                     gamma         = c(1, 1),
                     dt            = dt,
                     def_file      = "euler_odin_def.R",
                     strain2_delay = strain2_delay,
                     stop_cond     = stop_cond,
                     nsim          = n_sim,
                     platform      = "odin")
plan(sequential)

out <- dplyr::bind_cols(row, as.data.frame(as.list(sumfun_discrete(runs, strain2_delay = strain2_delay))))

dir.create("outputs", showWarnings = FALSE)
saveRDS(out, sprintf("outputs/%s_task_%06d.rds", base_fn, task_id))
if (task_id == 1L)
  saveRDS(runs, sprintf("outputs/%s_runs_task_%06d.rds", base_fn, task_id))

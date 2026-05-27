## HPC array version of euler_onestrain_run.R
## Submit via submit_euler_onestrain.sh (full), submit_euler_onestrain_mini.sh (mini),
##   or submit_euler_onestrain_batch2.sh (batch2: higher-alpha extension).
## Combine outputs with euler_onestrain_combine.R

library(plagueMetapop)
library(future)
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-m", "--mini"), action = "store_true", default = FALSE,
              help = "run mini test grid (coarser, fewer sims)"),
  make_option(c("-b", "--batch2"), action = "store_true", default = FALSE,
              help = "batch 2: same full grid but log10(alpha) = seq(-3, -2, by=0.5)")
)))

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) stop("SLURM_ARRAY_TASK_ID not set")

base_fn <- "euler_onestrain"
if (opt$mini)   base_fn <- paste0(base_fn, "_mini")
if (opt$batch2) base_fn <- paste0(base_fn, "_batch2")

if (opt$mini) {
  dd <- expand.grid(R0    = seq(1.1, 3, by = 0.5),
                    K     = 10^seq(4, 6, by = 1),
                    alpha = 10^seq(-5, -4))
  n_patch <- 50
  n_sim   <- 50L
  dt      <- 0.2
} else if (opt$batch2) {
  dd <- expand.grid(R0    = seq(1.1, 5, by = 0.1),
                    K     = 10^seq(3, 6, by = 0.5),
                    alpha = 10^seq(-3, -2, by = 0.5))
  n_patch <- 200
  n_sim   <- 200L
  dt      <- 0.1
} else {
  dd <- expand.grid(R0    = seq(1.1, 5, by = 0.1),
                    K     = 10^seq(3, 6, by = 0.5),
                    alpha = 10^seq(-5.5, -3.5, by = 0.5))
  n_patch <- 200
  n_sim   <- 200L
  dt      <- 0.1
}

row <- dd[task_id, ]

plan(multicore(workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1L))))
runs <- discrete_run(beta_vec  = c(row$R0, 0),
                     K         = row$K,
                     r         = 0.125,
                     n_patch   = n_patch,
                     nt        = round(200 / dt),
                     alpha     = row$alpha,
                     I_init    = c(10, 0),
                     gamma     = c(1, 1),
                     dt        = dt,
                     def_file  = "euler_odin_def.R",
                     stop_cond = NULL,
                     nsim      = n_sim,
                     platform  = "odin")
plan(sequential)

out <- dplyr::bind_cols(row, as.data.frame(as.list(sumfun_discrete(runs))))

dir.create("outputs", showWarnings = FALSE)
saveRDS(out, sprintf("outputs/%s_task_%04d.rds", base_fn, task_id))

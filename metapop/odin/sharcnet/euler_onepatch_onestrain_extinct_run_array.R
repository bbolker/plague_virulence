## HPC array version of euler_onepatch_onestrain_extinct_run.R
## Submit via submit_euler_extinct.sh; each task runs one (R0, K) grid point.
## Combine outputs with euler_onepatch_onestrain_extinct_combine.R

library(plagueMetapop)
library(future)

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) stop("SLURM_ARRAY_TASK_ID not set")

nsim <- 100L
dt   <- 0.1

dd <- expand.grid(R0 = seq(1.1, 5, by = 0.1),
                  K  = 10^(seq(3, 6, by = 0.5)))

row <- dd[task_id, ]

plan(multicore(workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1L))))
runs <- discrete_run(beta_vec  = c(row$R0, 0),
                     K         = row$K,
                     r         = 0.125,
                     n_patch   = 1,
                     nt        = round(200 / dt),
                     alpha     = 0,
                     I_init    = c(10, 0),
                     gamma     = c(1, 1),
                     dt        = dt,
                     def_file  = "euler_odin_def.R",
                     stop_cond = NULL,
                     nsim      = nsim,
                     platform  = "odin")
plan(sequential)

out <- dplyr::bind_cols(row, as.data.frame(as.list(sumfun_discrete(runs))))

dir.create("outputs", showWarnings = FALSE)
saveRDS(out, sprintf("outputs/euler_onepatch_onestrain_extinct_task_%04d.rds", task_id))

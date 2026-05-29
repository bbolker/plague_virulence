## HPC array version of euler_onepatch_onestrain_extinct_run.R
## Submit via submit_euler_extinct.sh; each task runs one (R0, K) grid point.
## Combine outputs with euler_onepatch_onestrain_extinct_combine.R

library(plagueMetapop)
library(future)
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-l", "--lineargrowth"), action = "store_true", default = FALSE,
              help = "use linear restoring force demography (logistic_growth=0)"),
  make_option(c("-r", "--reedfrost"), action = "store_true", default = FALSE,
              help = "use Reed-Frost (100%% removal per step) dynamics")
)))

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) stop("SLURM_ARRAY_TASK_ID not set")

logistic_growth <- if (opt$lineargrowth) 0 else 1
reedfrost       <- if (opt$reedfrost) 1 else 0
dt              <- if (reedfrost == 1) 1 else 0.1
nsim            <- 1000L

base_fn <- paste("euler_onepatch_onestrain_extinct",
                 if (opt$lineargrowth) "linear" else "logistic",
                 if (opt$reedfrost) "reedfrost" else "continuous",
                 sep = "_")

dd <- expand.grid(R0 = seq(1.1, 5, by = 0.1),
                  K  = 10^seq(3, 6, by = 0.25))

row <- dd[task_id, ]

plan(multicore(workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1L))))
runs <- discrete_run(beta_vec        = c(row$R0, 0),
                     K               = row$K,
                     n_patch         = 1,
                     nt              = round(200 / dt),
                     alpha           = 0,
                     I_init          = c(10, 0),
                     gamma           = c(1, 1),
                     dt              = dt,
                     logistic_growth = logistic_growth,
                     reedfrost       = reedfrost,
                     def_file        = "euler_odin_def.R",
                     stop_cond       = stop_both_extinct,
                     nsim            = nsim,
                     platform        = "odin")
plan(sequential)

out <- dplyr::bind_cols(row, as.data.frame(as.list(sumfun_discrete(runs))))

dir.create("outputs", showWarnings = FALSE)
saveRDS(out, sprintf("outputs/%s_task_%04d.rds", base_fn, task_id))

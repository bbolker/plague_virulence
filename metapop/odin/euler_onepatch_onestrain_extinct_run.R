library(plagueMetapop)
library(future)
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-m", "--mini"), action = "store_true", default = FALSE,
              help = "run a small test grid (coarser, fewer sims)"),
  make_option(c("-l", "--lineargrowth"), action = "store_true", default = FALSE,
              help = "use linear restoring force demography (logistic_growth=0)"),
  make_option(c("-r", "--reedfrost"), action = "store_true", default = FALSE,
              help = "use Reed-Frost (100%% removal per step) dynamics")
)))

logistic_growth <- if (opt$lineargrowth) 0 else 1
reedfrost       <- if (opt$reedfrost) 1 else 0
dt              <- if (reedfrost == 1) 1 else 0.1
nsim            <- 100L

base_fn <- paste("euler_onepatch_onestrain_extinct",
                 if (opt$lineargrowth) "linear" else "logistic",
                 if (opt$reedfrost) "reedfrost" else "continuous",
                 sep = "_")
if (opt$mini) base_fn <- paste0(base_fn, "_mini")

if (opt$mini) {
  dd <- expand.grid(R0 = seq(1.1, 3, by = 0.25),
                    K  = 10^seq(4, 6, by = 1))
} else {
  dd <- expand.grid(R0 = seq(1.1, 5, by = 0.1),
                    K  = 10^seq(3, 6, by = 0.5))
}

plan(multicore(workers = 14L))
set.seed(101)

result_list <- lapply(seq_len(nrow(dd)), function(i) {
  cat(dd$R0[i], dd$K[i], "\n")
  runs <- discrete_run(beta_vec        = c(dd$R0[i], 0),
                       K               = dd$K[i],
                       n_patch         = 1,
                       nt              = round(200 / dt),
                       alpha           = 0,
                       I_init          = c(10, 0),
                       gamma           = c(1, 1),
                       dt              = dt,
                       logistic_growth = logistic_growth,
                       reedfrost       = reedfrost,
                       def_file        = "euler_odin_def.R",
                       stop_cond       = NULL,
                       nsim            = nsim,
                       platform        = "odin")
  dplyr::bind_cols(dd[i, ], as.data.frame(as.list(sumfun_discrete(runs))))
})
plan(sequential)

out <- dplyr::bind_rows(result_list)
saveRDS(out, paste0(base_fn, ".rds"))

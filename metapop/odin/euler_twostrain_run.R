library(plagueMetapop)
library(future)
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-m", "--mini"), action = "store_true", default = FALSE,
              help = "run a small test grid (coarser, fewer sims)")
)))

base_fn <- "euler_twostrain"
if (opt$mini) base_fn <- paste0(base_fn, "_mini")

if (opt$mini) {
  dd <- expand.grid(R01   = seq(1.1, 3, by = 0.5),
                    R02   = seq(1.1, 3, by = 0.5),
                    K     = 10^seq(3, 5, by = 1),
                    alpha = 10^seq(-5, -4))
  n_patch <- 50
  n_sim   <- 50L
  dt      <- 0.2
} else {
  dd <- expand.grid(R01   = seq(1.1, 5, by = 0.1),
                    R02   = seq(1.1, 5, by = 0.1),
                    K     = 10^seq(3, 5),
                    alpha = 10^seq(-5, -3))
  n_patch <- 200
  n_sim   <- 200L
  dt      <- 0.1
}

plan(multicore(workers = 14L))
set.seed(101)

result_list <- lapply(seq_len(nrow(dd)), function(i) {
  cat(dd$R01[i], dd$R02[i], dd$K[i], dd$alpha[i], "\n")
  runs <- discrete_run(beta_vec      = c(dd$R01[i], dd$R02[i]),
                       K             = dd$K[i],
                       r             = 0.125,
                       n_patch       = n_patch,
                       nt            = round(200 / dt),
                       alpha         = dd$alpha[i],
                       I_init        = c(10, 10),
                       gamma         = c(1, 1),
                       dt            = dt,
                       def_file      = "euler_odin_def.R",
                       strain2_delay = round(100 / dt),
                       stop_cond     = stop_either_extinct(),
                       nsim          = n_sim,
                       platform      = "odin")
  dplyr::bind_cols(dd[i, ], as.data.frame(as.list(sumfun_discrete(runs))))
})
plan(sequential)

out <- dplyr::bind_rows(result_list)
saveRDS(out, paste0(base_fn, ".rds"))

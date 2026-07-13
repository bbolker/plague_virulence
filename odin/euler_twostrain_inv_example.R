library(plagueMetapop)
library(future)
library(dplyr)
library(here)
library(ggplot2); theme_set(theme_bw())

## lower middle panel of PIP plot,yellow area (resident fine, invader loses)
R01 <- 2
R02 <- 3.8
K <- 1e3
alpha <- 1e-3
n_patch   <- 200
#n_sim     <- 100L
n_sim <- 20L
dt        <- 0.1
stop_cond <- stop_both_extinct()
strain2_delay <- round(100 / dt)
nt <- round(500/dt)
I_init <- cbind(rep(10, n_patch),
                c(10, rep(0, n_patch - 1L)))

# plan(multicore(workers=20))
plan(sequential)
set.seed(101)
runs <- discrete_run(beta_vec      = c(R01, R02),
                     K             = K,
                     r             = 0.125,
                     n_patch       = n_patch,
                     nt            = round(500 / dt),
                     alpha         = alpha,
                     I_init        = I_init,
                     gamma         = c(1, 1),
                     dt            = dt,
                     def_file      = "euler_odin_def.R",
                     strain2_delay = strain2_delay,
                     stop_cond     = stop_cond,
                     chunk         = strain2_delay + 1L,
                     nsim          = n_sim,
                     platform      = "odin")

runs_x <- sum_runs(runs) |> dplyr::filter(var == "pop") |>
  mutate(across(value, ~ . /n_patch))
# 
# ggplot(runs_x, aes(step, value, colour = type)) +
#   geom_line(aes(group = interaction(run, type)), alpha = 0.4) +
#   scale_y_log10()
# 
# sumfun_discrete(runs, strain2_delay = strain2_delay)

p <- ggplot(runs_x, aes(step, value, colour = type)) +
  geom_line(aes(group = interaction(run, type)), alpha = 0.4) +
  scale_y_log10()

print(p)

fmt <- function(x) {
  x <- format(x, scientific = FALSE, trim = TRUE)
  x <- gsub("\\.", "p", x)
  x <- gsub("-", "m", x)
  x
}

out_file <- sprintf(
  "euler_twostrain_inv_example_R01_%s_R02_%s_K_%s_alpha_%s.pdf",
  fmt(R01), fmt(R02), fmt(K), fmt(alpha)
)

outdir <- here::here("odin/outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  file.path(outdir, out_file),
  p,
  width = 8,
  height = 5,
  device = "pdf"
)

message("Saved plot to: ", file.path(outdir, out_file))

sumfun_discrete(runs, strain2_delay = strain2_delay)
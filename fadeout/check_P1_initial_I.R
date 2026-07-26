## Check whether the single-patch persistence probability P1 changes strongly
## when the mean initial infected count is increased from 10 toward I*.
## Run from the repository root:
##   Rscript fadeout/check_P1_initial_I.R

library(future)
library(ggplot2)
library(here)
library(mgcv)
library(plagueMetapop)

R0 <- 2.5
K <- 10000
r <- 0.125
gamma <- 1
dt <- 0.1
t_max <- 200
nsim <- 1000L
I0_values <- c(10, 30, 100, 300)

single_patch_file <- here(
  "odin", "sharcnet", "outputs",
  "euler_onepatch_onestrain_extinct_logistic_continuous.rds"
)
output_file <- here("fadeout", "output", "P1_initial_I_effect.pdf")

if (!file.exists(single_patch_file)) {
  stop("Required single-patch data not found: ", single_patch_file)
}
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

single_patch <- readRDS(single_patch_file)
required_columns <- c("R0", "K", "ext_prob.I1")
if (!all(required_columns %in% names(single_patch))) {
  stop("Single-patch data must contain: ",
       paste(required_columns, collapse = ", "))
}

gam_fit <- mgcv::gam(
  ext_prob.I1 ~ te(R0, K, k = c(12, 12)),
  data = single_patch,
  method = "REML"
)
P1_gam <- 1 - as.numeric(predict(
  gam_fit,
  newdata = data.frame(R0 = R0, K = K),
  type = "response"
))

I_star <- unname(plagueMetapop::ode_eq(
  beta = R0, gamma = gamma, K = K, r = r, logistic_growth = 1
)[["eq_I"]])
I0_values[length(I0_values)] <- round(I_star)

## Match the existing P1 calibration workflow: Poisson initial infected count,
## logistic continuous-time Euler model, and extinction assessed by t = 200.
## Keep this sequential: on Windows an odin compiled object cannot reliably be
## serialized to multisession workers. The single-patch runs remain lightweight.
future::plan(future::sequential)

estimate_one <- function(I0, index) {
  runs <- plagueMetapop::discrete_run(
    beta_vec = c(R0, 0),
    K = K,
    r = r,
    n_patch = 1,
    nt = round(t_max / dt),
    alpha = 0,
    I_init = c(I0, 0),
    I_ini_method = "rpois",
    gamma = c(gamma, gamma),
    dt = dt,
    logistic_growth = 1,
    reedfrost = 0,
    def_file = "euler_odin_def.R",
    stop_cond = plagueMetapop::stop_both_extinct(require_seeded = FALSE),
    nsim = nsim,
    seed = 1200L + index,
    platform = "odin"
  )
  persistence <- 1 - unname(
    plagueMetapop::sumfun_discrete(runs)[["ext_prob.I1"]]
  )
  survived <- round(persistence * nsim)
  interval <- stats::binom.test(survived, nsim)$conf.int
  data.frame(
    I0 = I0,
    P1 = persistence,
    lower = interval[1],
    upper = interval[2]
  )
}

results <- do.call(
  rbind,
  Map(estimate_one, I0_values, seq_along(I0_values))
)

P1_at_10 <- results$P1[results$I0 == 10]
P1_at_Istar <- results$P1[results$I0 == round(I_star)]
difference <- P1_at_Istar - P1_at_10

annotation <- sprintf(
  "P1 at E[I(0)]=10: %.3f\nP1 at E[I(0)]=I*=%.0f: %.3f\nDifference: %+.3f",
  P1_at_10, I_star, P1_at_Istar, difference
)

p <- ggplot(results, aes(I0, P1)) +
  geom_hline(
    yintercept = P1_gam, colour = "#D55E00", linewidth = 0.7,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper), width = 0.08,
    colour = "#0072B2", linewidth = 0.6
  ) +
  geom_line(colour = "#0072B2", linewidth = 0.7) +
  geom_point(colour = "#0072B2", size = 2.5) +
  annotate(
    "label", x = max(I0_values), y = 0.16,
    label = annotation, hjust = 1, vjust = 0.5, size = 3.6
  ) +
  scale_x_log10(breaks = I0_values) +
  coord_cartesian(ylim = c(0, max(results$upper, P1_gam) + 0.035)) +
  labs(
    title = "Initial infected count has limited effect on P1",
    subtitle = sprintf(
      "R0 = %.1f; K = %g; r = %.3f; gamma = %g; %d simulations per point",
      R0, K, r, gamma, nsim
    ),
    x = "Mean initial infected count, E[I(0)] (Poisson)",
    y = "Persistence probability, P1",
    caption = sprintf(
      "Points: persistence to t = %g (95%% binomial CI). Dashed line: existing GAM prediction P1 = %.3f.",
      t_max, P1_gam
    )
  ) +
  theme_bw(base_size = 11) +
  theme(plot.caption = element_text(hjust = 0))

ggsave(output_file, p, width = 8, height = 5.4)

cat(sprintf(
  paste0(
    "Parameters: R0=%.1f, K=%g, r=%.3f, gamma=%g, I*=%.1f\n",
    "Existing GAM P1: %.3f\n",
    "Direct P1 at E[I(0)]=10: %.3f\n",
    "Direct P1 at E[I(0)]=I*: %.3f\n",
    "Difference (I* minus 10): %+.3f\n",
    "Output: %s\n"
  ),
  R0, K, r, gamma, I_star, P1_gam,
  P1_at_10, P1_at_Istar, difference, output_file
))

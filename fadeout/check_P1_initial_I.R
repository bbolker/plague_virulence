## Check whether the single-patch persistence probability P1 changes strongly
## when the mean initial infected count is increased from 10 toward I*.
## Run from the repository root:
##   Rscript fadeout/check_P1_initial_I.R

library(future)
library(ggplot2)
library(here)
library(plagueMetapop)

R0 <- 2.5
K <- 10000
r <- 0.125
gamma <- 1
dt <- 0.1
t_max <- 200
nsim <- 1000L
I0_values <- c(1, 2, 3, 4, 5, 10, 30, 100, 300)

output_dir <- here("fadeout", "output", "P1_initial_I")
cache_file <- file.path(output_dir, "P1_initial_I_results.csv")
output_file <- file.path(output_dir, "P1_by_initial_I_to300.pdf")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Match the existing P1 calibration workflow: Poisson initial infected count,
## logistic continuous-time Euler model, and extinction assessed by t = 200.
## Keep this sequential: on Windows an odin compiled object cannot reliably be
## serialized to multisession workers. The single-patch runs remain lightweight.
future::plan(future::sequential)

estimate_one <- function(I0) {
  seed <- 120000L + as.integer(I0)
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
    seed = seed,
    platform = "odin"
  )
  persistence <- 1 - unname(
    plagueMetapop::sumfun_discrete(runs)[["ext_prob.I1"]]
  )
  survived <- round(persistence * nsim)
  interval <- stats::binom.test(survived, nsim)$conf.int
  data.frame(
    I0 = I0,
    nsim = nsim,
    survived = survived,
    extinct = nsim - survived,
    P1 = persistence,
    lower = interval[1],
    upper = interval[2],
    seed = seed,
    R0 = R0,
    K = K,
    r = r,
    gamma = gamma,
    dt = dt,
    t_max = t_max,
    initialization = "Poisson"
  )
}

cache_columns <- c(
  "I0", "nsim", "survived", "extinct", "P1", "lower", "upper", "seed",
  "R0", "K", "r", "gamma", "dt", "t_max", "initialization"
)
cache <- if (file.exists(cache_file)) {
  cached <- read.csv(cache_file, stringsAsFactors = FALSE)
  missing <- setdiff(cache_columns, names(cached))
  if (length(missing)) {
    stop("Cache is missing required columns: ", paste(missing, collapse = ", "))
  }
  cached[cache_columns]
} else {
  as.data.frame(setNames(replicate(
    length(cache_columns), logical(0), simplify = FALSE
  ), cache_columns))
}

is_compatible <- function(dat, I0) {
  if (!nrow(dat)) return(logical(0))
  dat$I0 == I0 &
    dat$nsim == nsim &
    dat$seed == 120000L + as.integer(I0) &
    dat$R0 == R0 &
    dat$K == K &
    dat$r == r &
    dat$gamma == gamma &
    dat$dt == dt &
    dat$t_max == t_max &
    dat$initialization == "Poisson"
}

result_list <- vector("list", length(I0_values))
new_results <- list()
for (index in seq_along(I0_values)) {
  I0 <- I0_values[index]
  cached_rows <- cache[is_compatible(cache, I0), , drop = FALSE]
  if (nrow(cached_rows)) {
    result_list[[index]] <- cached_rows[nrow(cached_rows), , drop = FALSE]
    cat("Reusing cached result for E[I(0)] = ", I0, "\n", sep = "")
  } else {
    cat("Simulating missing result for E[I(0)] = ", I0, "\n", sep = "")
    result_list[[index]] <- estimate_one(I0)
    new_results[[length(new_results) + 1L]] <- result_list[[index]]
  }
}
results <- do.call(rbind, result_list)
results <- results[order(results$I0), , drop = FALSE]

if (length(new_results)) {
  cache <- rbind(cache, do.call(rbind, new_results))
  cache <- cache[!duplicated(
    cache[c(
      "I0", "nsim", "seed", "R0", "K", "r", "gamma", "dt", "t_max",
      "initialization"
    )],
    fromLast = TRUE
  ), , drop = FALSE]
  cache <- cache[order(cache$R0, cache$K, cache$r, cache$gamma,
                       cache$I0), , drop = FALSE]
  write.csv(cache, cache_file, row.names = FALSE)
}

p <- ggplot(results, aes(I0, P1)) +
  geom_linerange(
    aes(ymin = lower, ymax = upper),
    colour = "#0072B2", linewidth = 0.6
  ) +
  geom_line(colour = "#0072B2", linewidth = 0.7) +
  geom_point(colour = "#0072B2", size = 2.5) +
  scale_x_log10(breaks = I0_values) +
  coord_cartesian(ylim = c(0, max(results$upper) + 0.035)) +
  labs(
    title = "Sensitivity of P1 to the initial infected count",
    subtitle = sprintf(
      "R0 = %.1f; K = %g; r = %.3f; gamma = %g; %d simulations per point",
      R0, K, r, gamma, nsim
    ),
    x = "Mean initial infected count, E[I(0)] (Poisson)",
    y = "Persistence probability, P1",
    caption = sprintf(
      "Points: persistence to t = %g (95%% binomial CI).",
      t_max
    )
  ) +
  theme_bw(base_size = 11) +
  theme(plot.caption = element_text(hjust = 0))

ggsave(output_file, p, width = 8, height = 5.4)

cat(sprintf(
  paste0(
    "Parameters: R0=%.1f, K=%g, r=%.3f, gamma=%g\n",
    "Data: %s\n",
    "Figure: %s\n"
  ),
  R0, K, r, gamma, cache_file, output_file
))

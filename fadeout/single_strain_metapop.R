library(plagueMetapop)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

theme_set(theme_bw())

## ------------------------------------------------------------
## Parameters
## ------------------------------------------------------------

R0 <- 3
K <- 1e4
alpha <- 1e-4

r <- 0.125
gamma <- 1

n_patch <- 200
dt <- 0.1
t_max <- 500


seed <- 101

## ------------------------------------------------------------
## Run one stochastic trajectory
## ------------------------------------------------------------

# set.seed(seed)
# 
# runs <- discrete_run(
#   beta_vec  = c(R0, 0),
#   K         = K,
#   r         = r,
#   n_patch   = n_patch,
#   nt        = round(t_max / dt),
#   alpha     = alpha,
#   I_init    = c(I0, 0),
#   gamma     = c(gamma, gamma),
#   dt        = dt,
#   def_file  = "euler_odin_def.R",
#   stop_cond = NULL,
#   nsim      = 1,
#   platform  = "odin"
# )


set.seed(seed)

## Compute deterministic endemic equilibrium

eq <- ode_eq(
  beta = R0,
  gamma = gamma,
  K = K,
  r = r,
  logistic_growth = 1
)

S_star <- unname(eq["eq_S"])
I_star <- unname(eq["eq_I"])

cat("Deterministic equilibrium:\n")
cat("S* =", S_star, "\n")
cat("I* =", I_star, "\n")

## Every patch starts exactly at the deterministic equilibrium

# S_ini <- rep(round(S_star), n_patch)
#
# I_ini <- cbind(
#   rep(round(I_star), n_patch),
#   rep(0, n_patch)
# )

I_outbreak <- 10

u <- runif(n_patch, min = 0, max = 1)

S_ini <- round(
  K * u + S_star * (1 - u)
)

I1_ini <- round(
  I_outbreak * u + I_star * (1 - u)
)

I_ini <- cbind(
  I1_ini,
  rep(0, n_patch)
)

## Compile and run the same odin model with explicit S and I initial values

gen <- compile_odin("euler_odin_def.R")

sim <- gen$new(
  beta = c(R0, 0),
  gamma = c(gamma, gamma),
  dt = dt,
  r = rep(r, n_patch),
  K = rep(K, n_patch),
  S_ini = S_ini,
  I_ini = I_ini,
  I2_ini = rep(0, n_patch),
  alpha = alpha,
  n_patch = n_patch,
  strain2_delay = as.integer(.Machine$integer.max),
  logistic_growth = 1,
  reedfrost = 0
)

raw_runs <- sim$run(seq.int(0L, round(t_max / dt)))
raw_runs[, "step"] <- raw_runs[, "step"] * dt
runs <- conv_odin(raw_runs)

## Verify actual initial values

initial_check <- runs |>
  filter(step == 0, state %in% c("S", "I1")) |>
  group_by(state) |>
  summarise(
    min = min(value),
    mean = mean(value),
    max = max(value),
    .groups = "drop"
  )

print(initial_check)

## ------------------------------------------------------------
## Extract strain-1 prevalence
## ------------------------------------------------------------

infected <- runs |>
  filter(state == "I1") |>
  select(step, patch, I = value) |>
  arrange(patch, step)

## Binary patch infection state:
## occupied = 1 if I > 0, otherwise 0

infected <- infected |>
  mutate(occupied = as.integer(I > 0))

## ------------------------------------------------------------
## Identify local extinction and recolonization events
## ------------------------------------------------------------

infected <- infected |>
  group_by(patch) |>
  arrange(step, .by_group = TRUE) |>
  mutate(
    prev_occupied = lag(occupied),
    
    local_extinction = as.integer(
      prev_occupied == 1 & occupied == 0
    ),
    
    recolonization = as.integer(
      prev_occupied == 0 & occupied == 1
    )
  ) |>
  ungroup() |>
  mutate(
    local_extinction = replace_na(local_extinction, 0L),
    recolonization = replace_na(recolonization, 0L)
  )

## ------------------------------------------------------------
## Time-dependent metapopulation summaries
## ------------------------------------------------------------

meta_summary <- infected |>
  group_by(step) |>
  summarise(
    occupied_patches = sum(occupied),
    global_I = sum(I),
    local_extinctions = sum(local_extinction),
    recolonizations = sum(recolonization),
    .groups = "drop"
  )

## ------------------------------------------------------------
## Overall summary statistics
## ------------------------------------------------------------

overall_summary <- meta_summary |>
  summarise(
    globally_persistent = all(global_I > 0),
    min_occupied_patches = min(occupied_patches),
    mean_occupied_patches = mean(occupied_patches),
    max_occupied_patches = max(occupied_patches),
    total_local_extinctions = sum(local_extinctions),
    total_recolonizations = sum(recolonizations)
  )

print(overall_summary)

## Patch-level turnover statistics

patch_summary <- infected |>
  group_by(patch) |>
  summarise(
    extinction_events = sum(local_extinction),
    recolonization_events = sum(recolonization),
    fraction_time_infected = mean(occupied),
    .groups = "drop"
  )

print(summary(patch_summary$extinction_events))
print(summary(patch_summary$recolonization_events))

## ------------------------------------------------------------
## Plot 1: patch x time occupancy raster
## ------------------------------------------------------------

p_raster <- ggplot(
  infected,
  aes(
    x = step,
    y = factor(patch),
    fill = factor(occupied)
  )
) +
  geom_raster() +
  scale_fill_manual(
    values = c("0" = "white", "1" = "black"),
    breaks = c("0", "1"),
    labels = c("0" = "uninfected", "1" = "infected"),
    drop = FALSE,
    name = NULL
  ) +
  labs(
    x = "Time",
    y = "Patch",
    title = sprintf(
      "Single-strain metapopulation dynamics: R0 = %g, K = %g, alpha = %g",
      R0, K, alpha
    )
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top"
  )

print(p_raster)

## ------------------------------------------------------------
## Plot 2: number of occupied patches
## ------------------------------------------------------------

p_occupancy <- ggplot(
  meta_summary,
  aes(x = step, y = occupied_patches)
) +
  geom_line() +
  labs(
    x = "Time",
    y = "Number of infected patches",
    title = "Metapopulation occupancy"
  )

print(p_occupancy)

## ------------------------------------------------------------
## Plot 3: local extinction and recolonization events
## ------------------------------------------------------------

event_summary <- meta_summary |>
  select(
    step,
    local_extinctions,
    recolonizations
  ) |>
  pivot_longer(
    cols = c(local_extinctions, recolonizations),
    names_to = "event",
    values_to = "count"
  )

p_events <- ggplot(
  event_summary,
  aes(
    x = step,
    y = count,
    colour = event
  )
) +
  geom_line() +
  labs(
    x = "Time",
    y = "Number of patches",
    colour = NULL,
    title = "Local extinction and recolonization"
  )

print(p_events)

## ------------------------------------------------------------
## Save outputs
## ------------------------------------------------------------

outdir <- here::here("odin", "outputs", "onestrain_rescue")

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  file.path(outdir, "occupancy_raster.pdf"),
  p_raster,
  width = 10,
  height = 7
)

ggsave(
  file.path(outdir, "occupied_patches.pdf"),
  p_occupancy,
  width = 8,
  height = 5
)

ggsave(
  file.path(outdir, "extinction_recolonization.pdf"),
  p_events,
  width = 8,
  height = 5
)

saveRDS(
  runs,
  file.path(outdir, "trajectory.rds")
)

write.csv(
  meta_summary,
  file.path(outdir, "metapop_summary.csv"),
  row.names = FALSE
)

write.csv(
  patch_summary,
  file.path(outdir, "patch_summary.csv"),
  row.names = FALSE
)

write.csv(
  overall_summary,
  file.path(outdir, "overall_summary.csv"),
  row.names = FALSE
)

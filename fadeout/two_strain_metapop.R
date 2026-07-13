library(plagueMetapop)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

theme_set(theme_bw())

## ------------------------------------------------------------
## Parameters
## ------------------------------------------------------------

R01 <- 3
R02 <- 2.5

K <- 1e4
alpha <- 1e-4

r <- 0.125
gamma <- 1

n_patch <- 200

dt <- 0.1

## Internal Euler step:
## step 1000 corresponds to time 100 disease generations
## step 5000 corresponds to time 500 disease generations

strain2_delay <- 1000L
nt <- 5000L

seed <- 101

## ------------------------------------------------------------
## Strain-1 deterministic endemic equilibrium
## ------------------------------------------------------------

eq1 <- ode_eq(
  beta = R01,
  gamma = gamma,
  K = K,
  r = r,
  logistic_growth = 1
)

S1_star <- unname(eq1["eq_S"])
I1_star <- unname(eq1["eq_I"])

cat("Strain-1 deterministic equilibrium:\n")
cat("S* =", S1_star, "\n")
cat("I* =", I1_star, "\n")

## ------------------------------------------------------------
## Initial infection matrix
##
## Strain 1:
##   every patch has Poisson mean I1_star
##
## Strain 2:
##   patch 1 has Poisson mean 10
##   all other patches have mean 0
##
## Strain 2 is not present at time 0.
## The second column is used as the delayed seed at step 1000.
## ------------------------------------------------------------

I_init <- cbind(
  rep(I1_star, n_patch),
  c(10, rep(0, n_patch - 1L))
)

## ------------------------------------------------------------
## Run one stochastic trajectory
## ------------------------------------------------------------

set.seed(seed)

runs <- discrete_run(
  beta_vec      = c(R01, R02),
  K             = K,
  r             = r,
  n_patch       = n_patch,
  nt            = nt,
  alpha         = alpha,
  I_init        = I_init,
  gamma         = c(gamma, gamma),
  dt            = dt,
  def_file      = "euler_odin_def.R",
  strain2_delay = strain2_delay,
  stop_cond     = NULL,
  nsim          = 1,
  platform      = "odin"
)

## ------------------------------------------------------------
## Extract S, I1, and I2
## ------------------------------------------------------------

patch_state <- runs |>
  filter(state %in% c("S", "I1", "I2")) |>
  select(step, patch, state, value) |>
  pivot_wider(
    names_from = state,
    values_from = value
  ) |>
  arrange(patch, step) |>
  mutate(
    N = S + I1 + I2
  )

## ------------------------------------------------------------
## Four patch infection states
##
## white: neither strain
## blue : strain 1 only
## red  : strain 2 only
## black: both strains
## ------------------------------------------------------------

patch_state <- patch_state |>
  mutate(
    infection_state = case_when(
      I1 == 0 & I2 == 0 ~ "Neither",
      I1 > 0  & I2 == 0 ~ "Strain 1 only",
      I1 == 0 & I2 > 0  ~ "Strain 2 only",
      I1 > 0  & I2 > 0  ~ "Both"
    ),
    infection_state = factor(
      infection_state,
      levels = c(
        "Neither",
        "Strain 1 only",
        "Strain 2 only",
        "Both"
      )
    )
  )

## ------------------------------------------------------------
## Patch-level extinction and recolonization transitions
## ------------------------------------------------------------

patch_state <- patch_state |>
  group_by(patch) |>
  arrange(step, .by_group = TRUE) |>
  mutate(
    occupied_I1 = as.integer(I1 > 0),
    occupied_I2 = as.integer(I2 > 0),
    
    prev_I1 = lag(occupied_I1),
    prev_I2 = lag(occupied_I2),
    
    I1_local_extinction = as.integer(
      prev_I1 == 1 & occupied_I1 == 0
    ),
    
    I1_recolonization = as.integer(
      prev_I1 == 0 & occupied_I1 == 1
    ),
    
    I2_local_extinction = as.integer(
      prev_I2 == 1 & occupied_I2 == 0
    ),
    
    I2_recolonization = as.integer(
      prev_I2 == 0 & occupied_I2 == 1
    )
  ) |>
  ungroup() |>
  mutate(
    across(
      c(
        I1_local_extinction,
        I1_recolonization,
        I2_local_extinction,
        I2_recolonization
      ),
      ~ replace_na(.x, 0L)
    )
  )

## ------------------------------------------------------------
## Metapopulation summaries over time
## ------------------------------------------------------------

meta_summary <- patch_state |>
  group_by(step) |>
  summarise(
    patches_neither = sum(I1 == 0 & I2 == 0),
    patches_I1_only = sum(I1 > 0 & I2 == 0),
    patches_I2_only = sum(I1 == 0 & I2 > 0),
    patches_both = sum(I1 > 0 & I2 > 0),
    
    occupied_I1 = sum(I1 > 0),
    occupied_I2 = sum(I2 > 0),
    
    global_I1 = sum(I1),
    global_I2 = sum(I2),
    
    mean_N = mean(N),
    
    I1_local_extinctions = sum(I1_local_extinction),
    I1_recolonizations = sum(I1_recolonization),
    
    I2_local_extinctions = sum(I2_local_extinction),
    I2_recolonizations = sum(I2_recolonization),
    
    .groups = "drop"
  )

## ------------------------------------------------------------
## Overall summary
## ------------------------------------------------------------

strain2_intro_time <- strain2_delay * dt

overall_summary <- meta_summary |>
  summarise(
    strain1_globally_persistent =
      all(global_I1 > 0),
    
    strain2_persistent_after_intro =
      all(global_I2[step >= strain2_intro_time] > 0),
    
    min_I1_occupied_patches =
      min(occupied_I1),
    
    mean_I1_occupied_patches =
      mean(occupied_I1),
    
    min_I2_occupied_patches_after_intro =
      min(occupied_I2[step >= strain2_intro_time]),
    
    mean_I2_occupied_patches_after_intro =
      mean(occupied_I2[step >= strain2_intro_time]),
    
    total_I1_local_extinctions =
      sum(I1_local_extinctions),
    
    total_I1_recolonizations =
      sum(I1_recolonizations),
    
    total_I2_local_extinctions_after_intro =
      sum(
        I2_local_extinctions[
          step >= strain2_intro_time
        ]
      ),
    
    total_I2_recolonizations_after_intro =
      sum(
        I2_recolonizations[
          step >= strain2_intro_time
        ]
      )
  )

print(overall_summary)

## ------------------------------------------------------------
## Patch-level summary
## ------------------------------------------------------------

patch_summary <- patch_state |>
  group_by(patch) |>
  summarise(
    I1_extinction_events =
      sum(I1_local_extinction),
    
    I1_recolonization_events =
      sum(I1_recolonization),
    
    I2_extinction_events_after_intro =
      sum(
        I2_local_extinction[
          step >= strain2_intro_time
        ]
      ),
    
    I2_recolonization_events_after_intro =
      sum(
        I2_recolonization[
          step >= strain2_intro_time
        ]
      ),
    
    fraction_time_I1_present =
      mean(I1 > 0),
    
    fraction_time_I2_present_after_intro =
      mean(
        I2[step >= strain2_intro_time] > 0
      ),
    
    mean_population =
      mean(N),
    
    .groups = "drop"
  )

print(summary(patch_summary$I1_extinction_events))
print(summary(patch_summary$I1_recolonization_events))
print(summary(patch_summary$I2_extinction_events_after_intro))
print(summary(patch_summary$I2_recolonization_events_after_intro))

## ------------------------------------------------------------
## Plot 1: four-state patch x time raster
## ------------------------------------------------------------

p_raster <- ggplot(
  patch_state,
  aes(
    x = step,
    y = factor(patch),
    fill = infection_state
  )
) +
  geom_raster() +
  scale_fill_manual(
    values = c(
      "Neither" = "white",
      "Strain 1 only" = "blue",
      "Strain 2 only" = "red",
      "Both" = "black"
    ),
    drop = FALSE,
    name = NULL
  ) +
  geom_vline(
    xintercept = strain2_intro_time,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  labs(
    x = "Time",
    y = "Patch",
    title = sprintf(
      paste0(
        "Two-strain metapopulation dynamics: ",
        "R01 = %g, R02 = %g, K = %g, alpha = %g"
      ),
      R01,
      R02,
      K,
      alpha
    )
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top"
  )

print(p_raster)

## ------------------------------------------------------------
## Plot 2: occupied patches by strain
## ------------------------------------------------------------

occupancy_summary <- meta_summary |>
  select(
    step,
    occupied_I1,
    occupied_I2
  ) |>
  pivot_longer(
    cols = c(
      occupied_I1,
      occupied_I2
    ),
    names_to = "strain",
    values_to = "occupied_patches"
  ) |>
  mutate(
    strain = recode(
      strain,
      occupied_I1 = "Strain 1",
      occupied_I2 = "Strain 2"
    )
  )

p_occupancy <- ggplot(
  occupancy_summary,
  aes(
    x = step,
    y = occupied_patches,
    colour = strain
  )
) +
  geom_line() +
  geom_vline(
    xintercept = strain2_intro_time,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  labs(
    x = "Time",
    y = "Number of infected patches",
    colour = NULL,
    title = "Patch occupancy by strain"
  )

print(p_occupancy)

## ------------------------------------------------------------
## Plot 3: global infected population by strain
## ------------------------------------------------------------

global_summary <- meta_summary |>
  select(
    step,
    global_I1,
    global_I2
  ) |>
  pivot_longer(
    cols = c(
      global_I1,
      global_I2
    ),
    names_to = "strain",
    values_to = "global_I"
  ) |>
  mutate(
    strain = recode(
      strain,
      global_I1 = "Strain 1",
      global_I2 = "Strain 2"
    )
  )

p_global <- ggplot(
  global_summary,
  aes(
    x = step,
    y = global_I,
    colour = strain
  )
) +
  geom_line() +
  geom_vline(
    xintercept = strain2_intro_time,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_y_log10() +
  labs(
    x = "Time",
    y = "Global infected population",
    colour = NULL,
    title = "Global prevalence by strain"
  )

print(p_global)

## ------------------------------------------------------------
## Plot 4: extinction and recolonization events
## ------------------------------------------------------------

event_summary <- meta_summary |>
  select(
    step,
    I1_local_extinctions,
    I1_recolonizations,
    I2_local_extinctions,
    I2_recolonizations
  ) |>
  pivot_longer(
    cols = -step,
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
  geom_vline(
    xintercept = strain2_intro_time,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  labs(
    x = "Time",
    y = "Number of patch transitions",
    colour = NULL,
    title = "Local extinction and recolonization"
  )

print(p_events)

## ------------------------------------------------------------
## Save outputs
## ------------------------------------------------------------

outdir <- here::here(
  "odin",
  "outputs",
  "twostrain_rescue"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  file.path(
    outdir,
    "four_state_occupancy_raster.pdf"
  ),
  p_raster,
  width = 10,
  height = 7
)

ggsave(
  file.path(
    outdir,
    "occupied_patches_by_strain.pdf"
  ),
  p_occupancy,
  width = 8,
  height = 5
)

ggsave(
  file.path(
    outdir,
    "global_prevalence_by_strain.pdf"
  ),
  p_global,
  width = 8,
  height = 5
)

ggsave(
  file.path(
    outdir,
    "extinction_recolonization_by_strain.pdf"
  ),
  p_events,
  width = 9,
  height = 5
)

saveRDS(
  runs,
  file.path(
    outdir,
    "trajectory.rds"
  )
)

write.csv(
  patch_state,
  file.path(
    outdir,
    "patch_state.csv"
  ),
  row.names = FALSE
)

write.csv(
  meta_summary,
  file.path(
    outdir,
    "metapop_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  patch_summary,
  file.path(
    outdir,
    "patch_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  overall_summary,
  file.path(
    outdir,
    "overall_summary.csv"
  ),
  row.names = FALSE
)
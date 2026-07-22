library(plagueMetapop)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(odin)

theme_set(theme_bw())


## ------------------------------------------------------------
## Parameters
## ------------------------------------------------------------

dt <- 1

gamma <- 0.2

R0 <- 2.5

beta0 <- R0 * gamma

r <- 0.02

alpha <- 1e-4

K <- 3e3

n_patch <- 200


## Seasonal forcing

season_period <- 365

seasonal_amp <- 0.40

## Day 0 = January 1
## Maximum transmission around January 15

peak_day <- 15


## Run for ten years

t_max <- 10 * 365

nt <- round(t_max / dt)

seed <- 101


## ------------------------------------------------------------
## Mean-transmission deterministic equilibrium
## ------------------------------------------------------------

S_star <- gamma * K / beta0

I_star <- r * (K - S_star) / beta0


cat(
  "Mean-transmission deterministic equilibrium:\n"
)

cat(
  "S* =", S_star, "\n"
)

cat(
  "I* =", I_star, "\n"
)


cat(
  "\nSeasonal transmission:\n"
)

cat(
  "Mean beta =", beta0, "\n"
)

cat(
  "Minimum beta =",
  beta0 * (1 - seasonal_amp),
  "\n"
)

cat(
  "Maximum beta =",
  beta0 * (1 + seasonal_amp),
  "\n"
)

cat(
  "Minimum R0 =",
  R0 * (1 - seasonal_amp),
  "\n"
)

cat(
  "Maximum R0 =",
  R0 * (1 + seasonal_amp),
  "\n"
)


## ------------------------------------------------------------
## Initial conditions
## ------------------------------------------------------------

set.seed(seed)


S_ini <- rpois(
  n_patch,
  lambda = S_star
)


I_ini <- cbind(
  rpois(
    n_patch,
    lambda = I_star
  ),
  rep(
    0L,
    n_patch
  )
)


I2_ini <- rep(
  0L,
  n_patch
)


## ------------------------------------------------------------
## Compile seasonal odin model
## ------------------------------------------------------------

model_file <- here::here(
  "fadeout",
  "seasonal",
  "seasonal_model_metapop.R"
)


if (!file.exists(model_file)) {
  stop(
    "Odin model file does not exist: ",
    model_file
  )
}


gen <- suppressMessages(
  odin::odin(model_file)
)


## ------------------------------------------------------------
## Construct model
## ------------------------------------------------------------

mod <- gen$new(
  beta = c(
    beta0,
    0
  ),
  
  gamma = c(
    gamma,
    gamma
  ),
  
  dt = dt,
  
  I_ini = I_ini,
  
  S_ini = S_ini,
  
  I2_ini = I2_ini,
  
  alpha = alpha,
  
  strain2_delay =
    .Machine$integer.max,
  
  r = rep(
    r,
    n_patch
  ),
  
  K = rep(
    K,
    n_patch
  ),
  
  season_period =
    season_period,
  
  seasonal_amp =
    seasonal_amp,
  
  peak_day =
    peak_day,
  
  n_patch =
    n_patch
)


## ------------------------------------------------------------
## Run stochastic simulation
## ------------------------------------------------------------

raw <- mod$run(
  seq(
    0L,
    nt
  )
)


if (dt != 1) {
  raw[, "step"] <-
    raw[, "step"] * dt
}


runs <- conv_odin(raw)


## ------------------------------------------------------------
## Extract patch states
## ------------------------------------------------------------

patch_state <- runs |>
  filter(
    state %in%
      c(
        "S",
        "I1"
      )
  ) |>
  select(
    step,
    patch,
    state,
    value
  ) |>
  pivot_wider(
    names_from = state,
    values_from = value
  ) |>
  rename(
    I = I1
  ) |>
  arrange(
    patch,
    step
  ) |>
  mutate(
    N = S + I,
    
    occupied =
      as.integer(I > 0)
  )


## ------------------------------------------------------------
## Local extinction and recolonization
## ------------------------------------------------------------

patch_state <- patch_state |>
  group_by(patch) |>
  arrange(
    step,
    .by_group = TRUE
  ) |>
  mutate(
    prev_occupied =
      lag(occupied),
    
    local_extinction =
      as.integer(
        prev_occupied == 1 &
          occupied == 0
      ),
    
    recolonization =
      as.integer(
        prev_occupied == 0 &
          occupied == 1
      )
  ) |>
  ungroup() |>
  mutate(
    local_extinction =
      replace_na(
        local_extinction,
        0L
      ),
    
    recolonization =
      replace_na(
        recolonization,
        0L
      )
  )


## ------------------------------------------------------------
## Seasonal beta and R0
## ------------------------------------------------------------

seasonal_curve <- tibble(
  step = seq(
    0,
    t_max,
    by = dt
  )
) |>
  mutate(
    beta_t =
      beta0 * (
        1 +
          seasonal_amp *
          cos(
            2 * pi *
              (step - peak_day) /
              season_period
          )
      ),
    
    R0_t =
      beta_t / gamma,
    
    year =
      floor(step / 365) + 1,
    
    day_of_year =
      step %% 365
  )


## ------------------------------------------------------------
## Metapopulation summary
## ------------------------------------------------------------

meta_summary <- patch_state |>
  group_by(step) |>
  summarise(
    occupied_patches =
      sum(occupied),
    
    global_I =
      sum(I),
    
    mean_I =
      mean(I),
    
    mean_N =
      mean(N),
    
    local_extinctions =
      sum(local_extinction),
    
    recolonizations =
      sum(recolonization),
    
    .groups = "drop"
  ) |>
  left_join(
    seasonal_curve,
    by = "step"
  )


## ------------------------------------------------------------
## Overall summary
## ------------------------------------------------------------

burnin <- 365


analysis_summary <- meta_summary |>
  filter(
    step >= burnin
  )


overall_summary <- analysis_summary |>
  summarise(
    globally_persistent =
      all(global_I > 0),
    
    min_occupied_patches =
      min(occupied_patches),
    
    mean_occupied_patches =
      mean(occupied_patches),
    
    max_occupied_patches =
      max(occupied_patches),
    
    total_local_extinctions =
      sum(local_extinctions),
    
    total_recolonizations =
      sum(recolonizations),
    
    mean_patch_population =
      mean(mean_N)
  )


print(overall_summary)


## ------------------------------------------------------------
## Patch-level summary
## ------------------------------------------------------------

patch_summary <- patch_state |>
  filter(
    step >= burnin
  ) |>
  group_by(patch) |>
  summarise(
    extinction_events =
      sum(local_extinction),
    
    recolonization_events =
      sum(recolonization),
    
    fraction_time_infected =
      mean(occupied),
    
    mean_population =
      mean(N),
    
    .groups = "drop"
  )


print(
  summary(
    patch_summary$extinction_events
  )
)


print(
  summary(
    patch_summary$recolonization_events
  )
)


## ------------------------------------------------------------
## Plot 1: occupancy raster
## ------------------------------------------------------------

p_raster <- ggplot(
  patch_state,
  aes(
    x = step,
    y = factor(patch),
    fill = factor(occupied)
  )
) +
  geom_raster() +
  scale_fill_manual(
    values = c(
      "0" = "white",
      "1" = "black"
    ),
    labels = c(
      "uninfected",
      "infected"
    ),
    name = NULL
  ) +
  labs(
    x = "Time (days)",
    y = "Patch",
    
    title = sprintf(
      paste0(
        "Seasonal single-strain dynamics: ",
        "mean R0 = %.1f, amplitude = %.2f, ",
        "K = %g, alpha = %g"
      ),
      R0,
      seasonal_amp,
      K,
      alpha
    )
  ) +
  theme(
    axis.text.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank(),
    
    legend.position =
      "top"
  )


print(p_raster)


## ------------------------------------------------------------
## Plot 2: infected patch occupancy
## ------------------------------------------------------------

p_occupancy <- ggplot(
  meta_summary,
  aes(
    x = step,
    y = occupied_patches
  )
) +
  geom_line() +
  labs(
    x = "Time (days)",
    y = "Number of infected patches",
    title = "Metapopulation occupancy"
  )


print(p_occupancy)


## ------------------------------------------------------------
## Plot 3: occupancy and seasonal R0
## ------------------------------------------------------------

R0_scale <-
  n_patch /
  max(
    meta_summary$R0_t
  )


p_occupancy_season <- ggplot(
  meta_summary,
  aes(
    x = step
  )
) +
  geom_line(
    aes(
      y = occupied_patches
    )
  ) +
  geom_line(
    aes(
      y = R0_t * R0_scale
    ),
    linetype = "dashed"
  ) +
  scale_y_continuous(
    name =
      "Number of infected patches",
    
    sec.axis = sec_axis(
      ~ . / R0_scale,
      name = "Seasonal R0"
    )
  ) +
  labs(
    x = "Time (days)",
    title =
      "Patch occupancy and seasonal transmission"
  )


print(p_occupancy_season)


## ------------------------------------------------------------
## Plot 4: extinction and recolonization events
## ------------------------------------------------------------

event_summary <- meta_summary |>
  mutate(
    time_bin =
      floor(step / 30) * 30
  ) |>
  group_by(time_bin) |>
  summarise(
    local_extinctions =
      sum(local_extinctions),
    
    recolonizations =
      sum(recolonizations),
    
    .groups = "drop"
  ) |>
  pivot_longer(
    cols = c(
      local_extinctions,
      recolonizations
    ),
    
    names_to = "event",
    
    values_to = "count"
  )


p_events <- ggplot(
  event_summary,
  aes(
    x = time_bin,
    y = count,
    colour = event
  )
) +
  geom_line() +
  labs(
    x = "Time (days)",
    y = "Events per 30 days",
    colour = NULL,
    title =
      "Local extinction and recolonization"
  )


print(p_events)


## ------------------------------------------------------------
## Plot 5: mean patch population
## ------------------------------------------------------------

p_population <- ggplot(
  meta_summary,
  aes(
    x = step,
    y = mean_N
  )
) +
  geom_line() +
  labs(
    x = "Time (days)",
    y = "Mean patch population",
    title =
      "Mean host population per patch"
  )


print(p_population)


## ------------------------------------------------------------
## Save outputs
## ------------------------------------------------------------

outdir <- here::here(
  "fadeout",
  "output",
  "single_strain"
)


dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


ggsave(
  file.path(
    outdir,
    "occupancy_raster.pdf"
  ),
  p_raster,
  width = 10,
  height = 7
)


ggsave(
  file.path(
    outdir,
    "occupied_patches.pdf"
  ),
  p_occupancy,
  width = 8,
  height = 5
)


ggsave(
  file.path(
    outdir,
    "occupancy_seasonal_R0.pdf"
  ),
  p_occupancy_season,
  width = 9,
  height = 5
)


ggsave(
  file.path(
    outdir,
    "extinction_recolonization.pdf"
  ),
  p_events,
  width = 9,
  height = 5
)


ggsave(
  file.path(
    outdir,
    "mean_patch_population.pdf"
  ),
  p_population,
  width = 8,
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

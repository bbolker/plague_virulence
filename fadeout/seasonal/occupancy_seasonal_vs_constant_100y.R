## Seasonal vs constant-transmission patch occupancy over 100 years.
##
## Standalone companion to seasonal_single_strain_metapop.R, which runs the
## same model for ten years and writes to fadeout/output/single_strain.
## Nothing here overwrites those outputs.
##
## Two single realizations of the stochastic single-strain metapopulation,
## identical in every parameter except the seasonal amplitude:
##
##   seasonal   seasonal_amp = 0.40
##   constant   seasonal_amp = 0
##
## beta_eff[j](t) = beta[j] * (1 + seasonal_amp * cos(2*pi*(t - peak_day)/P)),
## so the cosine integrates to zero over a whole year and both runs share the
## same annual mean transmission, hence the same mean R0. The comparison is
## therefore seasonal forcing at fixed mean transmission, not a change in mean.
##
## The state is thinned to every THIN-th day before conv_odin(). One day is far
## below one output pixel at these axis lengths, and conv_odin() on the full
## 100-year state is the single slowest step in the pipeline (~32 s vs ~4 s).
##
## NOTE: thinning makes the one-step local-extinction and recolonization
## counts meaningless, because any episode shorter than THIN days disappears.
## Those event counts are deliberately NOT computed here; see
## seasonal_single_strain_metapop.R for the full-resolution event diagnostics.

library(plagueMetapop)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(odin)

theme_set(theme_bw())


## ------------------------------------------------------------
## Parameters (matching seasonal_single_strain_metapop.R)
## ------------------------------------------------------------

dt <- 1

gamma <- 0.2

R0 <- 2.5

beta0 <- R0 * gamma

r <- 0.02

alpha <- 1e-4

K <- 3e3

n_patch <- 200

season_period <- 365

peak_day <- 15

seed <- 101


## Run for one hundred years

n_years <- 100

t_max <- n_years * 365

nt <- round(t_max / dt)


## Keep every THIN-th day of output

THIN <- 10


## Amplitudes to compare; names become the scenario factor levels

amps <- c(
  seasonal = 0.40,
  constant = 0.00
)


## ------------------------------------------------------------
## Deterministic equilibrium at the mean transmission rate
## ------------------------------------------------------------

S_star <- gamma * K / beta0

I_star <- r * (K - S_star) / beta0


## ------------------------------------------------------------
## Compile the shared odin model
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


## odin::odin() substitutes its argument, so the path must be a variable.
gen <- suppressMessages(
  odin::odin(model_file)
)


## ------------------------------------------------------------
## One realization at a given seasonal amplitude
## ------------------------------------------------------------

## set.seed() is called inside this function so that every scenario starts
## from the same RNG state and the same initial conditions. Drawing the
## initial conditions once outside and then running several simulations in
## one session would advance the stream between runs and silently compare
## different realizations.

run_scenario <- function(seasonal_amp) {

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

  ## Thin before conv_odin(), which is where the time goes.
  keep <- seq(
    1,
    nrow(raw),
    by = THIN
  )

  runs <- conv_odin(raw[keep, , drop = FALSE])

  runs |>
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
    mutate(
      N = S + I,

      occupied =
        as.integer(I > 0)
    ) |>
    group_by(step) |>
    summarise(
      occupied_patches =
        sum(occupied),

      global_I =
        sum(I),

      mean_N =
        mean(N),

      .groups = "drop"
    ) |>
    mutate(
      years = step / 365
    )
}


## ------------------------------------------------------------
## Run both scenarios
## ------------------------------------------------------------

t_start <- Sys.time()

meta_summary <- lapply(
  names(amps),
  function(nm) {
    cat(
      "running scenario: ",
      nm,
      " (seasonal_amp = ",
      amps[[nm]],
      ")\n",
      sep = ""
    )

    run_scenario(amps[[nm]]) |>
      mutate(
        scenario = nm,
        seasonal_amp = amps[[nm]]
      )
  }
) |>
  bind_rows() |>
  mutate(
    scenario = factor(
      scenario,
      levels = names(amps)
    )
  )

cat(
  "elapsed: ",
  round(
    as.numeric(
      difftime(
        Sys.time(),
        t_start,
        units = "secs"
      )
    ),
    1
  ),
  " s\n",
  sep = ""
)


## ------------------------------------------------------------
## Seasonal R0 curve, for the zoomed panel only
## ------------------------------------------------------------

seasonal_curve <- expand_grid(
  step = seq(
    0,
    t_max,
    by = THIN
  ),
  scenario = factor(
    names(amps),
    levels = names(amps)
  )
) |>
  mutate(
    seasonal_amp =
      amps[as.character(scenario)],

    R0_t =
      R0 * (
        1 +
          seasonal_amp *
          cos(
            2 * pi *
              (step - peak_day) /
              season_period
          )
      ),

    years = step / 365
  )


## ------------------------------------------------------------
## Summaries
## ------------------------------------------------------------

burnin <- 365


overall_summary <- meta_summary |>
  filter(
    step >= burnin
  ) |>
  group_by(scenario) |>
  summarise(
    seasonal_amp =
      first(seasonal_amp),

    globally_persistent =
      all(global_I > 0),

    extinction_year =
      if (any(global_I == 0)) {
        min(years[global_I == 0])
      } else {
        NA_real_
      },

    min_occupied_patches =
      min(occupied_patches),

    mean_occupied_patches =
      mean(occupied_patches),

    max_occupied_patches =
      max(occupied_patches),

    mean_patch_population =
      mean(mean_N),

    .groups = "drop"
  )

print(overall_summary)


## Mean occupancy per decade, as a coarse trend check
decade_summary <- meta_summary |>
  filter(
    step >= burnin
  ) |>
  mutate(
    decade =
      floor(years / 10) * 10
  ) |>
  group_by(
    scenario,
    decade
  ) |>
  summarise(
    mean_occupied_patches =
      mean(occupied_patches),

    .groups = "drop"
  )

print(
  decade_summary |>
    pivot_wider(
      names_from = scenario,
      values_from = mean_occupied_patches
    )
)


## ------------------------------------------------------------
## Outputs
## ------------------------------------------------------------

outdir <- here::here(
  "fadeout",
  "output",
  "seasonal_vs_constant_100y"
)


dir.create(
  file.path(
    outdir,
    "data"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  file.path(
    outdir,
    "figures"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  meta_summary,
  file.path(
    outdir,
    "data",
    "occupancy_thinned.csv"
  ),
  row.names = FALSE
)


write.csv(
  overall_summary,
  file.path(
    outdir,
    "data",
    "overall_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  decade_summary,
  file.path(
    outdir,
    "data",
    "decade_summary.csv"
  ),
  row.names = FALSE
)


## Colour encodes the variable, linetype encodes the scenario.
occupancy_colour <- "black"

R0_colour <- "red"

scenario_linetypes <- c(
  seasonal = "solid",
  constant = "dashed"
)


## ------------------------------------------------------------
## Direct labels instead of legends
## ------------------------------------------------------------

## Label each series at the right-hand end of the panel, at the mean of its
## last `window` years. Both trajectories wander over the same range, so the
## positions are pushed apart to a minimum gap to keep the text legible.

label_positions <- function(dat,
                            x_max,
                            window = 5,
                            min_gap = 12) {

  out <- dat |>
    filter(
      years >= x_max - window
    ) |>
    group_by(scenario) |>
    summarise(
      y = mean(occupied_patches),
      .groups = "drop"
    ) |>
    arrange(y) |>
    mutate(
      x = x_max,
      label = as.character(scenario)
    )

  ## Push overlapping labels apart, lowest first.
  if (nrow(out) > 1) {
    for (i in 2:nrow(out)) {
      if (out$y[i] - out$y[i - 1] < min_gap) {
        out$y[i] <- out$y[i - 1] + min_gap
      }
    }
  }

  out
}


## ------------------------------------------------------------
## Plot 1: 100-year occupancy comparison
## ------------------------------------------------------------

## The seasonal R0 curve is omitted here on purpose: 100 annual cycles across
## the panel width render as a solid band. It appears in the zoomed panel.

labs_100y <- label_positions(
  meta_summary,
  x_max = n_years
)


p_100y <- ggplot(
  meta_summary,
  aes(
    x = years,
    y = occupied_patches,
    linetype = scenario
  )
) +
  geom_line(
    colour = occupancy_colour,
    linewidth = 0.3
  ) +
  geom_text(
    data = labs_100y,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    colour = occupancy_colour,
    hjust = 0,
    nudge_x = 1,
    size = 3.5
  ) +
  scale_linetype_manual(
    values = scenario_linetypes,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      n_years,
      by = 10
    ),

    ## Room for the direct labels
    limits = c(
      0,
      n_years * 1.1
    )
  ) +
  labs(
    x = "Time (years)",
    y = "Number of infected patches",

    title = sprintf(
      paste0(
        "Patch occupancy over %d years: ",
        "seasonal vs constant transmission"
      ),
      n_years
    ),

    subtitle = sprintf(
      paste0(
        "mean R0 = %.1f, K = %g, alpha = %g, ",
        "%d patches, single realization (seed %d), ",
        "every %dth day plotted"
      ),
      R0,
      K,
      alpha,
      n_patch,
      seed,
      THIN
    )
  )


ggsave(
  file.path(
    outdir,
    "figures",
    "occupancy_seasonal_vs_constant_100y.pdf"
  ),
  p_100y,
  width = 11,
  height = 5
)


## ------------------------------------------------------------
## Plot 2: faceted version, one panel per scenario
## ------------------------------------------------------------

p_facet <- ggplot(
  meta_summary,
  aes(
    x = years,
    y = occupied_patches
  )
) +
  geom_line(
    colour = occupancy_colour,
    linewidth = 0.3
  ) +
  facet_wrap(
    ~ scenario,
    ncol = 1
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      n_years,
      by = 10
    )
  ) +
  labs(
    x = "Time (years)",
    y = "Number of infected patches",

    title = sprintf(
      "Patch occupancy over %d years, by scenario",
      n_years
    )
  )


ggsave(
  file.path(
    outdir,
    "figures",
    "occupancy_seasonal_vs_constant_100y_facet.pdf"
  ),
  p_facet,
  width = 11,
  height = 7
)


## ------------------------------------------------------------
## Plot 3: first ten years, with the seasonal R0 curve
## ------------------------------------------------------------

zoom_years <- 10


zoom_occ <- meta_summary |>
  filter(
    years <= zoom_years
  )


zoom_R0 <- seasonal_curve |>
  filter(
    years <= zoom_years
  )


## Secondary axis for R0, as in seasonal_single_strain_metapop.R.
R0_scale <- n_patch /
  max(seasonal_curve$R0_t)


labs_zoom <- label_positions(
  zoom_occ,
  x_max = zoom_years,
  window = 1
)


## One red label for the R0 curves; linetype already distinguishes the
## scenarios, so the label is placed on the constant (flat) curve.
lab_R0 <- tibble(
  x = zoom_years,
  y = R0 * R0_scale,
  label = "R0"
)


p_zoom <- ggplot(
  zoom_occ,
  aes(
    x = years,
    linetype = scenario
  )
) +
  geom_line(
    aes(
      y = occupied_patches
    ),
    colour = occupancy_colour,
    linewidth = 0.3
  ) +
  geom_line(
    data = zoom_R0,
    aes(
      y = R0_t * R0_scale
    ),
    colour = R0_colour,
    linewidth = 0.4
  ) +
  geom_text(
    data = labs_zoom,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    colour = occupancy_colour,
    hjust = 0,
    nudge_x = 0.12,
    size = 3.5
  ) +
  geom_text(
    data = lab_R0,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    colour = R0_colour,
    hjust = 0,
    nudge_x = 0.12,
    size = 3.5
  ) +
  scale_linetype_manual(
    values = scenario_linetypes,
    guide = "none"
  ) +
  scale_y_continuous(
    name =
      "Number of infected patches",

    sec.axis = sec_axis(
      ~ . / R0_scale,
      name = "Seasonal R0"
    )
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      zoom_years,
      by = 1
    ),

    ## Room for the direct labels
    limits = c(
      0,
      zoom_years * 1.1
    )
  ) +
  labs(
    x = "Time (years)",

    title = sprintf(
      paste0(
        "First %d years: occupancy (black) and R0 (red); ",
        "solid = seasonal, dashed = constant"
      ),
      zoom_years
    )
  ) +
  theme(
    axis.title.y.right =
      element_text(colour = R0_colour),

    axis.text.y.right =
      element_text(colour = R0_colour)
  )


ggsave(
  file.path(
    outdir,
    "figures",
    "occupancy_seasonal_vs_constant_first10y.pdf"
  ),
  p_zoom,
  width = 10,
  height = 5
)


cat(
  "Output: ",
  outdir,
  "\n",
  sep = ""
)

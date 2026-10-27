## Shared settings and functions for the pairwise invasion runs
## (pairwise_invasion_100y.R locally, and
## odin/sharcnet/pairwise_invasion_factorial_run_array.R on Nibi).
##
## A run is a 100-year resident phase from the ratchet initial conditions of
## occupancy_factorial_100y_batch.R, then 10 invader infecteds in one patch
## chosen at random among those the resident occupies (any patch if the
## resident is extinct), then 50 years of two-strain dynamics. See
## pairwise_invasion_100y.R for the outcome definitions.
##
## The caller loads dplyr and compiles the model with make_generator().


## ------------------------------------------------------------
## Fixed settings
## ------------------------------------------------------------

dt <- 1                   ## days
gamma <- 0.2              ## per day
r <- 0.125 * gamma        ## 0.125 per infectious period
n_patch <- 200
season_period <- 365
peak_day <- 15

## Infected count of a freshly invaded patch in the ratchet initial conditions
I_outbreak <- 10

## With susc_eq = TRUE, initial S is this range of multiples of S* + I*,
## capped at K
susc_eq_range <- c(0.7, 1.3)

## Invader infecteds put into the seed patch
I_invader <- 10

resident_years <- 100
invasion_years <- 50
nt_res <- round(resident_years * 365 / dt)
nt_inv <- round(invasion_years * 365 / dt)

## Keep every THIN-th day of the saved occupancy series, from save_from_year on
THIN <- 10
save_from_year <- 80

## Occupancy is summarized over the last summary_years years of each phase
summary_years <- 10

I1_cols <- sprintf("I[%d,1]", seq_len(n_patch))
I2_cols <- sprintf("I[%d,2]", seq_len(n_patch))
S_cols <- sprintf("S[%d]", seq_len(n_patch))


## ------------------------------------------------------------
## Model
## ------------------------------------------------------------

## odin::odin() substitutes its argument, so the path must be a variable.
## The compiled model reaches forked workers, hence plan(multicore).
make_generator <- function(model_file) {
  if (!file.exists(model_file)) {
    stop("Odin model file does not exist: ", model_file)
  }
  suppressMessages(odin::odin(model_file))
}

new_model <- function(gen, beta, K, alpha, seasonal_amp,
                      S_ini, I1_ini, I2_ini) {
  gen$new(
    beta = beta,
    gamma = c(gamma, gamma),
    dt = dt,
    I_ini = cbind(I1_ini, rep(0, n_patch)),
    S_ini = S_ini,
    I2_ini = I2_ini,
    alpha = alpha,
    strain2_delay = .Machine$integer.max,
    r = rep(r, n_patch),
    K = rep(K, n_patch),
    season_period = season_period,
    seasonal_amp = seasonal_amp,
    peak_day = peak_day,
    n_patch = n_patch
  )
}

## first day with zero global I (NA if never), counting from row 1 = day 0
first_zero <- function(x, day) {
  i <- match(0, x)
  if (is.na(i)) NA_integer_ else as.integer(day[i])
}


## ------------------------------------------------------------
## One realization
## ------------------------------------------------------------

## The two phases are separate model instances: the seed patch depends on the
## resident's state at year 100. 100 years is a whole number of seasons, so
## restarting the model clock at 0 keeps the seasonal phase.
run_one <- function(run_id, R0_res, R0_inv, K, alpha, seasonal_amp, gen,
                    susc_eq = FALSE, ...) {

  t0 <- Sys.time()
  beta <- c(R0_res, R0_inv) * gamma
  mk <- function(S_ini, I1_ini, I2_ini) {
    new_model(gen, beta, K, alpha, seasonal_amp, S_ini, I1_ini, I2_ini)
  }

  ## Deterministic resident equilibrium at the mean transmission rate
  S_star <- gamma * K / beta[1]
  I_star <- r * (K - S_star) / beta[1]

  ## Ratchet initial conditions: u = 1 is a fresh outbreak in a fully
  ## susceptible patch, u = 0 is the endemic equilibrium
  u <- runif(n_patch, min = 0, max = 1)
  S_ini <- round(K * u + S_star * (1 - u))
  I_ini <- round(I_outbreak * u + I_star * (1 - u))
  ## Drawn after u, so the default runs' random streams are unchanged
  if (susc_eq) {
    S_ini <- round(pmin(K, (S_star + I_star) *
                          runif(n_patch, susc_eq_range[1], susc_eq_range[2])))
  }

  ## Resident phase (odin steps daily; only every THIN-th day is returned)
  steps_res <- seq(0L, nt_res, by = THIN)
  raw <- mk(S_ini, I_ini, rep(0, n_patch))$run(steps_res)
  I1_res <- raw[, I1_cols]
  S_end <- raw[nrow(raw), S_cols]
  I1_end <- I1_res[nrow(I1_res), ]
  rm(raw)

  day_res <- steps_res * dt
  occ1_res <- rowSums(I1_res > 0)
  glob1_res <- rowSums(I1_res)
  late_res <- day_res >= (resident_years - summary_years) * 365

  ## Seed patch: random among patches the resident occupies
  occupied <- which(I1_end > 0)
  if (!length(occupied)) occupied <- seq_len(n_patch)
  seed_patch <- occupied[sample.int(length(occupied), 1)]
  I2_ini <- replace(rep(0, n_patch), seed_patch, I_invader)

  ## Invasion phase, daily; day counts from the invasion
  raw <- mk(S_end, I1_end, I2_ini)$run(seq(0L, nt_inv))
  I1 <- raw[, I1_cols]
  I2 <- raw[, I2_cols]
  rm(raw)

  day <- seq(0L, nt_inv) * dt
  occ1 <- rowSums(I1 > 0)
  occ2 <- rowSums(I2 > 0)
  occ_both <- rowSums(I1 > 0 & I2 > 0)
  glob1 <- rowSums(I1)
  glob2 <- rowSums(I2)
  late <- day >= (invasion_years - summary_years) * 365
  n <- length(day)

  resident_absent <- glob1[1] == 0
  res_end <- glob1[n] > 0
  inv_end <- glob2[n] > 0

  summary <- tibble(
    run_id = run_id,
    seed_patch = seed_patch,
    seed_patch_I1 = I1_end[seed_patch],
    resident_extinct_pre = resident_absent,
    resident_extinction_day_pre = first_zero(glob1_res, day_res),
    occ1_pre = mean(occ1_res[late_res]) / n_patch,
    occ1_invasion = occ1[1] / n_patch,
    resident_extinction_day = first_zero(glob1, day),
    invader_extinction_day = first_zero(glob2, day),
    invader_peak_occ = max(occ2) / n_patch,
    occ1_final = occ1[n] / n_patch,
    occ2_final = occ2[n] / n_patch,
    occ_both_final = occ_both[n] / n_patch,
    occ1_late = mean(occ1[late]) / n_patch,
    occ2_late = mean(occ2[late]) / n_patch,
    occ_both_late = mean(occ_both[late]) / n_patch,
    outcome = case_when(
      resident_absent ~ "resident_absent",
      res_end & inv_end ~ "coexist",
      res_end ~ "invader_lost",
      inv_end ~ "resident_replaced",
      TRUE ~ "both_extinct"
    ),
    elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )

  ## Saved series: resident phase from save_from_year, then the invasion
  ## phase; days count from the start of the resident phase
  keep_res <- day_res >= save_from_year * 365 & day_res < nt_res * dt
  keep <- seq(1, n, by = THIN)

  list(
    summary = summary,
    thinned = tibble(
      run_id = run_id,
      day = as.integer(c(day_res[keep_res], nt_res * dt + day[keep])),
      occ1 = as.integer(c(occ1_res[keep_res], occ1[keep])),
      occ2 = as.integer(c(rep(0L, sum(keep_res)), occ2[keep])),
      occ_both = as.integer(c(rep(0L, sum(keep_res)), occ_both[keep])),
      global_I1 = as.integer(c(glob1_res[keep_res], glob1[keep])),
      global_I2 = as.integer(c(rep(0L, sum(keep_res)), glob2[keep]))
    )
  )
}


## ------------------------------------------------------------
## Summaries across replicates
## ------------------------------------------------------------

## One row per group (e.g. pair, or environment x pair). Outcome
## probabilities exclude runs where the resident was already extinct at the
## invasion.
summarise_pairs <- function(run_summary, ...) {
  run_summary |>
    group_by(...) |>
    summarise(
      n_runs = n(),
      n_resident_absent = sum(outcome == "resident_absent"),
      n_valid = sum(outcome != "resident_absent"),
      p_invader_lost = sum(outcome == "invader_lost") / n_valid,
      p_replace = sum(outcome == "resident_replaced") / n_valid,
      p_coexist = sum(outcome == "coexist") / n_valid,
      p_both_extinct = sum(outcome == "both_extinct") / n_valid,
      p_invader_persist = p_replace + p_coexist,
      median_invader_extinction_years = median(
        invader_extinction_day[outcome != "resident_absent"] / 365,
        na.rm = TRUE
      ),
      occ1_pre = mean(occ1_pre),
      occ1_final = mean(occ1_final[outcome != "resident_absent"]),
      occ2_final = mean(occ2_final[outcome != "resident_absent"]),
      occ1_late = mean(occ1_late[outcome != "resident_absent"]),
      occ2_late = mean(occ2_late[outcome != "resident_absent"]),
      .groups = "drop"
    )
}

## Parallel analytical comparison in which transient and persistent source
## patches have distinct infected-host contributions.  This script reads only
## existing metapopulation simulations and never overwrites earlier analyses.

library(deSolve)
library(dplyr)
library(ggplot2)
library(here)
library(mgcv)
library(plagueMetapop)
library(tidyr)

source(here::here("fadeout", "two_state_functions.R"))
source(here::here("fadeout", "occupancy_functions.R"))
theme_set(theme_bw())

hybrid_duration <- "--hybrid-duration" %in% commandArgs(trailingOnly = TRUE)
c_transient <- 0.5
probability_tolerance <- 1e-6
single_patch_file <- here::here("odin", "sharcnet", "outputs",
  "euler_onepatch_onestrain_extinct_logistic_continuous_demoggrid.rds")
raw_file <- here::here("fadeout", "output", "stochastic_patch_occupancy", "data",
  "stochastic_patch_occupancy_results.rds")
established_file <- here::here("fadeout", "output",
  "stochastic_episode_occupancy", "data", "established_occupancy_results.rds")
transient_file <- here::here("fadeout", "output", "stochastic_episode_occupancy",
  "data", "transient_source_pressure_timeseries.csv")
outdir <- here::here(
  "fadeout", "output",
  if (hybrid_duration) "two_state_occupancy_hybrid_TT" else
    "two_state_occupancy"
)
datadir <- file.path(outdir, "data")
figdir <- file.path(outdir, "figures")
for (f in c(single_patch_file, raw_file, established_file, transient_file)) {
  if (!file.exists(f)) stop("Required input file not found: ", f)
}
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

single_patch <- readRDS(single_patch_file)
raw_results <- readRDS(raw_file)
established_results <- readRDS(established_file)
transient_series <- read.csv(transient_file)
required_transient <- c("result_id", "time", "transient_patches",
  "persistent_patches", "infected_patches")
if (any(!c("R0", "K", "r", "ext_prob.I1") %in% names(single_patch)))
  stop("Single-patch data lack R0, K, r, or ext_prob.I1")
if (any(!required_transient %in% names(transient_series)))
  stop("Transient analysis lacks required episode-classified columns")
if (length(raw_results) != length(established_results))
  stop("Raw and established result counts differ")

estimate_P1 <- make_P1_demoggrid_estimator(single_patch)

metric <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  if (!any(ok)) return(c(RMSE = NA, max_abs = NA, correlation = NA))
  x <- observed[ok]; y <- predicted[ok]
  c(RMSE = sqrt(mean((x-y)^2)), max_abs = max(abs(x-y)),
    correlation = if (sd(x) == 0 || sd(y) == 0) NA_real_ else cor(x,y))
}
format_parameter <- function(name, value) {
  if (name == "K") sprintf("%s = %s", name,
    format(value, big.mark = ",", scientific = FALSE)) else
    sprintf("%s = %g", name, value)
}
fixed_subtitle <- function(group, p) {
  fixed <- setdiff(c("R0", "K", "r", "alpha"), group)
  paste0(paste(vapply(fixed, function(x) format_parameter(x, p[[x]]),
    character(1)), collapse = "; "),
    if (hybrid_duration)
      sprintf("; gamma = %g; hybrid T_T; tau = 50\n", p$gamma) else
      sprintf("; gamma = %g; c = %g; tau = 50\n", p$gamma, c_transient),
    sprintf("n_patch = %d; dt = %g; t_max = %g; S(0)=K-%g; I(0)=%g",
      p$n_patch, p$dt, p$t_max, p$I_outbreak, p$I_outbreak))
}

analyse_one <- function(i) {
  raw <- raw_results[[i]]; est <- established_results[[i]]
  if (!identical(raw$varied_parameter, est$varied_parameter) ||
      !isTRUE(all.equal(raw$parameter_value, est$parameter_value)))
    stop("Raw/established mismatch for result ", i)
  params <- raw$params
  raw_curve <- raw$raw_meta_summary |>
    select(step, stochastic_total = occupancy_fraction)
  persistent_curve <- est$occupancy_summary |>
    filter(tau == 50, !is.na(occupancy_fraction_established)) |>
    select(step, stochastic_persistent = occupancy_fraction_established)
  shift_index <- which.min(raw_curve$stochastic_total)
  time_shift <- raw_curve$step[shift_index]
  tr <- transient_series |>
    filter(varied_parameter == raw$varied_parameter,
      abs(parameter_value - raw$parameter_value) < 1e-12,
      time >= time_shift) |>
    transmute(step = time,
      stochastic_transient = transient_patches / params$n_patch,
      episode_persistent = persistent_patches / params$n_patch)
  simulation <- inner_join(raw_curve, persistent_curve, by = "step") |>
    filter(step >= time_shift) |>
    inner_join(tr, by = "step") |>
    mutate(time_post_burnout = step - time_shift)
  if (!nrow(simulation)) stop("No aligned transient data for result ", i)

  P1_estimate <- estimate_P1(params$R0, params$K, params$r)
  P1_raw <- P1_estimate$P1_raw
  P1 <- pmin(pmax(P1_raw, probability_tolerance), 1-probability_tolerance)
  I_star <- unname(plagueMetapop::ode_eq(params$R0, params$gamma,
    params$K, params$r, 1)[["eq_I"]])
  transient <- if (hybrid_duration) {
    compute_transient_outbreak_summary_I1(
      params$R0, params$K, params$r, params$I_outbreak
    )
  } else {
    compute_transient_outbreak_summary(
      params$R0, params$K, params$r, params$gamma,
      params$I_outbreak, c_transient
    )
  }
  if (hybrid_duration && !transient$endpoint_found)
    stop("No hybrid transient endpoint for result ", i,
      ": ", transient$status)
  if (!hybrid_duration && !transient$oscillatory)
    stop("Non-oscillatory closure for result ", i,
      ": argument=", transient$exact_argument)

  ## q0 is observed, not fitted: episode-classified transient occupancy at the
  ## same raw-occupancy minimum used by the original comparison.  We retain
  ## p0=P1 for comparability with the analytical workflow.
  p0 <- P1
  q0_observed <- simulation$stochastic_transient[1]
  if (p0 + q0_observed > 1) q0_observed <- 1-p0
  new <- solve_two_state_transient_Ibar(simulation$time_post_burnout,
    p0, q0_observed, P1, params$alpha, I_star, transient$Ibar_T,
    transient$T_T)
  new_zero_q0 <- solve_two_state_transient_Ibar(
    simulation$time_post_burnout, p0, 0, P1, params$alpha, I_star,
    transient$Ibar_T, transient$T_T)
  lambda <- params$alpha * I_star * P1
  one_p <- 1/(1 + ((1-P1)/P1)*exp(-lambda*simulation$time_post_burnout))

  cmp <- bind_cols(simulation,
    new |> select(new_p=p, new_q=q, new_total=total)) |>
    mutate(one_state_p = one_p)
  m_one <- metric(cmp$stochastic_persistent, cmp$one_state_p)
  m_new_p <- metric(cmp$stochastic_persistent, cmp$new_p)
  m_new_total <- metric(cmp$stochastic_total, cmp$new_total)
  m_new_q <- metric(cmp$stochastic_transient, cmp$new_q)
  m_new_zero_p <- metric(cmp$stochastic_persistent, new_zero_q0$p)
  m_new_zero_total <- metric(cmp$stochastic_total, new_zero_q0$total)
  diag <- data.frame(
    result_id=i, varied_parameter=raw$varied_parameter,
    parameter_value=raw$parameter_value, R0=params$R0, K=params$K,
    r=params$r, alpha=params$alpha, gamma=params$gamma,
    I0=params$I_outbreak, P1=P1, P1_raw=P1_raw,
    P1_exact_grid=P1_estimate$P1_exact_grid,
    P1_source=P1_estimate$source,
    I_star=I_star, I_peak=transient$I_peak,
    time_of_I_peak=transient$time_of_I_peak,
    duration_method=if (hybrid_duration) "I1_or_first_trough" else
      "half_Tosc",
    endpoint_type=if (hybrid_duration) transient$transient_endpoint else
      "half_Tosc",
    endpoint_status=if (hybrid_duration) transient$status else
      "success_half_Tosc",
    I_first_trough=if (hybrid_duration) transient$I_first_trough else
      NA_real_,
    time_of_first_trough=if (hybrid_duration)
      transient$time_of_first_trough else NA_real_,
    T_osc=if (hybrid_duration) transient$T_osc_old else transient$T_osc,
    T_osc_approx=if (hybrid_duration) NA_real_ else
      transient$T_osc_approx,
    c=if (hybrid_duration) NA_real_ else c_transient,
    T_T=transient$T_T, integral_I=transient$integral_I,
    Ibar_T=transient$Ibar_T,
    Ibar_T_over_I_star=transient$Ibar_T/I_star,
    I_peak_over_I_star=transient$I_peak/I_star,
    oscillation_argument=if (hybrid_duration) NA_real_ else
      transient$exact_argument,
    exact_oscillation_valid=if (hybrid_duration)
      transient$endpoint_found else transient$oscillatory,
    turnover_mapping=if (hybrid_duration)
      "normalized logistic deterministic invasion; hybrid endpoint" else
      transient$mapping,
    p0=p0, q0=q0_observed, q0_source="observed episode-classified transient occupancy at raw minimum",
    time_shift=time_shift,
    stochastic_persistent_at_shift=cmp$stochastic_persistent[1],
    stochastic_transient_at_shift=cmp$stochastic_transient[1],
    stochastic_total_at_shift=cmp$stochastic_total[1],
    final_p=tail(cmp$new_p,1), final_q=tail(cmp$new_q,1),
    final_total=tail(cmp$new_total,1),
    RMSE_one_state_p=m_one[["RMSE"]],
    RMSE_new_two_state_p=m_new_p[["RMSE"]],
    RMSE_new_two_state_total=m_new_total[["RMSE"]],
    RMSE_new_two_state_q=m_new_q[["RMSE"]],
    RMSE_new_two_state_zero_q0_p=m_new_zero_p[["RMSE"]],
    RMSE_new_two_state_zero_q0_total=m_new_zero_total[["RMSE"]],
    max_abs_error_one_state_p=m_one[["max_abs"]],
    max_abs_error_new_two_state_p=m_new_p[["max_abs"]],
    correlation_new_two_state_p=m_new_p[["correlation"]]
  )
  curves <- cmp |>
    transmute(result_id=i, varied_parameter=raw$varied_parameter,
      parameter_value=raw$parameter_value,
      parameter_label=format_parameter(raw$varied_parameter,
        raw$parameter_value), step, time_post_burnout,
      stochastic_persistent, stochastic_transient, stochastic_total,
      one_state_p, new_p, new_q, new_total)
  list(curves=curves, diagnostics=diag, params=params,
    transient_trajectory=mutate(transient$trajectory, result_id=i))
}

analyses <- lapply(seq_along(raw_results), analyse_one)
curves <- bind_rows(lapply(analyses, `[[`, "curves"))
diagnostics <- bind_rows(lapply(analyses, `[[`, "diagnostics"))
transient_trajectories <- bind_rows(lapply(analyses, `[[`,
  "transient_trajectory"))
write.csv(curves, file.path(datadir, "two_state_occupancy_curves.csv"),
  row.names=FALSE)
write.csv(diagnostics,
  file.path(datadir, "two_state_occupancy_diagnostics.csv"),
  row.names=FALSE)
write.csv(transient_trajectories,
  file.path(datadir, "deterministic_transient_invasion_trajectories.csv"),
  row.names=FALSE)

palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
for (group_name in unique(diagnostics$varied_parameter)) {
  rows <- diagnostics |> filter(varied_parameter==group_name) |> arrange(result_id)
  lev <- vapply(rows$parameter_value,
    function(x) format_parameter(group_name,x), character(1))
  pd <- curves |> filter(varied_parameter==group_name) |>
    mutate(parameter_label=factor(parameter_label,levels=lev)) |>
    select(time_post_burnout,parameter_label,stochastic_persistent,
      stochastic_transient,stochastic_total,new_p,new_q,new_total) |>
    pivot_longer(-c(time_post_burnout,parameter_label),
      names_to="series",values_to="occupancy") |>
    mutate(panel=case_when(
      grepl("persistent|_p$",series) ~ "Persistent occupancy",
      grepl("transient|_q$",series) ~ "Transient occupancy",
      TRUE ~ "Total current occupancy"),
      model=if_else(grepl("^stochastic",series),
        "Stochastic trajectory", "Deterministic approximation"))
  p <- ggplot(pd,aes(time_post_burnout,occupancy,colour=parameter_label,
    linetype=model,group=interaction(parameter_label,model,panel)))+
    geom_line(linewidth=.55)+facet_wrap(~panel,ncol=1)+
    scale_colour_manual(values=palette[seq_along(lev)])+
    scale_linetype_manual(values=c("Stochastic trajectory"="solid",
      "Deterministic approximation"="22"))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2),
      labels=scales::percent_format(accuracy=1))+
    labs(x="Time since each simulated raw-occupancy minimum",
      y="Fraction of patches",colour=paste("Varied parameter:",group_name),
      linetype=NULL,
      title=paste(
        if (hybrid_duration)
          "Hybrid-duration two-state approximation: varying" else
          "Transient infected-load two-state approximation: varying",
        group_name
      ),
      subtitle=fixed_subtitle(group_name,analyses[[rows$result_id[1]]]$params),
      caption=paste0(
        "p(0)=P1; q(0)=observed episode-classified transient occupancy at ",
        "the aligned minimum. Ibar_T is deterministic, not fitted.",
        if (hybrid_duration)
          " T_T uses the I=1-or-first-trough hybrid endpoint." else ""
      ))+
    theme(legend.position="right",plot.subtitle=element_text(size=9),
      plot.caption=element_text(hjust=0,size=8),
      plot.margin=margin(8,8,8,160,unit="pt"))
  ggsave(file.path(figdir,sprintf("two_state_occupancy_compare_%s.pdf",
    group_name)),p,width=12,height=10)
}

cat("Scenarios analysed: ",nrow(diagnostics),"\n",sep="")
cat("Valid oscillatory closures: ",sum(diagnostics$exact_oscillation_valid),
  "/",nrow(diagnostics),"\n",sep="")
cat("Ibar_T/Istar range: ",paste(signif(range(diagnostics$Ibar_T_over_I_star),4),collapse=" to "),"\n",sep="")
cat("Median one-state/new two-state persistent RMSE: ",
  signif(median(diagnostics$RMSE_one_state_p),4)," / ",
  signif(median(diagnostics$RMSE_new_two_state_p),4),"\n",sep="")
cat("Median new two-state total RMSE: ",
  signif(median(diagnostics$RMSE_new_two_state_total),4),"\n",sep="")
cat("Output: ",outdir,"\n",sep="")

if (!requireNamespace("burnout", quietly = TRUE)) {
  stop("Please install the burnout package: remotes::install_github('davidearn/burnout')", call. = FALSE)
}
if (!requireNamespace("gsl", quietly = TRUE)) {
  stop("Please install the gsl package: install.packages('gsl')", call. = FALSE)
}

source("two_strain_burnout.R")

ans <- two_strain_wave_burnout(
  R01 = 2.2,
  R02 = 1.8,
  epsilon = 0.01,
  N = 1e6,
  k1 = 1,
  k2 = 1,
  dt = 0.01,
  t_max = 200,
  n_x=4000L
)

cat("\nTwo-strain burnout summary\n")
cat("boundary event: ", ans$boundary_event, "\n", sep = "")
cat("entry case:     ", ans$entry_case, "\n", sep = "")
cat("first strain:   ", ans$first_strain, "\n", sep = "")
cat("event time:     ", signif(ans$event_time, 8), "\n", sep = "")

cat("\nCommon-boundary one-infective extinction probabilities\n")
print(
  data.frame(
    strain = c("strain 1", "strain 2"),
    qC = unname(ans$qC_by_strain),
    one_minus_qC = 1 - unname(ans$qC_by_strain)
  ),
  row.names = FALSE,
  digits = 10
)

cat("\nWave burnout probabilities\n")
print(
  data.frame(
    strain = c("strain 1", "strain 2"),
    Q = unname(ans$strain_burnout),
    one_minus_Q = 1 - unname(ans$strain_burnout)
  ),
  row.names = FALSE,
  digits = 10
)

cat("\nJoint outcomes\n")
print(
  data.frame(
    outcome = names(ans$outcomes_by_strain),
    probability = unname(ans$outcomes_by_strain)
  ),
  row.names = FALSE,
  digits = 10
)

cat("\nEntry diagnostics\n")
diag <- data.frame(
  quantity = c("x_event", "xa", "xb", "G1_given_2", "G2_given_1"),
  value = c(
    unname(ans$event_state["x"]),
    if (!is.null(ans$xa)) ans$xa else NA_real_,
    if (!is.null(ans$xb)) ans$xb else NA_real_,
    if (!is.null(ans$G1_given_2)) ans$G1_given_2 else NA_real_,
    if (!is.null(ans$G2_given_1)) ans$G2_given_1 else NA_real_
  )
)
print(diag[!is.na(diag$value), ], row.names = FALSE, digits = 10)

saveRDS(ans, file = "output/two_strain_burnout_result.rds")
cat("\nSaved full result to two_strain_burnout_result.rds\n")

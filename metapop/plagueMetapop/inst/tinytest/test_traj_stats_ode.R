library(plagueMetapop)

## Single-strain, single-patch ODE run matching ode_example.R
dt   <- 0.01
run1 <- discrete_run(
  beta_vec     = c(4, 0),
  K            = 1,
  r            = 0.125,
  n_patch      = 1,
  nt           = round(50 / dt),
  alpha        = 0,
  I_init       = c(0.001, 0),
  I_ini_method = "fixed",
  gamma        = 1,
  dt           = dt,
  def_file     = "ode_odin_def.R",
  stop_cond    = NULL,
  nsim         = 1,
  platform     = "odin"
)

## step column is in disease-generation units (0 to nt * dt = 5)
expect_equal(max(run1$step), 50,
             info = "ODE step column is scaled to disease-generation units")

## output carries the metapop_run class
expect_true(inherits(run1, "metapop_run"),
            info = "discrete_run() output has class 'metapop_run'")

## eval(attr(run, "call")) reproduces the same trajectory
run2 <- eval(attr(run1, "call"))
t1   <- traj_stats_ode(run1)
t2   <- traj_stats_ode(run2)
expect_equal(t1, t2, tolerance = 1e-3,
             info = "eval(attr(run, 'call')) reproduces identical traj_stats_ode output")

## return value has the expected names
expect_equal(
  names(t1),
  c("eq_S", "eq_I", "t_enter.boundary", "t_Imin", "Imin", "t_Smin", "Smin",
    "t_eqS", "t_leave.boundary", "trough_area", "any_small"),
  info = "traj_stats_ode returns correctly named vector"
)

## all statistics are finite for a well-behaved endemic trajectory
expect_true(all(is.finite(t1)),
            info = "traj_stats_ode: all statistics are finite for endemic trajectory")

## temporal ordering is consistent
expect_true(unname(t1["t_enter.boundary"]) < unname(t1["t_Imin"]),
            info = "t_enter.boundary precedes t_Imin")
expect_true(unname(t1["t_Imin"]) < unname(t1["t_leave.boundary"]),
            info = "t_Imin precedes t_leave.boundary")

## snapshot: full numeric values
expect_equal_to_reference(t1, "traj_stats_ode_snapshot.rds",
                          info = "traj_stats_ode snapshot")

## any_small == 1 for parameters that drive I1 below .Machine$double.xmin
## (ode_example.R parameters: high beta, very slow host growth → deep prolonged trough)
dt_ex <- 0.01; r_ex <- 0.003
run_ex <- suppressWarnings(discrete_run(
  beta_vec     = c(7, 0),
  K            = 1,
  r            = r_ex,
  n_patch      = 1,
  nt           = round(10 / r_ex / dt_ex),
  alpha        = 0,
  I_init       = c(0.001, 0),
  I_ini_method = "fixed",
  dt           = dt_ex,
  def_file     = "ode_odin_def.R",
  logistic_growth = 1,
  stop_cond    = NULL,
  nsim         = 1,
  platform     = "odin"
)) |> dplyr::filter(state != "I2")
expect_equal(traj_stats_ode(run_ex)[["any_small"]], 1,
             info = "any_small == 1 when I1 falls below .Machine$double.xmin")

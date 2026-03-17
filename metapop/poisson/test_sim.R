# Load required simulation functions
source("simulation_funs.R")

# Define parameters for a test run
test_params <- list(
  n_patches = 100,
  n_years = 500,
  K = 1e6,
  r = 0.5,
  c0 = 0.2,
  nu = 5,
  rho = 3,
  alpha = 5e-6,
  D = 1,
  R0 = 2.5,
  early_stop = FALSE 
)

# Run multiple stochastic simulations
nsim_test <- 100
cat(sprintf("Running %d simulations...\n", nsim_test))
multi_res <- mult_sim_mp(nsim = nsim_test, params = test_params, ncores = 1, seed = 42)

# Calculate summary statistics (quasi-equilibrium based on the last 100 years)
summary_stats <- sumfun(multi_res, nsteps = 100)
qe_mean <- summary_stats["infected_patches"]
ext_rate <- summary_stats["extinction_rate"]
n_ext <- summary_stats["n_extinct"]

cat("\nSummary Statistics:\n")
print(summary_stats)

# Extract infection matrix and compute mean trajectory for surviving sims
inf_mat <- multi_res$infected_patches
n_years <- nrow(inf_mat)
mean_trajectory <- rowMeans(inf_mat, na.rm = TRUE)

# Plot individual stochastic trajectories
matplot(1:n_years, inf_mat, type = "l", lty = 1, 
        col = adjustcolor("grey", alpha.f = 0.5),
        xlab = "Year", ylab = "Number of Infected Patches",
        main = sprintf("Metapopulation Dynamics (R0 = %.1f, rho = %.1f)", 
                       test_params$R0, test_params$rho),
        ylim = c(0, test_params$n_patches))

# Overlay mean trajectory of surviving sims
lines(1:n_years, mean_trajectory, col = "royalblue", lwd = 3)

# Overlay quasi-equilibrium reference line (last 100 years)
start_qe <- n_years - 100 + 1
segments(x0 = start_qe, x1 = n_years, y0 = qe_mean, y1 = qe_mean,
         col = "firebrick", lwd = 3, lty = 2)

# Add legend for trajectories
legend("topleft",
       legend = c("Individual Sims", 
                  "Mean Trajectory (including extincted)", 
                  sprintf("Quasi-Eq Mean: %.1f", qe_mean)),
       col = c("grey", "royalblue", "firebrick"),
       lty = c(1, 1, 2), lwd = c(1, 3, 3), bty = "n")

# Add legend for extinction metrics
legend("topright",
       legend = sprintf("Extinction Rate: %.0f%% (%d/%d)", 
                        ext_rate * 100, n_ext, nsim_test),
       text.col = ifelse(ext_rate > 0.5, "firebrick", "black"), 
       text.font = 2, bty = "n", cex = 1.1)
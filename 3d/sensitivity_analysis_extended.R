# Source the original script to get getDFE and getEE functions
source("sensitivity_analysis.R")

# Define parameter ranges: c0 and r
c0_values <- c(0.01, 0.25, 0.5, 0.75, 1)
r_values <- c(0.01, 0.25, 0.5, 0.75, 1)

# R0 and alpha ranges 
R0_values <- seq(1.0, 3.0, by = 0.1)
alpha_values <- 5 * 10^seq(-6, -2, by = 0.5)

# Color palette
colors <- rainbow(length(alpha_values))

# Create output directory for plots
if (! dir.exists("sensitivity_plots")) {
  dir.create("sensitivity_plots")
}

# Open PDF device
pdf("all_sensitivity_plots.pdf", width = 10, height = 7)

# Store all results
all_extended_results <- list()

# Loop through c0 values
for (c0_val in c0_values) {
  
  # Loop through r values
  for (r_val in r_values) {
    
    cat("Processing c0 =", c0_val, ", r =", r_val, "\n")
    
    # Initialize plot
    plot(NULL, xlim = range(R0_values), ylim = c(0, 1),
         xlab = "R0", ylab = "p",
         main = sprintf("Equilibrium p vs R0 (c0 = %.2f, r = %.2f)", c0_val, r_val),
         las = 1)
    grid()
    
    # Store results for this combination
    combo_results <- list()
    legend_labels <- list()
    legend_colors <- list()
    
    # Loop through each alpha value
    for (j in seq_along(alpha_values)) {
      alpha_val <- alpha_values[j]
      
      # Initialize results dataframe for this alpha
      results_R0 <- data.frame(R0 = R0_values, p_EE = NA, NI_EE = NA, NS_EE = NA)
      
      cat("  Processing alpha =", alpha_val, "\n")
      
      # Loop through R0 values
      for (i in seq_along(R0_values)) {
        tryCatch({
          ee_temp <- getEE(R0 = R0_values[i], alpha = alpha_val, 
                           c0 = c0_val, r = r_val)
          if (ee_temp$convergence < 1e-2) { 
            results_R0[i, c("p_EE", "NI_EE", "NS_EE")] <- c(ee_temp$p, ee_temp$ni, ee_temp$ns)
          }
        }, error = function(e) cat("    EE failed for R0 =", R0_values[i], "\n"))
      }
      
      combo_results[[j]] <- results_R0
      
      # Plot valid data points
      valid_data <- ! is.na(results_R0$p_EE)
      if (sum(valid_data) > 0) {
        R0_valid <- results_R0$R0[valid_data]
        p_EE_valid <- results_R0$p_EE[valid_data]
        
        # Add points
        points(R0_valid, p_EE_valid, pch = 19, col = colors[j], cex = 0.8)
        
        # Add smooth line
        if (length(R0_valid) > 3) {
          lines(smooth.spline(R0_valid, p_EE_valid), col = colors[j], lwd = 2)
        } else {
          lines(R0_valid, p_EE_valid, col = colors[j], lwd = 2)
        }
        
        # Add to legend
        legend_labels[[length(legend_labels) + 1]] <- sprintf("alpha = %.0e", alpha_val)
        legend_colors[[length(legend_colors) + 1]] <- colors[j]
      }
    }
    
    # Add legend
    if (length(legend_labels) > 0) {
      legend("topleft", 
             legend = unlist(legend_labels), 
             col = unlist(legend_colors), 
             lwd = 2, pch = 19, cex = 0.8, bg = "white")
    }
    
    # Store results
    param_key <- sprintf("c0_%.2f_r_%.2f", c0_val, r_val)
    names(combo_results) <- paste0("alpha_", formatC(alpha_values, format = "e", digits = 1))
    all_extended_results[[param_key]] <- combo_results
  }
}

# Close PDF device
dev.off()
cat("All plots saved to: all_sensitivity_plots.pdf\n")

# Save results
save(all_extended_results, R0_values, alpha_values, c0_values, r_values,
     file = "sensitivity_analysis_extended_results.RData")
cat("Data saved to: sensitivity_analysis_extended_results.RData\n")
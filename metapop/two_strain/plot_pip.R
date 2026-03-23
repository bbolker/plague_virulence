library(ggplot2)

fn_sum <- "outputs_pip/sim_pip_alpha_rho.rds"
dd_out <- readRDS(fn_sum)

# Calculate Invasion Probability (1 - Extinction Rate of Strain 2)
dd_out$invasion_prob <- 1 - dd_out$extinction_rate_2

# Extract all unique combinations of alpha and rho
params_grid <- unique(dd_out[, c("alphavec", "rhovec")])

output_pdf <- "outputs_pip/pip_alpha_rho_multipage.pdf"

# Open PDF device
pdf(output_pdf, width = 8, height = 7)

for (i in 1:nrow(params_grid)) {
  cur_alpha <- params_grid$alphavec[i]
  cur_rho <- params_grid$rhovec[i]
  
  # Subset data for the current parameter combination
  sub_dat <- subset(dd_out, alphavec == cur_alpha & rhovec == cur_rho)
  
  # Generate PIP for the subset
  p <- ggplot(sub_dat, aes(x = R01vec, y = R02vec, fill = invasion_prob)) +
    geom_tile() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) + 
    theme_bw() +
    labs(
      title = sprintf("PIP: alpha = %e, rho = %.1f", cur_alpha, cur_rho),
      x = "Resident R0",
      y = "Invader R0",
      fill = "Invasion\nProb"
    ) +
    coord_fixed() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Print plot to create a new page in the PDF
  print(p)
}

# Close PDF device
invisible(dev.off())

cat("Multi-page PIP PDF saved to", output_pdf, "\n")
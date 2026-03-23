library(ggplot2)

# Load simulation results
fn_sum <- "outputs_pip/sim_pip_alpha_rho.rds"
dd_out <- readRDS(fn_sum)

# Invasion starts at year 100, simulation ends at year 500
# Total possible duration after invasion is 400 years
invade_start <- 100
max_duration <- 400

# Compute Hybrid Metric with Log-transformed Time for the extinction region:
# 1. Extinction Rate >= 90%: Metric 0-1 (Log-scale of years survived after invasion)
# 2. Extinction Rate < 90%: Metric 1.1-2 (Linear Invasion Probability 10%-100%)
dd_out$time_after_invasion <- dd_out$mean_extinct_time_2 - invade_start

# Apply log10 to spread out early extinctions (1yr, 10yrs, etc.)
# We use log10(time + 1) to handle the 0-year case
dd_out$hybrid_metric <- ifelse(
  dd_out$extinction_rate_2 >= 0.9,
  log10(pmax(1, dd_out$time_after_invasion)) / log10(max_duration), # Maps [1, 400] to [0, 1]
  1 + (1 - dd_out$extinction_rate_2)                               # Maps [0.1, 1.0] prob to (1.1, 2]
)

params_grid <- unique(dd_out[, c("alphavec", "rhovec")])
output_pdf <- "outputs_pip/pip_hybrid.pdf"

pdf(output_pdf, width = 8.5, height = 10)

for (i in 1:nrow(params_grid)) {
  cur_alpha <- params_grid$alphavec[i]
  cur_rho <- params_grid$rhovec[i]
  
  sub_dat <- subset(dd_out, alphavec == cur_alpha & rhovec == cur_rho)
  
  p <- ggplot(sub_dat, aes(x = R01vec, y = R02vec, fill = hybrid_metric)) +
    geom_tile() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "white", alpha = 0.7) +
    # Custom scale with log-spaced time breaks and linear prob breaks
    scale_fill_viridis_c(
      option = "magma",
      limits = c(0, 2),
      breaks = c(
        0,                               # 1 year survived
        log10(11)/log10(400),           # 10 years survived
        log10(101)/log10(400),          # 100 years survived
        1.0,                             # 400 years / 10% Probability
        1.5,                             # 50% Probability
        2.0                              # 100% Probability
      ),
      labels = c(
        "Immediate\n(<1 yr)", 
        "10 yrs", 
        "100 yrs", 
        "400 yrs / 10% Prob\n(Threshold)", 
        "50% Prob", 
        "100% Prob\n(Invasion Success)"
      )
    ) + 
    theme_bw() +
    labs(
      title = sprintf("Hybrid PIP: alpha = %e, rho = %.1f", cur_alpha, cur_rho),
      x = "Resident R0",
      y = "Invader R0",
      fill = "Invasion Status",
      caption = "Extinction region (0-1) is log-scaled to emphasize early extinction dynamics."
    ) +
    coord_fixed() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.key.height = unit(2.5, "cm"),
      legend.text = element_text(size = 7),
      panel.grid = element_blank()
    )
  
  print(p)
}

invisible(dev.off())
cat("PIP PDF with log-time scale saved to:", output_pdf, "\n")
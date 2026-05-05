#!/usr/bin/env Rscript
library(optparse)
op.parser <- OptionParser(prog="plot_pip",
                          option_list = list(
                            make_option(c("-i", "--input"), "store",
                                        help = "input data file",
                                        default = "outputs_pip/sim_pip_alpha_rho.rds"),
                            make_option(c("-o", "--output"), "store",
                                        help = "output plot file",
                                        default = "outputs_pip/pip_hybrid.pdf"))
                          )

opt <- parse_args(op.parser)
library(ggplot2)
dd_out <- readRDS(opt$input)

invade_start <- 100
max_duration <- 400

dd_out$time_after_invasion <- dd_out$mean_extinct_time_2 - invade_start

dd_out$hybrid_metric <- ifelse(
  dd_out$extinction_rate_2 >= 0.9,
  log10(pmax(1, dd_out$time_after_invasion)) / log10(max_duration),
  1 + (1 - dd_out$extinction_rate_2)
)

params_grid <- unique(dd_out[, c("alphavec", "rhovec")])
output_pdf <- opt$output

## FIXME: optionally facet instead?
plot_fun <- function(sub_dat, facet = FALSE, title = NULL) {
  p <- ggplot(sub_dat, aes(x = R01vec, y = R02vec, fill = hybrid_metric)) +
    geom_raster() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "white", alpha = 0.7) +
    scale_fill_viridis_c(
      option = "magma",
      limits = c(0, 2),
      breaks = c(
        0,
        log10(11)/log10(400),
        log10(101)/log10(400),
        1.0,
        1.5,
        2.0
      ),
      labels = c(
        "Immediate\n(<1 yr)(Mean extinction time)", 
        "10 yrs", 
        "100 yrs", 
        "400 yrs / 10% Prob\n", 
        "50% ", 
        "\n100% (persistence probability)\n"
      )
    ) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(
      title = title,
      x = "Resident R0",
      y = "Invader R0",
      fill = "Invasion capability"
    ) +
    coord_fixed() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.key.height = unit(2.5, "cm"),
      legend.text = element_text(size = 7),
      panel.grid = element_blank()
    )

  if (facet) {
    p <- p + facet_grid(alphavec ~ rhovec, labeller = label_both)
  }
  return(p)
}

plot_fun(dd_out, facet = TRUE)
ggsave(output_pdf, width = 8.5, height = 10)

## alternative approach
if (FALSE) {
  pdf(output_pdf, width = 8.5, height = 10)
  for (i in 1:nrow(params_grid)) {
    cur_alpha <- params_grid$alphavec[i]
    cur_rho <- params_grid$rhovec[i]
    title <- sprintf("PIP: alpha = %e, rho = %.1f", cur_alpha, cur_rho)
    sub_dat <- subset(dd_out, alphavec == cur_alpha & rhovec == cur_rho)
    print(plot_fun(sub_dat, title = title))
  }
  invisible(dev.off())
  cat("PIP PDF with log-time scale saved to:", output_pdf, "\n")
}

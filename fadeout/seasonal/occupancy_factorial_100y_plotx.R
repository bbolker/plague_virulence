library(ggplot2)
library(here)
library(dplyr)
library(readr)

amps <- c(seasonal = 0.40, constant = 0.00)
summary_years <- 50
facet_labels <- labeller(
  R0 = function(x) paste("R0 =", x),
  K = function(x) paste("K =", format(as.numeric(x), big.mark = ",")),
  alpha = function(x) paste0("alpha = ", x, "/day")
)

scenario_labels <- c(
  seasonal = sprintf("seasonal\n(amplitude = %g)", amps[["seasonal"]]),
  constant = "constant"
)

theme_set(theme_bw(base_size=14))
zmargin <- theme(panel.spacing = grid::unit(0, "pt"))
## Okabe-Ito minus black and yellow
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

outdir <- here::here("fadeout", "output", "occupancy_factorial_100y")

save_figure <- function(plot, name, width, height) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path(outdir, "figures", paste0(name, ".", ext)),
           plot, width = width, height = height, dpi = 300)
  }
}

late_summary <- read_csv(file.path(outdir, "data", "late_summary.csv")) |>
    filter(K < 3e4) |>
    mutate(across(c(lwr, median, upr), ~ . /200)) |>
    mutate(across(alpha,
                  ~ factor(., labels = c("slow", "medium", "fast"))))

p_late <- ggplot(
  late_summary,
  aes(
    x = factor(K),
    y = median,
    ymin = lwr,
    ymax = upr,
    colour = factor(R0),
    shape = factor(alpha)
  )
) +
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.3) +
  facet_grid(. ~ scenario, labeller = facet_labels) +
  scale_colour_manual(values = okabe_ito, labels = scenario_labels, name = expression("transmission "*(R[0]))) +
    scale_shape_manual(values = c(16, 17, 15),
                       name = "colonization rate") +
  labs(
    x = "Carrying capacity per patch",
    y = "Proportion of patches occupied",
    title = sprintf(
      "Occupancy over the last %d years: median and 95%% range",
      summary_years
    )
  ) +
  zmargin +
  theme(legend.position = "top")

p_late
save_figure(p_late, "occupancy_factorial_late_summary", width = 8, height = 4)

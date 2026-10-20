## Figure for s0_persistence_check.R: time to global extinction and minimum
## patch occupancy against the initial susceptible density, one column per R0.
##
## Reads data/s0_persistence.rds. Runs that persist the full 100 years are
## drawn at 100 years as open points.
##
## Usage:
##   Rscript fadeout/seasonal/s0_persistence_check_plots.R

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

theme_set(theme_bw())

## Facet columns abut, facet rows are separated
zmargin <- theme(
  panel.spacing.x = grid::unit(0, "pt"),
  panel.spacing.y = grid::unit(6, "pt")
)

outdir <- here::here("fadeout", "output", "s0_persistence_check")
dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)

years <- 100

dat <- readRDS(file.path(outdir, "data", "s0_persistence.rds")) |>
  mutate(
    persisted = is.na(ext_day),
    ext_years = ifelse(persisted, years, ext_day / 365)
  )

long <- dat |>
  select(R0, f, persisted, ext_years, min_occ) |>
  pivot_longer(c(ext_years, min_occ)) |>
  mutate(name = factor(name, levels = c("ext_years", "min_occ"),
                       labels = c("years to global\nextinction (max 100)",
                                  "minimum fraction of\npatches occupied")))

gg <- ggplot(long, aes(f, value)) +
  geom_point(aes(shape = persisted), size = 1, alpha = 0.5,
             colour = "#2a78d6",
             position = position_jitter(width = 0.02, height = 0, seed = 1)) +
  stat_summary(fun = median, geom = "line", linewidth = 0.7,
               colour = "grey15") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                     labels = c("went extinct", "persisted 100 y"),
                     name = NULL) +
  ## Breaks kept off the panel edges so labels of abutting columns don't collide
  scale_x_continuous(breaks = c(0.6, 1, 1.4)) +
  facet_grid(name ~ R0, scales = "free_y", switch = "y",
             labeller = labeller(R0 = \(x) paste0("R0 = ", x))) +
  labs(x = expression("Initial S as a multiple of " * S^"*" + I^"*"),
       y = NULL,
       subtitle = paste("K = 5000, alpha = 5e-05/day, seasonal amplitude 0.2;",
                        "homogeneous start, 20 runs per value; line = median")) +
  theme(strip.placement = "outside", legend.position = "top") +
  zmargin

for (ext in c("pdf", "png")) {
  ggsave(file.path(outdir, "figures", paste0("s0_persistence.", ext)),
         gg, width = 12, height = 5, dpi = 300)
}

cat("Figures: ", file.path(outdir, "figures"), "\n", sep = "")

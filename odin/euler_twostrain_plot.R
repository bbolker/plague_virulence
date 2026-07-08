library(dplyr)
library(ggplot2); theme_set(theme_bw())
zmargin <- theme(panel.spacing=grid::unit(0,"lines"))

input_dir <- here::here("odin/sharcnet/outputs")
##input_fn <- "euler_twostrain_mini2_singlepatchintro.rds"
input_fn <- "euler_twostrain_singlepatchintro.rds"

## FIXME: should retrieve from run attributes
dt <- 0.1
invade_start <- 100
max_duration <- 500

log10alpha_labeller <- function(labels) {
  list(lapply(sprintf("%s == 10^{%1.1f}", names(labels)[1],
                      log10(labels[[1]])), function(x) parse(text = x)[[1]]))
}
class(log10alpha_labeller) <- c("function", "labeller")

hybrid_metric_fun <- function(extinction_rate,
                          time_after_invasion,
                          max_duration = 200,
                          brk = 0.9) {
  ifelse(
    extinction_rate >= brk,
    log10(pmax(1, time_after_invasion)) / log10(max_duration),
    1 + (1 - extinction_rate) ## ??
  )
}

resident_ok_cutoff <- 0.9

dd <- readRDS(file.path(input_dir, input_fn)) |>
  mutate(
    time_after_invasion = mean_ext_time.I2 * dt - invade_start,
    resident_valid = resident_alive_at_intro_prob >= resident_ok_cutoff,
    hybrid_metric_raw = hybrid_metric_fun(ext_prob.I2, time_after_invasion),
    hybrid_metric = ifelse(resident_valid, hybrid_metric_raw, NA_real_)
  )

params_grid <- unique(dd[, c("alpha", "K")])

## FIXME: optionally facet instead?
plot_fun <- function(sub_dat, facet = FALSE, title = NULL) {
  p <- ggplot(sub_dat, aes(x = R01, y = R02, fill = hybrid_metric)) +
    geom_raster() +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "white", alpha = 0.7) +
    scale_fill_viridis_c(
      option = "magma",
      na.value = "grey80",
      limits = c(0, 2),
      breaks = c(
        0,
        log10(11)/log10(max_duration),
        log10(101)/log10(max_duration),
        1.0,
        1.5,
        2.0
      ),
      labels = c(
        "Immediate\n(<1 yr)(Mean extinction time)", 
        "10 yrs", 
        "100 yrs", 
        sprintf("%d yrs / 10%% Prob\n", max_duration),
        "50% ", 
        "\n100% (persistence probability)\n"
      )
    ) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(
      title = title,
      x = expression("Resident "*R[0]),
      y = expression("Invader "*R[0]),
      fill = "Invasion capability"
    ) +
    coord_fixed() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.key.height = unit(2.5, "cm"),
      legend.text = element_text(size = 7),
      panel.grid = element_blank()
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, alpha = 0.8)

  if (facet) {
    p <- p + facet_grid(alpha ~ K, labeller = log10alpha_labeller)
  }
  return(p)
}

summary(dd$mean_ext_time.I2) ## uh-oh
summary(dd$ext_prob.I2) ## ALL extinct ... ??? is this expected?
summary(dd$ext_prob.I1) ## sometimes extinct

ggplot(dd, aes(R01, R02, fill = mean_ext_time.I2)) +
  geom_raster() +
  facet_grid(alpha ~ K) +
  scale_fill_viridis_c(trans = "log10") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))
  
gg0 <- (plot_fun(dd, facet = TRUE)
  + zmargin
  + theme(strip.text.y = element_text(angle = 0))
)

print(gg0)

ggsave("euler_twostrain_pip.png", width = 10, height = 8)
ggsave("euler_twostrain_pip.pdf", width = 10, height = 8)

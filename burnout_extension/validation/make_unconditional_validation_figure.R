source('validation/R/theory.R')

stopifnot(
  requireNamespace('data.table', quietly = TRUE),
  requireNamespace('ggplot2', quietly = TRUE),
  requireNamespace('scales', quietly = TRUE)
)

library(data.table)
library(ggplot2)

out_dir <- 'validation/figures/paper'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theta_values <- seq(0, 1, by = .25)
pal_theta <- c(`0` = '#0072B2', `0.25` = '#56B4E9', `0.5` = '#009E73',
               `0.75` = '#E69F00', `1` = '#D55E00')
clip_prob <- function(x, eps = 1e-3) pmin(1 - eps, pmax(eps, x))
prob_breaks <- c(.001, .01, .05, .25, .5, .75, .95, .99, .999)
theme_paper <- theme_bw(base_size = 10.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = .25, colour = 'grey88'),
        strip.background = element_rect(fill = 'grey93', colour = 'grey55'),
        strip.text = element_text(face = 'bold'),
        legend.position = 'bottom',
        legend.box = 'horizontal',
        plot.margin = margin(5.5, 7, 5.5, 5.5))
x_label_R0 <- bquote(.(intToUtf8(0x211B))[0] - 1 ~ '(log scale)')
fig1_math_family <- 'serif'
if (capabilities('aqua')) {
  grDevices::quartzFonts(STIX = grDevices::quartzFont(
    c('STIXGeneral-Regular', 'STIXGeneral-Bold',
      'STIXGeneral-Italic', 'STIXGeneral-BoldItalic')))
  fig1_math_family <- 'STIX'
}

# Figure 1: direct validation in the largest, cleanest common stochastic scan.
source('validation/R/validation_plot_data.R')
validation_input <- load_validation_plot_data()
sim <- validation_input$sim
curve_grid <- CJ(rho = .01, theta = theta_values,
                 K = c(1e6, 1e7, 1e8, 1e9),
                 R0 = 1 + exp(seq(log(.03), log(validation_input$Rmax-1),
                                  length.out = if(validation_input$extended)800 else 400)))
curve_grid[, prediction := vapply(seq_len(.N), function(i)
  stable_validation_probability(R0[i],rho[i],theta[i],K[i],conditional=FALSE),numeric(1))]
validation_x_breaks <- if(validation_input$extended) c(.03,.1,.3,1,3,10,30,100) else
  c(.03,.05,.1,.2,.5,1,2,5)
validation_x_labels <- if(validation_input$extended) c('0.03','0.1','0.3','1','3','10','30','100') else
  c('0.03','0.05','0.10','0.20','0.50','1','2','5')

sim_plot <- copy(sim)
sim_plot[, `:=`(estimate_display = clip_prob(P_unconditional),
                low_display = clip_prob(uncond_low),
                high_display = clip_prob(uncond_high),
                theta_label = factor(theta, levels = theta_values))]
curve_plot <- copy(curve_grid)
curve_plot[, `:=`(prediction_display = clip_prob(prediction),
                  theta_label = factor(theta, levels = theta_values))]

nonfizzle <- unique(curve_grid[, .(R0, K)])
nonfizzle[, probability := 1 - 1/R0]
fig1 <- ggplot() +
  geom_line(data = nonfizzle,
            aes(R0 - 1, qlogis(clip_prob(probability)), linetype = 'Non-fizzle'),
            colour = 'grey35', linewidth = .65) +
  scale_linetype_manual(values = c('Non-fizzle' = 'dashed'), name = NULL) +
  geom_line(data = curve_plot,
            aes(R0 - 1, qlogis(prediction_display), colour = theta_label),
            linewidth = .9) +
  geom_errorbar(data = sim_plot,
                aes(R0 - 1, ymin = qlogis(low_display), ymax = qlogis(high_display),
                    colour = theta_label),
                width = 0, linewidth = .35, alpha = .75) +
  geom_point(data = sim_plot,
             aes(R0 - 1, qlogis(estimate_display), colour = theta_label),
             shape = 16, size = 1.8, stroke = .45) +
  facet_wrap(~K, ncol = 2,
             labeller = as_labeller(function(x)
               paste0('italic(K) == 10^', round(log10(as.numeric(x)))),
               default = label_parsed)) +
  scale_x_log10(breaks = validation_x_breaks, labels = validation_x_labels) +
  scale_y_continuous(breaks = qlogis(prob_breaks),
                     labels = format(prob_breaks, trim = TRUE),
                     limits = qlogis(c(.001, .999)),
                     expand = expansion(mult = c(.015, .015))) +
  scale_colour_manual(values = pal_theta, name = expression(theta)) +
  labs(x = x_label_R0,
       y = 'Unconditional persistence probability') +
  guides(colour = guide_legend(override.aes = list(linewidth = .9, size = 2.2,
                                                 shape = 16))) + theme_paper +
  theme(axis.title.x = element_text(family = fig1_math_family))

# Quartz embeds native macOS fonts without requiring XQuartz/Cairo libraries.
fig1_device <- if (capabilities('aqua')) {
  function(filename, width, height, ...) {
    grDevices::quartz(file = filename, type = 'pdf',
                     width = width, height = height, ...)
  }
} else grDevices::cairo_pdf
fig1_file <- file.path(out_dir, 'unconditional_persistence_validation.pdf')
fig1_tmp <- tempfile(pattern = '.fig1-', tmpdir = out_dir, fileext = '.pdf')
ggsave(fig1_tmp, fig1, width = 7.15, height = 5.45, device = fig1_device)
stopifnot(file.rename(fig1_tmp, fig1_file))

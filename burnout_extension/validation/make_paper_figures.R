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

pal_theta <- c(`0` = '#0072B2', `0.5` = '#E69F00', `1` = '#009E73')
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
x_label_R0 <- paste0(intToUtf8(0x211B), intToUtf8(0x2080),
                     ' - 1 (log scale)')

# Figure 1: direct validation in the largest, cleanest common stochastic scan.
sim <- fread('validation/data/fig11_scan_results.csv')
curve_grid <- unique(fread('validation/data/fig11_scan_analytic.csv')[,
  .(rho, theta, K, R0)])
curve_grid[, prediction := vapply(seq_len(.N), function(i)
  bi_quantities(R0[i], rho[i], theta[i], K[i])$P_conditional, numeric(1))]

sim_plot <- copy(sim)
sim_plot[, `:=`(estimate_display = clip_prob(P_conditional),
                low_display = clip_prob(cond_low),
                high_display = clip_prob(cond_high),
                theta_label = factor(theta, levels = c(0, .5, 1)))]
curve_plot <- copy(curve_grid)
curve_plot[, `:=`(prediction_display = clip_prob(prediction),
                  theta_label = factor(theta, levels = c(0, .5, 1)))]

fig1 <- ggplot() +
  geom_line(data = curve_plot,
            aes(R0 - 1, qlogis(prediction_display), colour = theta_label),
            linewidth = .9) +
  geom_errorbar(data = sim_plot,
                aes(R0 - 1, ymin = qlogis(low_display), ymax = qlogis(high_display),
                    colour = theta_label),
                width = 0, linewidth = .35, alpha = .75) +
  geom_point(data = sim_plot,
             aes(R0 - 1, qlogis(estimate_display), colour = theta_label,
                 shape = theta_label),
             size = 1.8, stroke = .45) +
  facet_wrap(~K, ncol = 2,
             labeller = labeller(K = function(x)
               paste0('K = ', format(as.numeric(x), scientific = TRUE)))) +
  scale_x_log10(breaks = c(.03, .05, .1, .2, .5, 1, 2, 5),
                labels = c('.03', '.05', '.10', '.20', '.50', '1', '2', '5')) +
  scale_y_continuous(breaks = qlogis(prob_breaks),
                     labels = format(prob_breaks, trim = TRUE),
                     limits = qlogis(c(.001, .999)),
                     expand = expansion(mult = c(.015, .015))) +
  scale_colour_manual(values = pal_theta, name = expression(theta)) +
  scale_shape_manual(values = c(16, 17, 15), name = expression(theta)) +
  labs(x = x_label_R0,
       y = 'Conditional persistence probability') +
  guides(colour = guide_legend(override.aes = list(linewidth = .9, size = 2.2)),
         shape = 'none') + theme_paper +
  theme(axis.title.x = element_text(family = 'Cambria Math'))

ggsave(file.path(out_dir, 'main_persistence_validation.pdf'), fig1,
       width = 7.15, height = 5.45, device = cairo_pdf)

# Figure 2: a directly interpretable population-size threshold.  Because
# B_BI is linear in K, K_50 = log(2)/(B_BI/K).
theta_values <- seq(0, 1, by = .25)
threshold <- CJ(rho = c(.005, .01, .02), theta = theta_values,
                R0 = exp(seq(log(1.05), log(6), length.out = 260)))
threshold[, K50 := vapply(seq_len(.N), function(i) {
  unit <- bi_quantities(R0[i], rho[i], theta[i], 1)$B
  log(2) / unit
}, numeric(1))]
threshold[, K50_display := fifelse(K50 >= 1e2 & K50 <= 1e14, K50, NA_real_)]
threshold[, theta_label := factor(theta, levels = theta_values)]
threshold[, rho_label := paste0('rho == ', rho)]

pal5 <- c('#0072B2', '#56B4E9', '#009E73', '#E69F00', '#D55E00')
fig2 <- ggplot(threshold, aes(R0 - 1, K50_display, colour = theta_label)) +
  annotate('rect', xmin = .05, xmax = 5, ymin = 1e6, ymax = 1e9,
           fill = 'grey75', alpha = .22) +
  geom_line(linewidth = .9, na.rm = TRUE) +
  facet_wrap(~rho_label, nrow = 1, labeller = label_parsed) +
  scale_x_log10(breaks = c(.05, .1, .2, .5, 1, 2, 5),
                labels = c('.05', '.10', '.20', '.50', '1', '2', '5')) +
  scale_y_log10(breaks = 10^(2:14), labels = scales::label_math(10^.x),
                limits = c(1e2, 1e14)) +
  scale_colour_manual(values = pal5, name = expression(theta)) +
  labs(x = x_label_R0,
       y = expression('50% persistence threshold, '~K[50])) +
  theme_paper + theme(legend.key.width = unit(1.1, 'cm'),
                      axis.title.x = element_text(family = 'Cambria Math'))

ggsave(file.path(out_dir, 'theta_persistence_threshold.pdf'), fig2,
       width = 7.15, height = 3.35, device = cairo_pdf)

# Supplementary diagnostic: isolate deterministic entry and Laplace errors on
# the cached ODE grid, and contrast them with the stochastic error by rho.
ode <- fread('validation/data/K_R0_scan.csv')[boundary_choice == 'sqrt']
ode[, BI := vapply(seq_len(.N), function(i)
  bi_quantities(R0[i], rho[i], theta[i], K[i])$P_conditional, numeric(1))]
ode_ok <- ode[status == 'OK' & is.finite(P1_ref)]
diag_ode <- melt(ode_ok[, .(
  `BI closed form` = mean(abs(BI - P1_ref), na.rm = TRUE),
  `Laplace step` = mean(abs(P1_L_ODE - P1_ref), na.rm = TRUE),
  `First-order entry` = mean(abs(P1_K_first - P1_ref), na.rm = TRUE),
  `Second-order entry` = mean(abs(P1_K_second - P1_ref), na.rm = TRUE)),
  by = rho], id.vars = 'rho', variable.name = 'component', value.name = 'MAE')

sim_all <- rbindlist(list(
  fread('validation/data/stochastic_scan_results.csv'),
  fread('validation/data/fig11_scan_results.csv')), fill = TRUE)
sim_all[, prediction := vapply(seq_len(.N), function(i)
  bi_quantities(R0[i], rho[i], theta[i], K[i])$P_conditional, numeric(1))]
diag_sim <- sim_all[, .(MAE = mean(abs(prediction - P_conditional))), by = rho]

figS1a <- ggplot(diag_sim, aes(rho, MAE)) +
  geom_line(linewidth = .8, colour = '#7A0177') +
  geom_point(size = 2, colour = '#7A0177') +
  scale_x_continuous(breaks = sort(unique(diag_sim$rho))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .08))) +
  labs(x = expression(rho), y = 'Mean absolute probability error',
       title = '(a) Complete stochastic comparison') + theme_paper +
  theme(legend.position = 'none', plot.title = element_text(size = 10.5))

figS1b <- ggplot(diag_ode, aes(rho, MAE, colour = component, shape = component)) +
  geom_line(linewidth = .75) + geom_point(size = 1.9) +
  scale_colour_manual(values = c('#CC79A7', '#D55E00', '#56B4E9', '#009E73')) +
  scale_x_continuous(breaks = sort(unique(diag_ode$rho))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .08))) +
  labs(x = expression(rho), y = 'Mean absolute probability difference',
       colour = NULL, shape = NULL,
       title = '(b) ODE/Kendall error isolation') + theme_paper +
  theme(plot.title = element_text(size = 10.5),
        legend.position = 'inside', legend.position.inside = c(.98, .98),
        legend.justification = c(1, 1), legend.direction = 'vertical',
        legend.background = element_rect(fill = scales::alpha('white', .88),
                                         colour = 'grey70'),
        legend.key.width = unit(.55, 'cm'), legend.text = element_text(size = 7.5))

cairo_pdf(file.path(out_dir, 'accuracy_diagnostics_supplement.pdf'),
          width = 7.15, height = 3.5)
grid::grid.newpage()
layout <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = layout))
print(figS1a, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(figS1b, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

# Reproducible numerical summary used in the manuscript text.
summary_rows <- rbindlist(list(
  sim_all[, .(quantity = 'simulation_MAE', value = mean(abs(prediction - P_conditional))),
          by = .(rho)],
  sim_all[rho <= .02 & R0 >= 1.15,
          .(quantity = 'simulation_MAE_separated',
            value = mean(abs(prediction - P_conditional))), by = .(rho, theta)],
  ode_ok[, .(quantity = 'second_order_entry_MAE',
             value = mean(abs(P1_K_second - P1_ref), na.rm = TRUE)), by = rho]
), fill = TRUE)
fwrite(summary_rows, 'validation/data/paper_numerical_summary.csv')

message('Paper figures and numerical summary written to ', out_dir)

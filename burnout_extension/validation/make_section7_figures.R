source('validation/R/theory.R')
source('validation/R/kendall.R')
source('validation/R/extrema.R')

stopifnot(requireNamespace('data.table', quietly = TRUE),
          requireNamespace('ggplot2', quietly = TRUE),
          requireNamespace('scales', quietly = TRUE))
library(data.table)
library(ggplot2)

out_dir <- 'validation/figures/paper'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
scan_version <- 'section7_analytic_extrema_v2'
theta5 <- seq(0, 1, by = .25)
theta3 <- c(0, .5, 1)
pal5 <- c('#0072B2', '#56B4E9', '#009E73', '#E69F00', '#D55E00')
x_label_R0 <- paste0(intToUtf8(0x211B), intToUtf8(0x2080))

theme_paper <- theme_bw(base_size = 9.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = .22, colour = 'grey88'),
        strip.background = element_rect(fill = 'grey93', colour = 'grey55'),
        strip.text = element_text(face = 'bold'), legend.position = 'bottom',
        plot.margin = margin(4, 5, 4, 4))

log_conditional <- function(log_B) {
  ans <- numeric(length(log_B)); small <- log_B < -20; large <- log_B > 20
  ans[small] <- log_B[small]; ans[large] <- 0; mid <- !(small | large)
  ans[mid] <- log(-expm1(-exp(log_B[mid]))); ans
}
prob_conditional <- function(log_B) {
  ans <- numeric(length(log_B)); large <- log_B > 20; ans[large] <- 1
  ans[!large] <- -expm1(-exp(log_B[!large])); ans
}
core_one <- function(R0, theta) {
  xf <- x_final(R0); xs <- 1 / R0
  C <- C_explicit(R0, theta); delta <- 1 - R0 * xf
  data.table(R0 = R0, theta = theta, xf = xf, delta = delta,
             hf = h_theta(xf, theta), C = C, logC = log(C),
             hstar = h_theta(xs, theta),
             DeltaA = action_diff(xf, xs, R0, theta))
}
core_grid <- function(R0_values, theta_values) {
  g <- CJ(R0 = R0_values, theta = theta_values)
  rbindlist(lapply(seq_len(nrow(g)), function(i) core_one(g$R0[i], g$theta[i])))
}
add_probabilities <- function(d, K) {
  z <- copy(d)
  z[, log_B := log(K) + logC + .5 * log(R0 * rho * hstar / (2 * pi)) -
      DeltaA / rho]
  z[, P_conditional := prob_conditional(log_B)]
  z[, log_P_unconditional := log1p(-1 / R0) + log_conditional(log_B)]
  z[, P_unconditional := exp(log_P_unconditional)]
  z
}
add_entry_diagnostic <- function(d, K) {
  z <- copy(d)
  z[, ystar := rho * hstar]
  z[, ybl := sqrt(ystar / K)]
  z[, `:=`(a = hf / delta, b = R0 * xf / delta)]
  z[, M1 := rho * log(C / ybl) + (b / a) * ybl]
  z[, finite_entry := is.finite(M1) & M1 > 0 & M1 < DeltaA]
  z
}
read_cache <- function(path) {
  if (!file.exists(path)) return(NULL)
  z <- fread(path)
  if (!'scan_version' %in% names(z) ||
      !identical(unique(z$scan_version), scan_version)) return(NULL)
  z
}
cross_join_rho <- function(core, rho_values) {
  a <- copy(core); b <- data.table(rho = rho_values)
  a[, join_key__ := 1L]; b[, join_key__ := 1L]
  z <- merge(a, b, by = 'join_key__', allow.cartesian = TRUE)
  z[, join_key__ := NULL]
  z
}

# Figure 1: stochastic validation with five analytical recruitment shapes.
sim_large <- fread('validation/data/fig11_scan_results.csv')
R_validation <- sort(unique(c(seq(1.03, 1.3, length.out = 90),
                              seq(1.3, 6, length.out = 170))))
validation_core <- core_grid(R_validation, theta5)
validation_curves <- rbindlist(lapply(sort(unique(sim_large$K)), function(Kv) {
  z <- copy(validation_core); z[, rho := .01]; z <- add_probabilities(z, Kv)
  z[, K := Kv]; z
}))
clip_prob <- function(x, eps = 1e-3) pmin(1 - eps, pmax(eps, x))
prob_breaks <- c(.001, .01, .05, .25, .5, .75, .95, .99, .999)
validation_curves[, `:=`(probability_display = clip_prob(P_conditional),
                         theta_label = factor(theta, levels = theta5))]
sim_large[, `:=`(estimate_display = clip_prob(P_conditional),
                 low_display = clip_prob(cond_low), high_display = clip_prob(cond_high),
                 theta_label = factor(theta, levels = theta5))]
fig1 <- ggplot() +
  geom_line(data = validation_curves,
            aes(R0, qlogis(probability_display), colour = theta_label), linewidth = .78) +
  geom_errorbar(data = sim_large,
                aes(R0, ymin = qlogis(low_display), ymax = qlogis(high_display),
                    colour = theta_label), width = 0, linewidth = .32, alpha = .75) +
  geom_point(data = sim_large,
             aes(R0, qlogis(estimate_display), colour = theta_label, shape = theta_label),
             size = 1.65, stroke = .4) +
  facet_wrap(~K, ncol = 2, labeller = labeller(K = function(x)
    paste0('K = ', format(as.numeric(x), scientific = TRUE)))) +
  scale_x_continuous(breaks = c(1.05, 2, 3, 4, 6)) +
  scale_y_continuous(breaks = qlogis(prob_breaks), labels = format(prob_breaks, trim = TRUE),
                     limits = qlogis(c(.001, .999)), expand = expansion(mult = c(.012, .012))) +
  scale_colour_manual(values = pal5, name = expression(theta), drop = FALSE) +
  scale_shape_manual(values = c(`0` = 16, `0.5` = 17, `1` = 15),
                     name = expression(theta), na.translate = FALSE) +
  labs(x = x_label_R0, y = 'Conditional persistence probability') +
  guides(colour = guide_legend(override.aes = list(linewidth = .85, size = 2)),
         shape = 'none') +
  theme_paper + theme(axis.title.x = element_text(family = 'Cambria Math'))
ggsave(file.path(out_dir, 'main_persistence_validation.pdf'), fig1,
       width = 7.15, height = 5.35, device = cairo_pdf)

# Figure 2: dense conditional-persistence landscape.
conditional_cache <- 'validation/data/paper_conditional_landscape.csv'
conditional <- read_cache(conditional_cache)
if (is.null(conditional)) {
  R_land <- seq(1.05, 8, length.out = 300)
  rho_land <- exp(seq(log(.001), log(.05), length.out = 170))
  core <- core_grid(R_land, theta3)
  conditional <- rbindlist(lapply(c(1e6, 1e9), function(Kv) {
    z <- cross_join_rho(core, rho_land)
    z <- add_probabilities(z, Kv)
    z <- add_entry_diagnostic(z, Kv)
    z[, `:=`(K = Kv, scan_version = scan_version)]
    z
  }))
  fwrite(conditional, conditional_cache)
}
conditional[, `:=`(
  theta_label = factor(sprintf('theta == %g', theta), levels = sprintf('theta == %g', theta3)),
  K_label = factor(sprintf('K == 10^%d', round(log10(K))),
                   levels = sprintf('K == 10^%d', c(6, 9))),
  probability_display = fifelse(finite_entry, P_conditional, NA_real_))]
fig2 <- ggplot(conditional, aes(R0, rho, fill = probability_display)) +
  geom_raster(interpolate = FALSE) +
  facet_grid(K_label ~ theta_label, labeller = label_parsed) +
  scale_y_log10(breaks = c(.001, .002, .005, .01, .02, .05),
                labels = c('.001', '.002', '.005', '.01', '.02', '.05')) +
  scale_x_continuous(breaks = 1:8, expand = c(0, 0)) +
  scale_fill_viridis_c(option = 'C', limits = c(0, 1), na.value = 'grey72',
                       name = 'Conditional\npersistence') +
  labs(x = x_label_R0, y = expression(rho)) + theme_paper +
  theme(axis.title.x = element_text(family = 'Cambria Math'),
        legend.position = 'right', panel.spacing = unit(.45, 'lines'))
ggsave(file.path(out_dir, 'conditional_persistence_landscape.pdf'), fig2,
       width = 7.15, height = 4.65, device = cairo_pdf)

# Figure 3: unconditional landscape and analytically continued extremum structure.
struct_R <- seq(1.05, 8, length.out = 360)
struct_core <- core_grid(struct_R, theta5)
rho_uncond <- exp(seq(log(.001), log(.04), length.out = 155))
unconditional_cache <- 'validation/data/paper_unconditional_landscape.csv'
unconditional <- read_cache(unconditional_cache)
if (is.null(unconditional)) {
  unconditional <- cross_join_rho(struct_core, rho_uncond)
  unconditional <- add_probabilities(unconditional, 1e6)
  unconditional <- add_entry_diagnostic(unconditional, 1e6)
  unconditional[, `:=`(K = 1e6, scan_version = scan_version)]
  fwrite(unconditional, unconditional_cache)
}
display_extrema <- rbindlist(lapply(theta5, function(th) {
  rbindlist(lapply(exp(seq(log(.0012), log(.03), length.out = 90)), function(rv) {
    ans <- as.data.table(find_extrema_analytic(rv, th, 1e6, Rmax = 40))
    if (nrow(ans)) {
      ans <- ans[R0_ext <= 8]
      ans[, `:=`(theta = th, rho = rv)]
    }
    ans
  }), fill = TRUE)
}), fill = TRUE)
unconditional[, `:=`(
  theta_label = factor(sprintf('theta == %g', theta), levels = sprintf('theta == %g', theta5)),
  probability_display = fifelse(finite_entry, P_unconditional, NA_real_))]
display_extrema[, theta_label := factor(sprintf('theta == %g', theta),
                                        levels = sprintf('theta == %g', theta5))]
landscape3 <- ggplot(unconditional, aes(R0, rho, fill = probability_display)) +
  geom_raster(interpolate = FALSE) +
  geom_line(data = display_extrema, aes(R0_ext, rho, linetype = type),
            inherit.aes = FALSE, colour = 'white', linewidth = .55) +
  facet_wrap(~theta_label, ncol = 3, labeller = label_parsed) +
  scale_y_log10(breaks = c(.001, .003, .01, .03), labels = c('.001', '.003', '.01', '.03')) +
  scale_x_continuous(breaks = c(1, 2, 4, 6, 8), expand = c(0, 0)) +
  scale_fill_viridis_c(option = 'C', limits = c(0, 1), na.value = 'grey72',
                       name = 'Unconditional\npersistence') +
  scale_linetype_manual(values = c('Local maximum' = 'dashed',
                                   'Local minimum' = 'solid'), name = NULL) +
  labs(x = x_label_R0, y = expression(rho), title = '(a) Unconditional persistence') +
  guides(linetype = guide_legend(override.aes = list(colour = 'black'))) +
  theme_paper + theme(axis.title.x = element_text(family = 'Cambria Math'),
                      legend.position = 'right', plot.title = element_text(size = 10, face = 'bold'),
                      panel.spacing = unit(.35, 'lines'))

branch_cache <- 'validation/data/paper_extrema_branches.csv'
branch_version <- 'one_sided_continuation_v2'
branches <- read_cache(branch_cache)
if (!is.null(branches) && (!'branch_version' %in% names(branches) ||
    !identical(unique(branches$branch_version), branch_version))) branches <- NULL
if (is.null(branches)) {
  critical_seed <- fread('validation/data/paper_theta_critical.csv')[K == 1e6]
  # Both branches move one-sidedly on these four rho slices. Bracketing from
  # the previous root tracks the same stationary point and stops the minimum
  # cleanly when it leaves the displayed R0 <= 8 window.
  continue_one_branch <- function(rv, th0, r0, type, step = .004) {
    th_grid <- unique(c(th0, seq(th0 + step, 1, by = step), 1))
    th_grid <- th_grid[th_grid <= 1]
    out <- data.table(theta = th0, R0_ext = r0, type = type, rho = rv)
    r_prev <- r0; th_prev <- th0
    for (th in th_grid[-1]) {
      fun <- function(r) unname(extremum_quantities(r, rv, th, 1e6)['G_scaled'])
      interval <- if (type == 'Local maximum') c(1.01, r_prev) else c(r_prev, 8)
      values <- vapply(interval, fun, 0.)
      if (!all(is.finite(values)) || prod(sign(values)) >= 0) {
        if (type == 'Local minimum' && is.finite(values[2])) {
          th_edge <- uniroot(function(t) unname(
            extremum_quantities(8, rv, t, 1e6)['G_scaled']),
            c(th_prev, th), tol = 2e-9)$root
          out <- rbind(out, data.table(theta = th_edge, R0_ext = 8,
                                       type = type, rho = rv))
        }
        break
      }
      r_new <- uniroot(fun, interval, tol = 2e-9)$root
      gr <- unname(extremum_quantities(r_new, rv, th, 1e6)['GR_scaled'])
      if ((type == 'Local maximum' && gr >= 0) ||
          (type == 'Local minimum' && gr <= 0)) stop('Branch type changed')
      out <- rbind(out, data.table(theta = th, R0_ext = r_new,
                                   type = type, rho = rv))
      r_prev <- r_new; th_prev <- th
    }
    out
  }
  branches <- rbindlist(lapply(c(.005, .01, .02, .03), function(rv) {
    if (rv <= .01) {
      th0 <- 0
      seed <- as.data.table(find_extrema_analytic(rv, th0, 1e6,
                                                   Rmax = 8, n_bracket = 600))
      stopifnot(nrow(seed) == 2)
      rmax <- seed[type == 'Local maximum', R0_ext]
      rmin <- seed[type == 'Local minimum', R0_ext]
    } else {
      saddle <- critical_seed[abs(rho - rv) < 1e-10]
      stopifnot(nrow(saddle) == 1)
      th0 <- saddle$theta_c; rmax <- rmin <- saddle$Rc
    }
    rbind(continue_one_branch(rv, th0, rmax, 'Local maximum'),
          continue_one_branch(rv, th0, rmin, 'Local minimum'))
  }), fill = TRUE)
  branches[, `:=`(scan_version = scan_version, branch_version = branch_version)]
  fwrite(branches, branch_cache)
}
branches[, rho_label := factor(sprintf('rho == %g', rho),
                                levels = sprintf('rho == %g', c(.005, .01, .02, .03)))]
branch_plot <- ggplot(branches, aes(theta, R0_ext, colour = type,
                                    linetype = type, group = type)) +
  geom_line(linewidth = .8) +
  facet_wrap(~rho_label, ncol = 2, labeller = label_parsed) +
  scale_colour_manual(values = c('Local maximum' = '#D55E00',
                                 'Local minimum' = '#0072B2'), name = NULL) +
  scale_linetype_manual(values = c('Local maximum' = 'dashed',
                                   'Local minimum' = 'solid'), name = NULL) +
  scale_x_continuous(breaks = c(0, .5, 1)) +
  scale_y_continuous(breaks = c(1, 2, 4, 6, 8), limits = c(1, 8),
                     oob = scales::oob_censor) +
  labs(x = expression(theta), y = paste0('Extremum location, ', x_label_R0),
       title = '(b) Continuous-theta extremum branches') + theme_paper +
  theme(legend.position = 'bottom', plot.title = element_text(size = 10, face = 'bold'),
        axis.title.y = element_text(family = 'Cambria Math'))

critical_cache <- 'validation/data/paper_theta_critical.csv'
critical_version <- 'direct_G_GR_v1'
critical <- if (file.exists(critical_cache)) fread(critical_cache) else NULL
if (!is.null(critical) && (!'critical_version' %in% names(critical) ||
    !identical(unique(critical$critical_version), critical_version))) critical <- NULL
if (is.null(critical)) {
  critical <- rbindlist(lapply(c(1e5, 1e6, 1e7), function(Kv) {
    rbindlist(lapply(seq(.014, .03, by = .002), function(rv) {
      ans <- as.data.table(solve_saddle_node(rv, Kv))
      if (!nrow(ans)) return(data.table())
      ans <- ans[theta_c >= 0 & theta_c <= 1]
      if (!nrow(ans)) return(data.table())
      ans <- ans[which.min(residual)]
      ans[, `:=`(K = Kv, rho = rv, scan_version = scan_version,
                 critical_version = critical_version)]
      ans
    }), fill = TRUE)
  }))
  fwrite(critical, critical_cache)
}
critical[, K_label := factor(sprintf('10^%d', round(log10(K))),
                             levels = sprintf('10^%d', 5:7))]
critical_plot <- ggplot(critical, aes(rho, theta_c, colour = K_label)) +
  geom_line(linewidth = .75) + geom_point(size = 1.25) +
  scale_colour_manual(values = c('#56B4E9', '#0072B2', '#004C6D'), name = expression(K)) +
  scale_x_continuous(breaks = c(.015, .02, .025, .03),
                     labels = c('.015', '.020', '.025', '.030')) +
  scale_y_continuous(limits = c(0, .5), breaks = seq(0, .5, .1)) +
  labs(x = expression(rho), y = expression(theta[c]), title = '(c) Saddle-node boundary') +
  theme_paper + theme(legend.position = 'bottom', plot.title = element_text(size = 10, face = 'bold'))

cairo_pdf(file.path(out_dir, 'unconditional_persistence_structure.pdf'),
          width = 7.15, height = 8.15)
grid::grid.newpage()
lay <- grid::grid.layout(2, 2, heights = grid::unit(c(1.45, 1), 'null'),
                         widths = grid::unit(c(1.75, 1), 'null'))
grid::pushViewport(grid::viewport(layout = lay))
print(landscape3, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2))
print(branch_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
print(critical_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

checks <- rbindlist(lapply(c(0, .25, .5), function(th) {
  row <- as.data.table(find_extrema_analytic(.01, th, 1e6, Rmax = 40))[
    type == 'Local minimum']
  if (!nrow(row)) return(NULL)
  r <- row$R0_ext[1]; m <- matching(r, .01, th, 1e6)
  kq <- kendall_quantities(m$x_second, r, .01, th, 1e6, yline = m$ybl)
  pbi <- bi_quantities(r, .01, th, 1e6)$P_conditional
  data.table(theta = th, rho = .01, K = 1e6, R0_ext = r,
             P_conditional_BI = pbi, P_conditional_exact = kq$P1,
             relative_difference = abs(kq$P1 / pbi - 1))
}), fill = TRUE)
fwrite(checks, 'validation/data/paper_extrema_exact_checks.csv')

# A separate check just above each confirmed saddle-node: evaluate the
# finite-entry, exact-Kendall unconditional probability on a local R0 mesh.
exact_branch_checks <- rbindlist(lapply(list(
  list(rho = .02, theta = .15, R0 = seq(1.6, 3.2, by = .08)),
  list(rho = .03, theta = .34, R0 = seq(2.1, 4.3, by = .1))
), function(p) {
  z <- data.table(rho = p$rho, theta = p$theta, K = 1e6, R0 = p$R0)
  z[, P_exact_unconditional := vapply(R0, function(r) {
    m <- matching(r, p$rho, p$theta, 1e6)
    if (!is.finite(m$x_second)) return(NA_real_)
    (1 - 1 / r) * kendall_quantities(m$x_second, r, p$rho,
                                     p$theta, 1e6, yline = m$ybl)$P1
  }, numeric(1))]
  z
}))
fwrite(exact_branch_checks, 'validation/data/paper_exact_kendall_branch_checks.csv')

# Figure 4: stochastic error and validity map at K = 30,000.
sim_standard <- fread('validation/data/stochastic_scan_results.csv')[K == 30000]
sim_standard[, prediction := vapply(seq_len(.N), function(i)
  bi_quantities(R0[i], rho[i], theta[i], K[i])$P_conditional, numeric(1))]
sim_standard[, abs_error := abs(prediction - P_conditional)]
sim_standard[, error_display := pmax(abs_error, 1e-3)]
sim_standard[, no_entry := vapply(seq_len(.N), function(i) {
  m <- matching(R0[i], rho[i], theta[i], K[i], D = 0); !is.finite(m$x_first)
}, logical(1))]
sim_standard[, theta_label := factor(sprintf('theta == %g', theta),
                                     levels = sprintf('theta == %g', theta3))]
fig4 <- ggplot(sim_standard, aes(R0, rho)) +
  geom_point(aes(fill = error_display), shape = 22, size = 4.5,
             colour = 'grey25', stroke = .25) +
  geom_point(data = sim_standard[no_entry == TRUE], shape = 4,
             size = 1.45, colour = 'grey65', stroke = .42) +
  facet_wrap(~theta_label, nrow = 1, labeller = label_parsed) +
  scale_y_log10(breaks = c(.01, .02, .05, .1), labels = c('.01', '.02', '.05', '.10')) +
  scale_x_continuous(breaks = c(1.05, 2, 3, 4, 6)) +
  scale_fill_viridis_c(option = 'B', trans = 'log10', limits = c(1e-3, .6),
                       breaks = c(1e-3, .01, .05, .1, .3, .6),
                       labels = c('<.001', '.01', '.05', '.10', '.30', '.60'),
                       name = 'Absolute\nprobability error') +
  labs(x = x_label_R0, y = expression(rho)) + theme_paper +
  theme(axis.title.x = element_text(family = 'Cambria Math'), legend.position = 'right')
ggsave(file.path(out_dir, 'stochastic_error_validity_map.pdf'), fig4,
       width = 7.15, height = 3.15, device = cairo_pdf)

# Supplementary approximation-layer diagnostics.
ode <- fread('validation/data/K_R0_scan.csv')[boundary_choice == 'sqrt']
ode[, BI := vapply(seq_len(.N), function(i)
  bi_quantities(R0[i], rho[i], theta[i], K[i])$P_conditional, numeric(1))]
ode_ok <- ode[status == 'OK' & is.finite(P1_ref)]
diag_ode <- melt(ode_ok[, .(`BI closed form` = mean(abs(BI - P1_ref), na.rm = TRUE),
  `Laplace step` = mean(abs(P1_L_ODE - P1_ref), na.rm = TRUE),
  `First-order entry` = mean(abs(P1_K_first - P1_ref), na.rm = TRUE),
  `Second-order entry` = mean(abs(P1_K_second - P1_ref), na.rm = TRUE)), by = rho],
  id.vars = 'rho', variable.name = 'component', value.name = 'MAE')
sim_all <- rbindlist(list(fread('validation/data/stochastic_scan_results.csv'),
                          fread('validation/data/fig11_scan_results.csv')), fill = TRUE)
sim_all[, prediction := vapply(seq_len(.N), function(i)
  bi_quantities(R0[i], rho[i], theta[i], K[i])$P_conditional, numeric(1))]
diag_sim <- sim_all[, .(MAE = mean(abs(prediction - P_conditional))), by = rho]
figS2a <- ggplot(diag_sim, aes(rho, MAE)) +
  geom_line(linewidth = .8, colour = '#7A0177') + geom_point(size = 2, colour = '#7A0177') +
  labs(x = expression(rho), y = 'Mean absolute probability error',
       title = '(a) Complete stochastic comparison') + theme_paper +
  theme(legend.position = 'none', plot.title = element_text(size = 10))
figS2b <- ggplot(diag_ode, aes(rho, MAE, colour = component, shape = component)) +
  geom_line(linewidth = .7) + geom_point(size = 1.7) +
  scale_colour_manual(values = c('#CC79A7', '#D55E00', '#56B4E9', '#009E73')) +
  labs(x = expression(rho), y = 'Mean absolute probability difference',
       colour = NULL, shape = NULL, title = '(b) ODE/Kendall error isolation') +
  theme_paper + theme(plot.title = element_text(size = 10), legend.position = 'right',
                      legend.text = element_text(size = 7))
cairo_pdf(file.path(out_dir, 'accuracy_diagnostics_supplement.pdf'),
          width = 7.15, height = 3.45)
grid::grid.newpage(); lay2 <- grid::grid.layout(1, 2)
grid::pushViewport(grid::viewport(layout = lay2))
print(figS2a, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(figS2b, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

summary_rows <- rbindlist(list(
  sim_all[, .(quantity = 'simulation_MAE', value = mean(abs(prediction - P_conditional))), by = rho],
  sim_all[rho <= .02 & R0 >= 1.15,
          .(quantity = 'simulation_MAE_separated',
            value = mean(abs(prediction - P_conditional))), by = .(rho, theta)],
  ode_ok[, .(quantity = 'second_order_entry_MAE',
             value = mean(abs(P1_K_second - P1_ref), na.rm = TRUE)), by = rho]
), fill = TRUE)
fwrite(summary_rows, 'validation/data/paper_numerical_summary.csv')

# Independent finite-difference diagnostics for the quadrature-level analytic
# derivatives.  These checks are not used in any production root calculation.
derivative_checks <- rbindlist(lapply(c(1.5, 2.5, 5), function(rv)
  rbindlist(lapply(c(0, .5, 1), function(th) {
    h <- 1e-3
    cd <- C_log_derivatives(rv, th)
    f <- function(r) log(C_explicit(r, th))
    data.table(R0 = rv, theta = th,
      logC_relative_error = exp(cd['logC']) / C_explicit(rv, th) - 1,
      d1_analytic = cd['d1'], d1_finite_difference = (f(rv + h) - f(rv - h)) / (2 * h),
      d2_analytic = cd['d2'],
      d2_finite_difference = (f(rv + h) - 2 * f(rv) + f(rv - h)) / h^2)
  }))))
fwrite(derivative_checks, 'validation/data/paper_extremum_derivative_checks.csv')
message('Section 7 figures and cached scans written successfully')

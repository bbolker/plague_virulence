source('validation/R/theory.R')
source(here::here('graphics_utils.R'))
source('validation/R/kendall.R')
source('validation/R/extrema.R')

stopifnot(requireNamespace('data.table', quietly = TRUE),
          requireNamespace('ggplot2', quietly = TRUE),
          requireNamespace('scales', quietly = TRUE))
library(data.table)
library(ggplot2)

out_dir <- 'validation/figures/paper'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
scan_version <- 'section7_landscape_R32_stable_v4'
theta5 <- seq(0, 1, by = .25)
theta3 <- c(0, .5, 1)
pal5 <- c('#0072B2', '#56B4E9', '#009E73', '#E69F00', '#D55E00')
x_label_R0 <- paste0(intToUtf8(0x211B), intToUtf8(0x2080))

## Output backend for the landscape figures: tikz when requested (--tikz
## argument or FIG_BACKEND=tikz), Quartz on macOS, otherwise cairo PDF.
## Quartz PDF output inverts rasters, so only the Quartz backend draws the
## landscapes as vector cells. tikz output must use rasters: thousands of
## vector cells make LaTeX very slow or exhaust its memory.
fig_backend <- if ('--tikz' %in% commandArgs(trailingOnly = TRUE) ||
                   Sys.getenv('FIG_BACKEND') == 'tikz') {
  'tikz'
} else if (capabilities('aqua')) 'quartz' else 'cairo'
use_tikz <- fig_backend == 'tikz'
## Landscape labels: LaTeX strings under tikz, plotmath otherwise.
theta_lab <- function(th) {
  if (use_tikz) sprintf('$\\theta = %g$', th) else sprintf('theta == %g', th)
}
K_lab <- function(K) {
  e <- round(log10(K))
  if (use_tikz) sprintf('$K = 10^{%d}$', e) else sprintf('K == 10^%d', e)
}
landscape_labeller <- if (use_tikz) label_value else label_parsed
landscape_xlab <- if (use_tikz) '$\\mathcal{R}_0$' else bquote(.(intToUtf8(0x211B))[0])
landscape_ylab <- if (use_tikz) '$\\rho$' else expression(rho)
## Write a landscape figure to out_dir/<name>.{pdf,tex} with the chosen backend.
save_landscape <- function(p, name, width = 10.5, height = 4.65) {
  switch(fig_backend,
    tikz = {
      tex_file <- file.path(out_dir, paste0(name, '.tex'))
      tikzDevice::tikz(tex_file, width = width, height = height)
      print(p + theme(axis.title.x = element_text(family = '')))
      dev.off()
      ## tikzDevice refers to its raster PNGs by bare file name, and \pgfimage
      ## ignores \graphicspath, so prefix the path relative to the directory
      ## the paper is compiled from (the same one these scripts run from).
      tex <- readLines(tex_file)
      writeLines(gsub(paste0('{', name, '_ras'),
                      paste0('{', out_dir, '/', name, '_ras'), tex, fixed = TRUE),
                 tex_file)
    },
    quartz = {
      quartzFonts(STIX = quartzFont(c('STIXGeneral-Regular', 'STIXGeneral-Bold',
                                     'STIXGeneral-Italic', 'STIXGeneral-BoldItalic')))
      ## Vector cells need a continuous gradient in the colour bar, too.
      p <- p + guides(fill = guide_colourbar(display = 'gradient', order = 1)) +
        theme(axis.title.x = element_text(family = 'STIX'))
      quartz(file = file.path(out_dir, paste0(name, '.pdf')),
             type = 'pdf', width = width, height = height)
      print(p)
      dev.off()
    },
    cairo = ggsave(file.path(out_dir, paste0(name, '.pdf')), p,
                   width = width, height = height, device = cairo_pdf))
  invisible(NULL)
}

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
  logC <- unname(C_log_derivatives(R0, theta)['logC']); C <- exp(logC); delta <- 1 - R0 * xf
  data.table(R0 = R0, theta = theta, xf = xf, delta = delta,
             hf = h_theta(xf, theta), C = C, logC = logC,
             hstar = h_theta(xs, theta),
             DeltaA = unname(action_barrier_derivatives(R0, theta)['value']))
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

# Figure 2: dense conditional-persistence landscape.
conditional_cache <- 'validation/data/paper_conditional_landscape_theta5_R32.csv'
conditional <- read_cache(conditional_cache)
if (is.null(conditional)) {
  R_land <- 1 + exp(seq(log(.05), log(31), length.out = 320))
  rho_land <- exp(seq(log(.001), log(.05), length.out = 130))
  core <- core_grid(R_land, theta5)
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
  theta_label = factor(theta_lab(theta), levels = theta_lab(theta5)),
  K_label = factor(K_lab(K), levels = K_lab(c(1e6, 1e9))),
  probability_display = pmax(1e-6, P_unconditional))]
land_dlogR <- min(diff(log(sort(unique(conditional$R0))-1)))
land_dlogrho <- min(diff(log(sort(unique(conditional$rho)))))
# The grid is evenly spaced in log(R0 - 1) and log(rho), so it is a regular
# raster on the log axes. Quartz PDF output inverts rasters, so the Quartz
# backend draws vector cells instead.
landscape_cells <- if (fig_backend == 'quartz') {
  geom_rect(aes(xmin = 1 + pmax(.05, (R0 - 1) * exp(-land_dlogR / 2)),
                xmax = 1 + pmin(31, (R0 - 1) * exp(land_dlogR / 2)),
                ymin = rho * exp(-land_dlogrho / 2),
                ymax = rho * exp(land_dlogrho / 2)))
} else {
  geom_raster(aes(R0, rho), interpolate = FALSE)
}
fig2 <- ggplot(conditional, aes(fill = probability_display)) +
  landscape_cells +
  facet_grid(K_label ~ theta_label, labeller = landscape_labeller) +
  scale_y_log10(breaks = c(.001, .002, .005, .01, .02, .05),
                labels = c('0.001', '0.002', '0.005', '0.01', '0.02', '0.05'),
                expand = c(0, 0)) +
  scale_x_log10_shifted(breaks = c(1.1, 1.3, 2, 4, 10, 30), expand = c(0, 0)) +
  scale_fill_viridis_c(option = 'C', limits = c(1e-6, 1),
                       transform = scales::transform_log10(),
                       breaks = c(1e-6, 1e-5, 1e-4, .001, .01, .1, 1),
                       labels = if (use_tikz) {
                         c('$\\le 10^{-6}$', sprintf('$10^{%d}$', -5:-1), '1')
                       } else {
                         expression('' <= 10^{-6}, 10^{-5}, 10^{-4}, 10^{-3}, 10^{-2}, 10^{-1}, 1)
                       },
                       name = 'Unconditional\npersistence\n(log scale)') +
  labs(x = landscape_xlab, y = landscape_ylab) + theme_paper +
  theme(axis.title.x = element_text(family = 'Cambria Math'),
        axis.title.y = element_text(angle = 0, vjust = .5),
        strip.text.y.right = element_text(angle = 0),
        legend.position = 'right', panel.spacing = unit(.45, 'lines'))
ext_path <- 'validation/data/unconditional_landscape_extrema_theta5_K2_R32.csv'
if (file.exists(ext_path)) {
  extrema_display <- fread(ext_path)
} else {
  extrema_display <- rbindlist(lapply(c(1e6, 1e9), function(kv)
    rbindlist(lapply(theta5, function(th)
      rbindlist(lapply(exp(seq(log(.001), log(.05), length.out = 180)), function(rv) {
        ans <- as.data.table(find_extrema_analytic(rv, th, kv, Rmax = 32, n_bracket = 500))
        if (!nrow(ans)) return(NULL)
        ans <- ans[R0_ext >= 1.05 & R0_ext <= 32]
        ans[, `:=`(theta = th, K = kv, rho = rv)]
        ans
      }), fill = TRUE)), fill = TRUE)), fill = TRUE)
  fwrite(extrema_display, ext_path)
}
extrema_display[, `:=`(
  theta_label = factor(theta_lab(theta), levels = theta_lab(theta5)),
  K_label = factor(K_lab(K), levels = K_lab(c(1e6, 1e9))))]
fig2 <- fig2 +
  geom_path(data = extrema_display[order(K, theta, type, rho)],
            aes(R0_ext, rho, linetype = type,
                group = interaction(K, theta, type)),
            inherit.aes = FALSE, colour = 'grey50', linewidth = .55) +
  scale_linetype_manual(values = c('Local maximum' = 'dashed', 'Local minimum' = 'solid'), name = NULL) +
  guides(fill = guide_colourbar(order = 1),
         linetype = guide_legend(order = 2, override.aes = list(colour = 'grey50')))

save_landscape(fig2, 'unconditional_persistence_landscape')

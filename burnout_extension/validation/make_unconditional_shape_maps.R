source('validation/R/theory.R')
source('validation/R/extrema.R')
stopifnot(requireNamespace('data.table', quietly = TRUE),
          requireNamespace('ggplot2', quietly = TRUE))
library(data.table)
library(ggplot2)

out_dir <- 'validation/figures/paper'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rho_values <- seq(.002, .04, length.out = 65)
theta_values <- seq(0, 1, by = .025)
K_values <- c(1e5, 1e6, 1e7)
R_grid <- unique(c(seq(1.001, 10, length.out = 900),
                   exp(seq(log(10), log(300), length.out = 220))))
classified_cache <- 'validation/data/paper_unconditional_shape_map.csv'
if (!file.exists(classified_cache)) {

# Cache derivatives of the analytical G condition's R0-only constituents.
base <- rbindlist(lapply(theta_values, function(th) {
  z <- lapply(R_grid, function(r) {
    cd <- C_log_derivatives(r, th)
    ad <- action_barrier_derivatives(r, th)
    c(logC = unname(cd['logC']), logC_d1 = unname(cd['d1']),
      DeltaA = unname(ad['value']), DeltaA_d1 = unname(ad['d1']))
  })
  q <- as.data.table(do.call(rbind, z))
  q[, `:=`(R0 = R_grid, theta = th)]
  q
}))

g_sign <- function(q, rho, K) {
  r <- q$R0; th <- q$theta[1]
  L <- q$logC_d1 + 1 / (2 * (r - 1)) - th / (2 * r) - q$DeltaA_d1 / rho
  logB <- log(K) + q$logC + .5 * (log(rho) + log(r - 1) -
    th * log(r) - log(2 * pi)) - q$DeltaA / rho
  B <- exp(pmin(logB, log(40)))
  exprel <- ifelse(B < 1e-5, 1 + B / 2 + B^2 / 6, expm1(B) / B)
  value <- L + exprel / (r * (r - 1))
  value[logB >= log(40)] <- 1
  sign(value)
}

# For theta < 1 the far tail has C_theta ~ exp(-gamma_E)/[(1-theta)R0].
# Solve R0 G/B=0 in log R0; multiplying by R0 avoids underflow at remote roots.
remote_maximum <- function(rho, theta, K) {
  delta <- 1 - theta
  stopifnot(delta > 0)
  f <- function(y) {
    r <- exp(y)
    powers <- seq.int(0L, 7L) + 2 - theta
    rL <- -1 + r / (2 * (r - 1)) - theta / 2 +
      sum(exp(-(powers - 1) * y) / powers) / rho
    action <- exp(-delta * y) / (delta * (1 + delta))
    logB <- log(K) - .5772156649015329 - log(delta) - y +
      .5 * (log(rho) + log(r - 1) - theta * y - log(2 * pi)) - action / rho
    if (logB > log(40)) return(1)
    B <- exp(logB)
    exprel <- if (B < 1e-5) 1 + B / 2 + B^2 / 6 else expm1(B) / B
    rL + exprel / (r - 1)
  }
  lo <- log(300); hi <- 300
  if (!(f(lo) > 0 && f(hi) < 0)) stop('Remote-root bracket failed')
  exp(uniroot(f, c(lo, hi), tol = 1e-8)$root)
}

classify_one <- function(q, rho, theta, K) {
  signs <- g_sign(q, rho, K)
  stopifnot(all(is.finite(signs)), signs[1] > 0)
  ii <- which(signs[-length(signs)] != signs[-1])
  roots <- if (length(ii)) vapply(ii, function(i) uniroot(function(r)
    unname(extremum_quantities(r, rho, theta, K)['G_scaled']),
    q$R0[c(i, i + 1)], tol = 2e-8)$root, 0.) else numeric()
  types <- if (length(roots)) vapply(roots, function(r) {
    gr <- unname(extremum_quantities(r, rho, theta, K)['GR_scaled'])
    if (gr < 0) 'maximum' else 'minimum'
  }, '') else character()
  if (theta == 1) {
    stopifnot(length(roots) == 1L, identical(types, 'maximum'), tail(signs, 1) < 0)
    remote <- NA_real_
  } else {
    stopifnot(length(roots) %in% c(0L, 2L), tail(signs, 1) > 0,
              identical(types, rep(c('maximum', 'minimum'), length(roots) / 2)))
    remote <- remote_maximum(rho, theta, K)
  }
  nmax <- sum(types == 'maximum') + as.integer(theta < 1)
  nmin <- sum(types == 'minimum')
  category <- if (nmax == 0 && nmin == 0) 'No extrema' else
    if (nmax == 1 && nmin == 0) 'One maximum' else
    if (nmax == 0 && nmin == 1) 'One minimum' else
    if (nmax == 1 && nmin == 1) 'One maximum + one minimum' else
    if (nmax == 2 && nmin == 1) 'Two maxima + one minimum' else
    stop('Unexpected extremum structure')
  data.table(rho = rho, theta = theta, K = K, category = category,
             n_max = nmax, n_min = nmin,
             first_root = if (length(roots)) roots[1] else NA_real_,
             minimum_root = if (length(roots) == 2) roots[2] else NA_real_,
             remote_maximum = remote)
}

classified <- rbindlist(lapply(K_values, function(k) {
  rbindlist(lapply(theta_values, function(th) {
    q <- base[theta == th]
    rbindlist(lapply(rho_values, function(rv) classify_one(q, rv, th, k)))
  }))
}))
fwrite(classified, classified_cache)
}
classified <- fread(classified_cache)
classified[, remote_maximum := as.numeric(remote_maximum)]
classified[category == 'Maximum only', category := 'One maximum']
classified[category == 'Both maximum and minimum',
           category := 'Two maxima + one minimum']
fwrite(classified, classified_cache)
print(classified[, .N, by = .(K, category)][order(K, category)])

# Continue the verified saddle-node solutions rather than interpolating
# disconnected root-search points. Curves are omitted where theta_c < 0.
critical_seed <- fread('validation/data/paper_theta_critical.csv')
critical <- rbindlist(lapply(K_values, function(k) {
  rbindlist(lapply(seq(.014, .04, by = .001), function(rv) {
    seed <- critical_seed[K == k][which.min(abs(rho - rv))]
    ans <- solve_saddle_node(rv, k, starts = data.frame(R0 = seed$Rc,
                                                        theta = seed$theta_c))
    if (!nrow(ans)) return(data.table())
    ans <- as.data.table(ans[which.min(ans$residual), ])
    if (ans$theta_c < 0 || ans$theta_c > 1) return(data.table())
    ans[, `:=`(rho = rv, K = k)]
    ans
  }), fill = TRUE)
}), fill = TRUE)
fwrite(critical, 'validation/data/paper_unconditional_shape_boundary.csv')

levels5 <- c('No extrema', 'One maximum', 'One minimum',
             'One maximum + one minimum', 'Two maxima + one minimum')
colors5 <- c('No extrema' = '#E7E7E7', 'One maximum' = '#4477AA',
             'One minimum' = '#DDAA33',
             'One maximum + one minimum' = '#AA4499',
             'Two maxima + one minimum' = '#228866')
classified[, category := factor(category, levels = levels5)]
theme_map <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        strip.background = element_rect(fill = 'grey94', colour = 'grey55'),
        legend.position = 'bottom', legend.title = element_blank())
map_plot <- function(data, boundary) {
  ggplot(data, aes(rho, theta, fill = category)) +
    geom_raster(interpolate = FALSE) +
    geom_line(data = boundary, aes(rho, theta_c), inherit.aes = FALSE,
              colour = 'white', linewidth = 1.5) +
    geom_line(data = boundary, aes(rho, theta_c), inherit.aes = FALSE,
              colour = 'black', linetype = 'dashed', linewidth = .65) +
    scale_fill_manual(values = colors5, drop = TRUE) +
    scale_x_continuous(breaks = c(.005, .01, .02, .03, .04),
                       labels = c('.005', '.010', '.020', '.030', '.040'),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, by = .2), expand = c(0, 0)) +
    labs(x = expression(rho), y = expression(theta)) + theme_map
}

main <- map_plot(classified[K == 1e6], critical[K == 1e6])
ggsave(file.path(out_dir, 'unconditional_shape_map.pdf'), main,
       width = 6.6, height = 4.65, device = cairo_pdf)

comparison <- map_plot(classified, critical) +
  facet_wrap(~K, nrow = 1, labeller = labeller(K = function(x)
    paste0('K = ', format(as.numeric(x), scientific = TRUE)))) +
  theme(panel.spacing = grid::unit(.35, 'lines'))
ggsave(file.path(out_dir, 'unconditional_shape_map_K_comparison.pdf'), comparison,
       width = 8.2, height = 3.55, device = cairo_pdf)

# Analytical extremum conditions for the boundary-layer-independent
# conditional and unconditional persistence approximations. Derivatives of C_theta are
# propagated through its regularized one-dimensional quadrature with
# second-order forward automatic differentiation; no persistence curve or
# action-barrier finite differences are used in production calculations.

stopifnot(requireNamespace('statmod', quietly = TRUE))

.jet <- function(v, d1 = 0, d2 = 0) list(v = v, d1 = d1, d2 = d2)
.jadd <- function(a, b) .jet(a$v + b$v, a$d1 + b$d1, a$d2 + b$d2)
.jneg <- function(a) .jet(-a$v, -a$d1, -a$d2)
.jsub <- function(a, b) .jadd(a, .jneg(b))
.jmul <- function(a, b) .jet(a$v * b$v,
  a$d1 * b$v + a$v * b$d1,
  a$d2 * b$v + 2 * a$d1 * b$d1 + a$v * b$d2)
.jinv <- function(a) .jet(1 / a$v, -a$d1 / a$v^2,
  2 * a$d1^2 / a$v^3 - a$d2 / a$v^2)
.jdiv <- function(a, b) .jmul(a, .jinv(b))
.jlog <- function(a) .jet(log(a$v), a$d1 / a$v,
  a$d2 / a$v - (a$d1 / a$v)^2)
.jexp <- function(a) {
  e <- exp(a$v)
  .jet(e, e * a$d1, e * (a$d2 + a$d1^2))
}
.jexpm1 <- function(a) {
  e <- exp(a$v)
  .jet(expm1(a$v), e * a$d1, e * (a$d2 + a$d1^2))
}
.jpow <- function(a, p) .jet(a$v^p, p * a$v^(p - 1) * a$d1,
  p * (p - 1) * a$v^(p - 2) * a$d1^2 + p * a$v^(p - 1) * a$d2)
.jconst <- function(x) .jet(x, 0, 0)

x_final_derivatives <- function(R0) {
  xf <- x_final(R0); den <- 1 - R0 * xf
  d1 <- xf * (xf - 1) / den
  num <- xf * (xf - 1)
  num1 <- d1 * (2 * xf - 1)
  den1 <- -xf - R0 * d1
  d2 <- (num1 * den - num * den1) / den^2
  c(xf = xf, d1 = d1, d2 = d2)
}

C_log_derivatives <- local({
  cache <- new.env(parent = emptyenv())
  quad <- statmod::gauss.quad(192, kind = 'legendre')
  tt <- (quad$nodes + 1) / 2
  ww <- quad$weights / 2
  function(R0, theta) {
    key <- sprintf('%.13g|%.13g', R0, theta)
    if (exists(key, cache, inherits = FALSE)) return(get(key, cache))
    fd <- x_final_derivatives(R0)
    Rj <- .jet(R0, 1, 0)
    fj <- .jet(unname(fd['xf']), unname(fd['d1']), unname(fd['d2']))
    one <- .jconst(1)
    xsj <- .jinv(Rj)
    span <- .jsub(one, fj)
    smax <- .jneg(.jlog(fj))
    sj <- .jmul(smax, .jconst(tt))
    em1 <- .jexpm1(sj)
    es <- .jadd(one, em1)
    uj <- .jmul(fj, es)
    Fj <- .jadd(.jneg(.jmul(fj, em1)), .jdiv(sj, Rj))
    # Cancel powers of x_f symbolically before evaluating the quadrature:
    # (h(u)/u)/a_tilde = (1-u) exp((theta-1)s)/(1-x_f).
    # The uncancelled form overflows/underflows at large R0.
    decay <- .jexp(.jmul(.jconst(theta - 1), sj))
    first <- .jdiv(.jmul(.jmul(.jsub(xsj, uj), .jsub(one, uj)), decay),
                   .jmul(.jsub(one, fj), Fj))
    second <- .jdiv(es, em1)
    qj <- .jsub(first, second)
    intj <- .jet(sum(ww * qj$v), sum(ww * qj$d1), sum(ww * qj$d2))
    J_over_at <- .jmul(smax, intj)
    logCj <- .jadd(.jsub(.jadd(.jlog(span), .jlog(.jsub(xsj, fj))), .jlog(fj)),
                   J_over_at)
    ans <- c(logC = logCj$v, d1 = logCj$d1, d2 = logCj$d2,
             J_over_at = J_over_at$v, xf = fd['xf'])
    assign(key, ans, cache); ans
  }
})

action_barrier_derivatives <- local({
 cache <- new.env(parent = emptyenv())
 function(R0, theta) {
  key <- sprintf('%.13g|%.13g', R0, theta)
  if (exists(key, cache, inherits = FALSE)) return(get(key, cache))
  xf <- x_final(R0); xs <- 1 / R0
  upper <- log(xs / xf)
  value <- integrate(function(t) {
    x <- xf * exp(t); (1 - R0 * x) * x^(1 - theta) / (1 - x)
  }, 0, upper, rel.tol = 2e-11, abs.tol = 2e-12,
  subdivisions = 700L, stop.on.error = TRUE)$value
  integ <- integrate(function(t) {
    x <- xf * exp(t); x^(2 - theta) / (1 - x)
  }, 0, upper, rel.tol = 2e-11, abs.tol = 2e-12,
  subdivisions = 700L, stop.on.error = TRUE)$value
  d1 <- xf^(1 - theta) - integ
  d2 <- R0^(theta - 2) / (R0 - 1) -
    xf^(1 - theta) * ((1 - theta) + theta * xf) / (1 - R0 * xf)
  ans <- c(value = value, d1 = d1, d2 = d2)
  assign(key, ans, cache); ans
 }
})

extremum_quantities <- function(R0, rho, theta, K) {
  cd <- C_log_derivatives(R0, theta)
  ad <- action_barrier_derivatives(R0, theta)
  logC <- unname(cd['logC']); logC_d1 <- unname(cd['d1'])
  logC_d2 <- unname(cd['d2']); DeltaA <- unname(ad['value'])
  DeltaA_d1 <- unname(ad['d1']); DeltaA_d2 <- unname(ad['d2'])
  logB <- log(K) + logC + .5 * (log(rho) + log(R0 - 1) -
    theta * log(R0) - log(2 * pi)) - DeltaA / rho
  B <- exp(logB)
  L <- logC_d1 + 1 / (2 * (R0 - 1)) - theta / (2 * R0) - DeltaA_d1 / rho
  Lp <- logC_d2 - 1 / (2 * (R0 - 1)^2) + theta / (2 * R0^2) - DeltaA_d2 / rho
  exprel <- if (B < 1e-5) 1 + B / 2 + B^2 / 6 else expm1(B) / B
  G_scaled <- L + exprel / (R0 * (R0 - 1))
  GR_scaled <- L^2 + Lp - (2 * R0 - 1) * exprel /
    (R0^2 * (R0 - 1)^2) + exp(B) * L / (R0 * (R0 - 1))
  G <- if (B < 700) B * G_scaled else sign(G_scaled) * Inf
  GR <- if (B < 700) B * GR_scaled else sign(GR_scaled) * Inf
  c(B = B, logB = logB, L = L, Lp = Lp, G = G, GR = GR,
    G_scaled = G_scaled, GR_scaled = GR_scaled,
    logC = logC, logC_d1 = logC_d1, logC_d2 = logC_d2,
    DeltaA = DeltaA, DeltaA_d1 = DeltaA_d1, DeltaA_d2 = DeltaA_d2)
}

# Fixed-theta large-R0 continuation of the conditional L equation. The
# exponentially small final-size term and amplitude correction are omitted;
# this is used only for remote roots beyond double-precision quadrature.
conditional_L_largeR <- function(R0, rho, theta, n_terms = 8L) {
  stopifnot(R0 > 1, rho > 0, theta >= 0, theta < 1)
  powers <- seq.int(0L, n_terms - 1L) + 2 - theta
  integral <- sum(R0^(-powers) / powers)
  -1 / R0 + 1 / (2 * (R0 - 1)) - theta / (2 * R0) + integral / rho
}

find_extrema_analytic <- function(rho, theta, K, Rmax = 20, n_bracket = 240) {
  rr <- exp(seq(log(1.01), log(Rmax), length.out = n_bracket))
  gg <- vapply(rr, function(r) extremum_quantities(r, rho, theta, K)['G_scaled'], 0.)
  ii <- which(is.finite(gg[-length(gg)]) & is.finite(gg[-1]) &
                gg[-length(gg)] * gg[-1] < 0)
  if (!length(ii)) return(data.frame())
  roots <- vapply(ii, function(i) uniroot(function(r)
    extremum_quantities(r, rho, theta, K)['G_scaled'], rr[c(i, i + 1)],
    tol = 2e-9)$root, 0.)
  gr <- vapply(roots, function(r) extremum_quantities(r, rho, theta, K)['GR_scaled'], 0.)
  data.frame(R0_ext = roots,
             type = ifelse(gr < 0, 'Local maximum', 'Local minimum'),
             GR_scaled = gr)
}

solve_saddle_node <- function(rho, K, starts = NULL) {
  stopifnot(requireNamespace('nleqslv', quietly = TRUE))
  if (is.null(starts)) starts <- expand.grid(R0 = c(1.5, 2, 2.3, 3, 4),
                                             theta = c(.1, .15, .25, .4, .7))
  fun <- function(z) {
    if (any(!is.finite(z)) || z[1] <= 1.0001 || z[1] > 25 ||
        z[2] < -1 || z[2] > 2) return(c(1e3, 1e3))
    q <- try(extremum_quantities(z[1], rho, z[2], K), silent = TRUE)
    if (inherits(q, 'try-error') || any(!is.finite(q[c('G_scaled', 'GR_scaled')])))
      return(c(1e3, 1e3))
    unname(q[c('G_scaled', 'GR_scaled')])
  }
  sols <- lapply(seq_len(nrow(starts)), function(i) {
    fit <- try(nleqslv::nleqslv(as.numeric(starts[i, ]), fun, method = 'Newton',
                                control = list(ftol = 1e-10, xtol = 1e-10,
                                               maxit = 150)), silent = TRUE)
    if (inherits(fit, 'try-error')) return(NULL)
    res <- fun(fit$x)
    if (fit$x[1] <= 1 || fit$x[2] < -0.5 || fit$x[2] > 1.5 ||
        max(abs(res)) > 2e-6) return(NULL)
    data.frame(Rc = fit$x[1], theta_c = fit$x[2],
               residual = max(abs(res)), termcd = fit$termcd)
  })
  sols <- do.call(rbind, sols)
  if (is.null(sols) || !nrow(sols)) return(data.frame())
  sols <- sols[order(sols$residual), ]
  sols[!duplicated(round(sols$theta_c, 6)), , drop = FALSE]
}

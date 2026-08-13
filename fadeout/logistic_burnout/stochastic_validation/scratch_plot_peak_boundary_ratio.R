## Scratch: base-R plot of y_peak/y_BL vs R0 at r=0.1, from just above the
## overdamped threshold R0_crit through R0=5. Companion to
## scratch_peak_boundary_ratio_vs_R0.R; not part of the committed pipeline.
out_dir <- file.path("fadeout", "logistic_burnout", "stochastic_validation", "outputs")
fig_dir <- file.path("fadeout", "logistic_burnout", "stochastic_validation", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

res <- read.csv(file.path(out_dir, "scratch_peak_boundary_ratio_vs_R0.csv"))
res <- res[res$status == "success", ]

r_fixed <- 0.1
R0_crit <- 2 * (1 - sqrt(1 - r_fixed)) / r_fixed

x <- res$R0 - 1
y <- res$ratio
flagged <- c(1.1, 2.26426, 3)
is_flag <- sapply(res$R0, function(v) any(abs(v - flagged) < 1e-3))

png(file.path(fig_dir, "scratch_peak_boundary_ratio_vs_R0.png"),
    width = 1700, height = 1150, res = 200)
par(mar = c(4.8, 4.8, 3.8, 1.2))
plot(x, y, log = "xy", type = "n",
     xlab = "R0 - 1 (log scale)", ylab = "y_peak / y_BL (log scale)",
     main = "Peak / boundary-layer ratio vs R0 at fixed r = 0.1 (logistic model)")
abline(h = 1, lty = 3, col = "grey50")
abline(v = R0_crit - 1, lty = 2, col = "steelblue", lwd = 1.5)
lines(x, y, lwd = 2.2)
points(x, y, pch = 16, cex = 0.9)
points(x[is_flag], y[is_flag], pch = 16, cex = 1.8, col = "firebrick")
dev.off()
cat("Saved figures/scratch_peak_boundary_ratio_vs_R0.png\n")
cat(sprintf("R0_crit(r=%.3f) = %.6f (dashed blue line)\n", r_fixed, R0_crit))
cat("Red points (documented large-error cells): R0 =", paste(flagged, collapse=", "), "\n")

library(odin)
library(dde) 

set.seed(101)

T1 <- 1000
T2 <- 200
I2i <- 5

## model parameters 
N <- 10000
beta <- c(0.2, 0.19)
gamma <- c(0.1, 0.1)
mu <- 0.02
I1_ini <- 5

## compile model
twostrain_generator <- odin::odin("./odin_twostrain0.R")

## resident-only run to QE 
x1 <- twostrain_generator$new(
  beta = beta,
  gamma = gamma,
  mu = mu,
  N = N,
  I_ini = c(I1_ini, 0),
  S_ini = N - I1_ini
)

res1 <- x1$run(0:T1)

## Pick restart time 
idx <- nrow(res1)

S_qe  <- res1[idx, "S"]
I1_qe <- res1[idx, "I[1]"]
R_qe  <- res1[idx, "R"]

cat("S_qe =", S_qe, " I1_qe =", I1_qe, " R_qe =", R_qe, "\n")

## invasion 
S_ini2 <- S_qe
I_ini2 <- c(I1_qe, I2i)

x2 <- twostrain_generator$new(
  beta = beta,
  gamma = gamma,
  mu = mu,
  N = N,
  I_ini = I_ini2,
  S_ini = S_ini2
)

res2 <- x2$run(0:T2)

## early growth rate of invader 
I2 <- res2[, "I[2]"]
t  <- res2[, "step"]

keep <- (I2 > 0) & (I2 <= 20)  
if (sum(keep) >= 10) {
  fit <- lm(log(I2[keep]) ~ t[keep])
  slope <- coef(fit)[2]
  cat("Estimated early log-growth slope for invader (I2):", slope, "\n")
} else {
  slope <- NA_real_
  cat("Not enough points\n")
}

## Plot
pdf("odin_invasion_run.Rout.pdf")
par(mfrow = c(2, 1), mar = c(6, 4, 2, 1))

## Phase 1 plot (resident-only)
matplot(res1[, "step"], res1[, c("S", "I[1]", "R")],
        type = "l", lwd = 2, lty = 1,
        col = c("grey30", "red3", "black"),
        main = "resident-only",
        xlab = "t", ylab = "count", log = "y")

t_restart <- res1[idx, "step"]
abline(v = t_restart, col = "purple4", lwd = 2, lty = 3)

legend("right",
       legend = c("S", "I1", "R", "restart point"),
       col = c("grey30", "red3", "black", "purple4"),
       lty = c(1, 1, 1, 3),
       lwd = c(2, 2, 2, 2),
       bty = "n")

## Phase 2 plot (invasion)
matplot(res2[, "step"], res2[, c("I[1]", "I[2]")],
        type = "l", lwd = 2.5, lty = 1,
        col = c("red3", "blue3"),
        main = "invasion",
        xlab = "t", ylab = "infected", log = "y")
legend("right",
       legend = c("I1 (resident)", "I2 (invader)"),
       col = c("red3", "blue3"),
       lty = 1, lwd = 2.5, bty = "n")

mtext(sprintf("Invasion summary: I2i=%d, slope=%s", I2i,
              ifelse(is.na(slope), "NA", format(signif(slope, 4)))),
      side = 1, line = 4)

dev.off()

cat("Figure saved to: odin_invasion_run.Rout.pdf\n")
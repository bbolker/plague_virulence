library(burnout)

P1_at_DFE <- function(R0 = 2,
                      epsilon = 0.02,
                      c0 = 0.5,
                      S = 1e5, #when it's 1e6 P1 almost constant...
                      D = 1,
                      r = 0.5) {
  z = final_size(R0)
  N_DFE = c0 * S / (1 - (1 + r) * (1 - z * D) * (1 - c0))
  P1 = P1_prob(R0, epsilon = epsilon, k = 1, N = N_DFE)
  return(P1)
}

R0_vals <- c(1, 1.5, 2, 2.5)
epsilon_vals <- c(0.005, 0.01, 0.02, 0.05)
D_values <- seq(0, 1, length.out = 1000)

par(mfrow = c(4, 4), mar = c(3, 3, 2, 1), oma = c(2, 2, 3, 1))

for (R0 in R0_vals) {
  for (epsilon in epsilon_vals) {
    P1_values <- tryCatch({
      sapply(D_values, function(D) P1_at_DFE(R0 = R0, epsilon = epsilon, D = D))
    }, error = function(e) {
      rep(NA, length(D_values))
    })

    plot(D_values, P1_values, type = "l", col = "blue", lwd = 2,
         xlab = "D", ylab = "P1", main = paste0("R0=", R0, ", eps=", epsilon),
         ylim = c(0, 1))
  }
}

mtext("P1 at DFE vs D for different R0 and epsilon", outer = TRUE, cex = 1.2, line = 1)

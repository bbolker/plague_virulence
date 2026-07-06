# Load required packages
library(burnout)
source("simulation_funs.R")

# parameters
params <- list(
  R0 = 2.5,
  alpha = 5e-5,
  K = 1e6,
  r = 0.5,
  c = 0.2,
  epsilon = 0.02,
  initial_infected = 10
)

# replicate
set.seed(216)
results <- replicate(20, {
  cat(".")
  simulate_metapopulation(
    R0 = params$R0,
    alpha = params$alpha,
    K = params$K,
    r = params$r,
    c = params$c,
    epsilon = params$epsilon,
    initial_infected = params$initial_infected
  )
}, simplify = FALSE)

total_pops <- sapply(results, function(res) {
  apply(res$N[, , 1], 2, sum)
})

infected_patches <- sapply(results, function(res) {
  apply(res$I[, , 1], 2, sum)
})

total_inf <- sapply(results, function(res) {
  res$total_inf
})

# plot
par(mfrow = c(3, 1), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

# plot total population
matplot(total_pops, type = "l", lty = 1, 
        col = "grey",
        main = "Total Population Over Time",
        xlab = "Year", ylab = "Total Population")
lines(rowMeans(total_pops), col = "blue", lwd = 3)

# plot number of infected patches
matplot(infected_patches, type = "l", lty = 1, 
        col = "grey",
        main = "Number of Infected Patches",
        xlab = "Year", ylab = "Infected Patches")
lines(rowMeans(infected_patches), col = "red", lwd = 3)

# plot total infections
matplot(total_inf, type = "l", lty = 1, 
        col = "grey",
        main = "Total Infections Over Time",
        xlab = "Year", ylab = "Total Infections")
lines(rowMeans(total_inf), col = "purple", lwd = 3)
# show parameters
# show parameters
mtext(sprintf("Parameters: R0=%.1f, α=%.0e, K=%.0e, r=%.1f, c=%.1f, ε=%.2f, initial_infected=%d", 
              params$R0, params$alpha, params$K, params$r, params$c, params$epsilon, params$initial_infected),
      side = 3, line = 0, outer = TRUE, cex = 0.9)
par(mfrow = c(1, 1))

filename <- sprintf("figures/R0_%.1f_alpha_%.0e_K_%.0e_r_%.1f_c_%.1f_epsilon_%.2f_init_%d.pdf",
                    params$R0, params$alpha, params$K, params$r, params$c, params$epsilon, params$initial_infected)

# save
dev.copy(pdf, file = filename, width = 8, height = 10)
dev.off()

cat("Figure saved to:", filename, "\n")

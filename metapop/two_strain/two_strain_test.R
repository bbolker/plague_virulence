library(burnout) 
source("sim_fun.R")

## if profiling
## Rprof()

test_params <- list(
  n_patches = 100,
  n_years = 500,    
  K = 1e6,
  r = 0.5,
  c0 = 0.2,              
  nu = 5,                
  rho = 3,               
  alpha = 5e-6,          
  R01 = 2.0,             
  R02 = 1.8,            
  invade_year = 100,     
  initial_inf_ratio_1 = 0.1,
  initial_inf_ratio_2 = 0.05,
  coinf_approx = "polyfit"
)

nsim <- 50

all_total_inf1 <- matrix(0, nrow = test_params$n_years, ncol = nsim)
all_total_inf2 <- matrix(0, nrow = test_params$n_years, ncol = nsim)
all_patches1 <- matrix(0, nrow = test_params$n_years, ncol = nsim)
all_patches2 <- matrix(0, nrow = test_params$n_years, ncol = nsim)

cat(sprintf("Running %d simulations...\n", nsim))

for (i in 1:nsim) {
  set.seed(123 + i) 
  result <- do.call(simulate_metapopulation_2strain, test_params)
  all_total_inf1[, i] <- result$total_inf1
  all_total_inf2[, i] <- result$total_inf2
  all_patches1[, i] <- colSums(result$I1)
  all_patches2[, i] <- colSums(result$I2)
}

extinct_years1 <- apply(all_patches1, 2, function(col) {
  idx <- which(col == 0)
  if (length(idx) > 0) return(idx[1]) else return(NA_integer_)
})
persist_rate1 <- sum(is.na(extinct_years1)) / nsim

extinct_years2 <- apply(all_patches2, 2, function(col) {
  idx <- which(col == 0 & seq_along(col) > test_params$invade_year)
  if (length(idx) > 0) return(idx[1]) else return(NA_integer_)
})
persist_rate2 <- sum(is.na(extinct_years2)) / nsim

cat(sprintf("Resident Persistence Rate: %.1f%%\n", persist_rate1 * 100))
cat(sprintf("Invader Persistence Rate:  %.1f%%\n", persist_rate2 * 100))

inf1_cond <- all_total_inf1
inf1_cond[all_patches1 == 0] <- NA 
avg_total_inf1 <- rowMeans(inf1_cond, na.rm = TRUE)

patches1_cond <- all_patches1
patches1_cond[all_patches1 == 0] <- NA
avg_patches1 <- rowMeans(patches1_cond, na.rm = TRUE)

inf2_cond <- all_total_inf2
inf2_cond[all_patches2 == 0] <- NA 
avg_total_inf2 <- rowMeans(inf2_cond, na.rm = TRUE)

patches2_cond <- all_patches2
patches2_cond[all_patches2 == 0] <- NA
avg_patches2 <- rowMeans(patches2_cond, na.rm = TRUE)

pdf("two_strain_test.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))

leg_text1 <- sprintf("Resident (%.0f%% Survive)", persist_rate1 * 100)
leg_text2 <- sprintf("Invader (%.0f%% Survive)", persist_rate2 * 100)

y_max_inf <- max(c(all_total_inf1, all_total_inf2), na.rm = TRUE) * 1.05

plot(NULL, xlim = c(1, test_params$n_years), ylim = c(0, y_max_inf),
     xlab = "Year", ylab = "Total Infections", main = "Total Infections")

for (i in 1:nsim) {
  lines(all_total_inf1[, i], col = rgb(0.7, 0.7, 0.7, 0.3)) 
  lines(all_total_inf2[, i], col = rgb(0.7, 0.7, 0.7, 0.3)) 
}

lines(avg_total_inf1, col = "blue", lwd = 2.5)
lines(avg_total_inf2, col = "red", lwd = 2.5)
abline(v = test_params$invade_year, col = "black", lty = 2)

legend("topright", 
       legend = c(leg_text1, leg_text2), 
       col = c("blue", "red"), 
       lwd = c(2.5, 2.5), 
       lty = c(1, 1),
       cex = 0.8, bg = "white")

plot(NULL, xlim = c(1, test_params$n_years), ylim = c(0, test_params$n_patches),
     xlab = "Year", ylab = "Infected Patches", main = "Infected Patches")

for (i in 1:nsim) {
  lines(all_patches1[, i], col = rgb(0.7, 0.7, 0.7, 0.3)) 
  lines(all_patches2[, i], col = rgb(0.7, 0.7, 0.7, 0.3)) 
}

lines(avg_patches1, col = "blue", lwd = 2.5)
lines(avg_patches2, col = "red", lwd = 2.5)
abline(v = test_params$invade_year, col = "black", lty = 2)

legend("topright", 
       legend = c(leg_text1, leg_text2), 
       col = c("blue", "red"), 
       lwd = c(2.5, 2.5), 
       lty = c(1, 1),
       cex = 0.8, bg = "white")

param_text <- paste0(
  "R01=", test_params$R01, ", R02=", test_params$R02,
  ", alpha=", test_params$alpha, ", rho=", test_params$rho,
  ", c0=", test_params$c0, ", r=", test_params$r,
  ", invade_year=", test_params$invade_year,
  ", init_inf1=", test_params$initial_inf_ratio_1,
  ", init_inf2=", test_params$initial_inf_ratio_2
)
mtext(param_text, side = 3, outer = TRUE, cex = 0.8, line = 0)

par(mfrow = c(1, 1))
dev.off()

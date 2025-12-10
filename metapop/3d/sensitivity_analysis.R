getDFE <- function(c0 = 0.01,
                   K = 1e6,
                   alpha = 5e-6,
                   R0 = 2.5,
                   epsilon = 0.05,
                   D = 1,
                   r = 0.5,
                   n = 100
) {
  if (!require("burnout")) stop(
    "please install the 'burnout' package: ",
    "`remotes::install_github('davidearn/burnout')`")
  
  if (!require("rootSolve")) stop(
    "Please install the 'rootSolve' package:\n",
    "install.packages('rootSolve')")
  
  f <- function(x) {
    ni <- x[1]
    ns <- x[2]
    
    p <- 0
    z <- final_size(R0)
    
    ## burnout
    ni1 <- ni * (1 - z * D)
    
    ## logistic growth
    ni2 <- ni1 + r * ni1 * (1 - ni1 / K)
    ns1 <- ns + r * ns * (1 - ns / K)
    
    ## host movement
    ni3 <- ni2 + c0 * (ns1 - ni2)
    
    ## infection between patches (with p=0)
    B <- alpha * c0 * n * p * z * ni  # This equals 0 when p=0
    ni4 <- (ni3 + B * ns1) / (1 + B)  # Simplifies to ni3 when B=0
    
    eq1 <- ni4 - ni
    eq2 <- ns1 - ns
    
    return(c(eq1, eq2))
  }
  
  x0 <- c(K * 0.5, K)  
  result <- multiroot(f, x0)
  
  return(list(ni = result$root[1], ns = result$root[2], p = 0, 
              convergence = result$estim.precis))
}

getEE <- function(c0 = 0.01,
                  K = 1e6,
                  alpha = 5e-6,
                  R0 = 2.5,
                  epsilon = 0.05,
                  D = 1,
                  r = 0.5,
                  n = 100
) {
  if (!require("burnout")) stop(
    "please install the 'burnout' package: ",
    "`remotes::install_github('davidearn/burnout')`")
  
  if (!require("nleqslv")) stop(
    "Please install the 'nleqslv' package:\n",
    "install.packages('nleqslv')")
  
  f <- function(x) {
    p <- x[1]
    ni <- x[2]
    ns <- x[3]
    
    ## burnout
    Pb <- burnout_prob(R0, epsilon, N = ni)
    z <- final_size(R0)
    ni1 <- ni * (1 - z * D)
    ns1 <- (ni * Pb * p + ns * (1 - p)) / (Pb * p + 1 - p)
    p1 <- (1 - Pb) * p
    
    ## logistic growth
    ns2 <- ns1 + r * ns1 * (1 - ns1 / K)
    ni2 <- ni1 + r * ni1 * (1 - ni1 / K)
    p2 <- p1
    
    ## host movement
    delta_N <- c0 * (ns2 - ni2)
    ni3 <- ni2 + delta_N * (1 - p2)
    ns3 <- ns2 - delta_N * p2
    p3 <- p2
    
    ## infection between patches
    delta_p <- (1 - exp(-alpha * c0 * n * p3 * z * ni)) * (1 - p3)
    ni4 <- (p3 * ni3 + delta_p * ns3) / (p3 + delta_p)
    p4 <- p3 + delta_p
    
    ## fizzle
    Pf <- 1 / R0
    ns5 <- ((1 - p4) * ns3 + Pf * p4 * ni4) / (1 - p4 + Pf * p4)
    p5 <- (1 - Pf) * p4
    ni5 <- ni4
    
    eq1 <- p5 - p
    eq2 <- ni5 - ni
    eq3 <- ns5 - ns
    
    return(c(eq1, eq2, eq3))
  }
  
  x0 <- c(0.5, K * 0.5, K * 0.5)  
  result <- multiroot(f, x0,positive = TRUE)

  return(list(p = result$root[1], ni = result$root[2], ns = result$root[3],
              convergence = result$estim.precis))
}

# Plot EE p values vs R0
valid_data <- ! is.na(results_R0$p_EE)
R0_valid <- results_R0$R0[valid_data]
p_EE_valid <- results_R0$p_EE[valid_data]

plot(R0_valid, p_EE_valid, 
     xlab = "R0", ylab = "p", 
     main = "p vs R0",
     pch = 19, col = "blue")
lines(smooth.spline(R0_valid, p_EE_valid), col = "blue", lwd = 2)
grid()

# R0 sensitivity analysis 
R0_values <- seq(1.0, 3.0, by = 0.1)
alpha_values <- 5* 10^seq(-6, -2, by = 0.5)

# color palette for different alpha values
colors <- rainbow(length(alpha_values))

# Initialize plot
plot(NULL, xlim = range(R0_values), ylim = c(0, 1),
     xlab = "R0", ylab = "p (proportion of infected patches)",
     main = "Equilibrium p vs R0",
     las = 1)
grid()

# Store results for legend 
legend_labels <- list()
legend_colors <- list()

# Loop through each alpha value
for (j in seq_along(alpha_values)) {
  alpha_val <- alpha_values[j]
  
  # Initialize results dataframe for this alpha
  results_R0 <- data.frame(R0 = R0_values, p_EE = NA, NI_EE = NA, NS_EE = NA)
  
  cat("Processing alpha =", alpha_val, "\n")
  
  # Loop through R0 values
  for (i in seq_along(R0_values)) {
    tryCatch({
      ee_temp <- getEE(R0 = R0_values[i], alpha = alpha_val)
      if (ee_temp$convergence < 1e-2) { 
        results_R0[i, c("p_EE", "NI_EE", "NS_EE")] <- c(ee_temp$p, ee_temp$ni, ee_temp$ns)
      }
    }, error = function(e) cat("EE failed for R0 =", R0_values[i], "\n"))
  }
  
  # Plot valid data points 
  valid_data <- ! is.na(results_R0$p_EE)
  if (sum(valid_data) > 0) {
    R0_valid <- results_R0$R0[valid_data]
    p_EE_valid <- results_R0$p_EE[valid_data]
    
    # Add points
    points(R0_valid, p_EE_valid, pch = 19, col = colors[j], cex = 0.8)
    
    # Add smooth line 
    if (length(R0_valid) > 3) {
      lines(smooth.spline(R0_valid, p_EE_valid), col = colors[j], lwd = 2)
    } else {
      lines(R0_valid, p_EE_valid, col = colors[j], lwd = 2)
    }
    
    # add legend
    legend_labels[[length(legend_labels) + 1]] <- sprintf("alpha = %.0e", alpha_val)
    legend_colors[[length(legend_colors) + 1]] <- colors[j]
  }
}

# Add legend 
if (length(legend_labels) > 0) {
  legend("topright", 
         legend = unlist(legend_labels), 
         col = unlist(legend_colors), 
         lwd = 2, pch = 19, cex = 0.8, bg = "white")
}
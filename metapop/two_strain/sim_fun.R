#' Two-Strain Metapopulation Simulation with Vectorized Logic
#' @param n_patches Number of patches in the metapopulation
#' @param n_years Number of years to simulate
#' @param K Carrying capacity per patch
#' @param r Growth rate of the host population
#' @param c0 Migration probability (proportion of hosts moving)
#' @param nu Scaling constant (dead hosts -> total environmental fleas)
#' @param rho Carrying capacity of fleas per migrating host (The Bottleneck)
#' @param alpha Transmission efficiency per flea (infectivity)
#' @param D Disease-induced mortality (usually 1 for plague)
#' @param R01 Basic reproduction number for Strain 1 (Resident)
#' @param R02 Basic reproduction number for Strain 2 (Invader)
#' @param invade_year Year when Strain 2 is introduced
#' @param initial_inf_ratio_1 Proportion of patches initially infected with Strain 1
#' @param initial_inf_ratio_2 Proportion of patches to seed Strain 2 during invasion
#' @param initial_pop_ratio Initial population as ratio of K
#' @param early_stop If TRUE, stops simulation when both strains are extinct
#' @return A list containing N (Host dynamics), I (Indicators), S (Mortality), and Total Infection counts
simulate_metapopulation_2strain <- function(
    n_patches = 100,
    n_years = 500,
    K = 1e6,
    r = 0.5,
    c0 = 0.2,
    nu = 5,        
    rho = 3,       
    alpha = 5e-6,   
    D = 1,
    R01 = 2.0,     
    R02 = 1.9,    
    invade_year = 100,      
    initial_inf_ratio_1 = 0.1, 
    initial_inf_ratio_2 = 0.02, 
    initial_pop_ratio = 1,
    early_stop = FALSE      
) {
  if (!require("burnout", quietly = TRUE)) stop("Please install the 'burnout' package")
  
  # Vectorized helper for the Final Size Equation (1 - z = exp(-R0 * z))
  z_func_vec <- function(R0_vec) {
    res <- numeric(length(R0_vec))
    idx <- R0_vec > 1
    if(any(idx)) {
      # Apply final_size to each valid R0 entry
      res[idx] <- sapply(R0_vec[idx], final_size)
    }
    return(res)
  }
  
  # Initialize recording structures
  N <- array(0, dim = c(n_patches, n_years, 4)) # Total host population across 4 stages
  I1 <- matrix(0, nrow = n_patches, ncol = n_years) # Indicator for Strain 1
  I2 <- matrix(0, nrow = n_patches, ncol = n_years) # Indicator for Strain 2
  S1 <- matrix(0, nrow = n_patches, ncol = n_years) # Deaths from Strain 1
  S2 <- matrix(0, nrow = n_patches, ncol = n_years) # Deaths from Strain 2
  total_inf1 <- rep(0, n_years) # Global infection count Strain 1
  total_inf2 <- rep(0, n_years) # Global infection count Strain 2
  
  # Initial Setup: Seed Resident Strain (Strain 1)
  N[, 1, 1] <- initial_pop_ratio * K
  infected_init_1 <- sample(1:n_patches, initial_inf_ratio_1 * n_patches)
  z_init <- final_size(R01)
  S1[infected_init_1, 1] <- z_init * D * N[infected_init_1, 1, 1]
  I1[infected_init_1, 1] <- 1 
  total_inf1[1] <- sum(S1[, 1])
  
  # Arithmetic mean of R0 for co-infection scenarios
  R0_mean <- (R01 + R02) / 2
  
  for (k in 2:n_years) {
    
    # --- Stage 1 -> 2: Host Growth ---
    # Retrieve previous year-end population (Stage 4)
    N_prev_end <- if(k == 2) N[, 1, 1] - S1[, 1] else N[, k-1, 4]
    
    # Logistic growth formula
    N_after_growth <- N_prev_end + r * N_prev_end * (1 - N_prev_end / K)
    N[, k, 2] <- N_after_growth
    
    # --- Stage 2 -> 3: Host Migration (Colonization) ---
    pop_total_network <- sum(N_after_growth)
    N_after_col <- (1 - c0) * N_after_growth + c0 * (pop_total_network / n_patches)
    N[, k, 3] <- N_after_col
    
    # --- Transmission Dynamics: Force of Infection ---
    S_last_total <- S1[, k-1] + S2[, k-1]
    
    # Proportions of strains in the flea pool (q-vectors)
    q1 <- ifelse(S_last_total > 0, S1[, k-1] / S_last_total, 0)
    q2 <- ifelse(S_last_total > 0, S2[, k-1] / S_last_total, 0)
    
    # Calculate Transport Bottleneck (Min function for fleas vs migration capacity)
    E_total <- pmin(nu * S_last_total, rho * c0 * N_after_growth)
    R_total <- nu * S_last_total - E_total # Retained fleas
    
    # Force of Infection 
    lambda1 <- alpha * (q1 * R_total + sum(q1 * E_total) / n_patches) * max(0, 1 - 1/R01)
    lambda2 <- alpha * (q2 * R_total + sum(q2 * E_total) / n_patches) * max(0, 1 - 1/R02)
    
    # Poisson sampling for initial infection seeds (a and b)
    a <- rpois(n_patches, lambda1)
    b <- rpois(n_patches, lambda2)
    
    # --- Invasion Event: Introduction of Strain 2 ---
    if (k == invade_year) {
      invade_idx <- sample(1:n_patches, initial_inf_ratio_2 * n_patches)
      b[invade_idx] <- b[invade_idx] + 1 # Force establishment in selected patches
    }
    # Zero out Strain 2 if it's before the invasion year
    if (k < invade_year) b <- rep(0, n_patches)
    
    # --- Stage 3 -> 4: Epidemic Outbreak & Resource Partitioning ---
    ind1 <- a > 0 # Indicator for successful establishment
    ind2 <- b > 0
    ind_any <- (a + b) > 0
    
    I1[, k] <- as.numeric(ind1)
    I2[, k] <- as.numeric(ind2)
    
    # Assign Effective R0 based on establishment indicators
    R_eff <- R01 * (ind1 & !ind2) + R02 * (!ind1 & ind2) + R0_mean * (ind1 & ind2)
    
    # Calculate total death toll Z using vectorized final size
    Z_total <- ind_any * z_func_vec(R_eff) * D * N_after_col
    
    # Partition Z based on initial seed ratios (Galton-Watson allocation)
    # denom is (a + b); we only divide when ind_any is TRUE to avoid 0/0
    denom <- a + b
    S1[, k] <- ifelse(ind_any, Z_total * a / denom, 0)
    S2[, k] <- ifelse(ind_any, Z_total * b / denom, 0)
    
    # Update end-of-year surviving population
    N[, k, 4] <- N_after_col - (S1[, k] + S2[, k])
    
    # Record global time-series statistics
    total_inf1[k] <- sum(S1[, k])
    total_inf2[k] <- sum(S2[, k])
    
    # Survival check: exit if both strains disappear after invasion starts
    if (early_stop && sum(ind_any) == 0 && k > invade_year) break
  }
  
  return(list(
    N = N, 
    I1 = I1, I2 = I2, 
    S1 = S1, S2 = S2, 
    total_inf1 = total_inf1, 
    total_inf2 = total_inf2
  ))
}